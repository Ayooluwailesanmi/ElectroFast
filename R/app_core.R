# Calculation layer for the Shiny app: works on uploaded files (temporary path
# plus original name), so it does not depend on the project's folder layout.
# Reports only the metrics common in electroanalytical journals: slope (mean +/-
# SD over days), R2, linear range, % of Nernstian, selectivity coefficients.

source("R/read_emf16_grf.R")
source("R/activity.R")

nernst_S = function(T_C = 25, z = 1) log(10) * 8.314462618 * (T_C + 273.15) / (96485.33212 * abs(z)) * 1000

# File name -> day and standard letter. Handles "7-4-2026_8.9.34 A.txt", doubled
# letters ("AA"), sodium-series names ("AN"), a repeat marked "C MAIN", and a
# missing space ("17.41.47C.txt").
parse_emf_name = function(name) {
  stem = sub("\\.txt$", "", name, ignore.case = TRUE)
  main = grepl("MAIN$", stem)
  stem = sub("\\s*MAIN$", "", stem)
  # A letter code is read only from EMF16-style names ("9-29-2026_17.13.17 A"), so an
  # arbitrary name such as "LiCl 1mM" is not mistaken for code "M"
  emf_style = grepl("^[0-9]{1,2}-[0-9]{1,2}-[0-9]{4}_[0-9]{1,2}\\.[0-9]{1,2}\\.[0-9]{1,2}\\s*[A-Za-z]+$", stem)
  code = if (emf_style) toupper(regmatches(stem, regexpr("[A-Za-z]+$", stem))) else character()
  if (!length(code)) code = NA_character_
  day = regmatches(name, regexpr("^[0-9]+-[0-9]+-[0-9]{4}", name))
  if (!length(day)) day = "run"
  data.frame(file = name, day = day, code = code, letter = substr(code, 1, 1), main = main)
}

# Standard codes and their default concentrations. "july": the July-August runs,
# A to I = 1e-9 to 1e-1 M in decade steps. "protocol": the five-day protocol
# (2026-09-28), A to O with mid-decade points up to 0.5 M, and X, Y, Z = the 1e-3 M
# check standard read at the start, middle and end of each day.
code_schemes = list(
  july = setNames(10^(-9:-1), LETTERS[1:9]),
  protocol = c(A = 1e-9, B = 1e-8, C = 1e-7, D = 1e-6, E = 1e-5, F = 3e-5, G = 1e-4, H = 3e-4,
               I = 1e-3, J = 3e-3, K = 1e-2, L = 3e-2, M = 1e-1, N = 3e-1, O = 5e-1,
               X = 1e-3, Y = 1e-3, Z = 1e-3)
)

# Concentration written in a file name, in mol/L: "1 mM", "0.1M", "5 uM" or "5 µM", "1e-3 M" or "1e-3", "10-3 M" or
# "10^-3 M". A unit is needed for plain numbers, so a date or a time in the name is not read as a concentration.
# NA when the name holds no concentration or two different ones.
conc_from_name = function(name) {
  s = sub("\\.[A-Za-z0-9]+$", "", name)
  unit = c(m = 1, mm = 1e-3, um = 1e-6, "\u00b5m" = 1e-6, nm = 1e-9, pm = 1e-12)
  v = numeric(0)
  m1 = regmatches(s, gregexpr("(?<![0-9.eE+-])[0-9]+(\\.[0-9]+)?\\s*(mM|uM|\u00b5M|nM|pM|M)(?![A-Za-z])", s, perl = TRUE))[[1]]
  for (x in m1) { u = tolower(sub("^[0-9.]+\\s*", "", x)); v = c(v, as.numeric(sub("\\s*[A-Za-z\u00b5]+$", "", x)) * unit[[u]]) }
  m2 = regmatches(s, gregexpr("(?<![0-9.])[0-9.]+[eE]-[0-9]+(\\s*M(?![A-Za-z]))?", s, perl = TRUE))[[1]]
  v = c(v, as.numeric(sub("\\s*M$", "", m2)))
  m3 = regmatches(s, gregexpr("(?<![0-9])10\\^?-([0-9]{1,2})\\s*M(?![A-Za-z])", s, perl = TRUE))[[1]]
  v = c(v, 10^-as.numeric(sub("^10\\^?-([0-9]{1,2}).*$", "\\1", m3)))
  v = unique(signif(v, 10))
  if (length(v) == 1) v else NA_real_
}

# Where concentrations come from when there is no sample sheet. EMF16-style names with letter codes keep their
# coding (letters J to O, or the X, Y, Z check standard, occur only in the five-day protocol); otherwise the
# concentration written in each name; otherwise a decade series in recording order, as a suggestion to check.
detect_scheme = function(names) {
  L = do.call(rbind, lapply(names, parse_emf_name))$letter
  if (any(!is.na(L))) return(if (any(L %in% c(LETTERS[10:15], "X", "Y", "Z"))) "protocol" else "july")
  if (all(!is.na(vapply(names, conc_from_name, 0)))) "names" else "order"
}
scheme_labels = c(names = "Concentration written in the file name", order = "Decade series in recording order (suggested)",
                  july = "Letter codes A to I", protocol = "Letter codes A to O with X, Y, Z check standards")

# Start time recorded inside an EMF16 file ("Started: 9/29/2026, 5:13:16 PM"), for the recording order
emf_start_time = function(path) {
  x = tryCatch(readLines(path, n = 40, warn = FALSE), error = function(e) character())
  s = grep("^\\s*Started:", x, value = TRUE)
  if (!length(s)) return(NA_real_)
  as.numeric(as.POSIXct(trimws(sub("^\\s*Started:\\s*", "", s[1])), format = "%m/%d/%Y, %I:%M:%S %p", tz = "UTC"))
}

# Day recorded inside an EMF16 file ("Started: 9/29/2026, 5:13:16 PM"), written as
# m-d-yyyy like the file-name date, so the day does not depend on the file name
emf_start_day = function(path) {
  x = tryCatch(readLines(path, warn = FALSE), error = function(e) character())
  s = grep("^\\s*Started:", x, value = TRUE)
  d = if (length(s)) regmatches(s[1], regexpr("[0-9]{1,2}/[0-9]{1,2}/[0-9]{4}", s[1])) else character()
  if (length(d)) gsub("/", "-", d) else NA_character_
}

truthy = function(v) toupper(trimws(as.character(v))) %in% c("TRUE", "T", "1", "YES", "Y")

# Optional sample sheet (CSV): one row per file, with at least the file name and its
# concentration, so files can be named anything. Optional columns: day, include, fit.
sample_sheet_template = function() {
  data.frame(file = c("9-29-2026_17.13.17 A.txt", "my standard 1mM.txt"), conc_M = c(1e-9, 1e-3),
             day = c("", "10-01-2026"), include = c("", ""), fit = c("", ""))
}
read_sample_sheet = function(path) {
  s = utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE, colClasses = "character")
  names(s) = gsub("[^a-z_]", "", tolower(names(s)))
  names(s)[names(s) %in% c("conc", "concentration", "concm", "concentrationm")] = "conc_m"
  if (!all(c("file", "conc_m") %in% names(s)))
    stop("The sample sheet needs the columns 'file' and 'conc_M' (download the template to see the layout).")
  s
}

# Which standards enter the fit: by letter when the file has a code, otherwise by
# concentration (default 1e-5 to 1e-1 M, the range used so far); a 'fit' column in
# the sample sheet overrides both.
fit_flags = function(st, fit_letters, conc_range = c(1e-5, 1e-1)) {
  f = ifelse(!is.na(st$letter), st$letter %in% fit_letters,
             !is.na(st$conc_M) & st$conc_M >= conc_range[1] * 0.999 & st$conc_M <= conc_range[2] * 1.001)
  if (!is.null(st$fit_sheet)) f = ifelse(is.na(st$fit_sheet), f, st$fit_sheet)
  st$include & f
}

# Standards table for uploaded files. Concentration comes from the sample sheet when
# one is given, otherwise from the letter code; the day from the sheet, otherwise from
# inside the file, otherwise from the file name. Anything still unknown is listed with
# a note, never silently dropped. A repeat marked MAIN replaces the earlier reading.
standards_table = function(files, fit_letters = c("E", "F", "G", "H", "I"), scheme = "july", sheet = NULL) {
  st = do.call(rbind, lapply(files$name, parse_emf_name))
  st$datapath = files$datapath
  inside = vapply(as.character(files$datapath), function(p) if (is.na(p)) NA_character_ else emf_start_day(p), "")
  st$day = ifelse(!is.na(inside), inside, st$day)
  suggested = rep(FALSE, nrow(st))
  if (scheme %in% names(code_schemes)) st$conc_M = unname(code_schemes[[scheme]][st$letter]) else {
    # no letter coding: the concentration written in each name, then (route "order") a decade series from 1e-9 M
    # upward in recording order within each day, marked as a suggestion
    st$letter = NA_character_
    st$conc_M = vapply(st$file, conc_from_name, 0, USE.NAMES = FALSE)
    if (scheme == "order") {
      tt = vapply(as.character(files$datapath), function(p) if (is.na(p)) NA_real_ else emf_start_time(p), 0)
      ord = ifelse(is.na(tt), seq_len(nrow(st)), tt)
      for (d in unique(st$day)) { k = which(st$day == d & is.na(st$conc_M)); k = k[order(ord[k])]
        n = min(length(k), 9); if (n) { st$conc_M[k[1:n]] = 10^(-9 + seq_len(n) - 1); suggested[k[1:n]] = TRUE } }
    }
  }
  st$fit_sheet = NA
  inc_sheet = rep(NA, nrow(st))
  in_sheet = rep(NA, nrow(st))
  if (!is.null(sheet)) {
    key = function(x) tolower(sub("\\.[A-Za-z0-9]+$", "", basename(trimws(x))))
    m = match(key(st$file), key(sheet$file))
    in_sheet = !is.na(m)
    conc = suppressWarnings(as.numeric(sheet$conc_m))[m]
    st$conc_M = ifelse(in_sheet & !is.na(conc), conc, st$conc_M)
    if ("day" %in% names(sheet)) { d = sheet$day[m]; use = in_sheet & !is.na(d) & nzchar(d); st$day[use] = d[use] }
    if ("include" %in% names(sheet)) { v = sheet$include[m]; use = in_sheet & !is.na(v) & nzchar(v); inc_sheet[use] = truthy(v[use]) }
    if ("fit" %in% names(sheet)) { v = sheet$fit[m]; use = in_sheet & !is.na(v) & nzchar(v); st$fit_sheet[use] = truthy(v[use]) }
  }
  # a file with no date anywhere joins the day of the other files when they share one;
  # otherwise it is flagged for the user to type the day
  no_day = st$day == "run"
  known = unique(st$day[!no_day])
  if (any(no_day) && length(known) == 1) st$day[no_day] = known
  key = paste(st$day, ifelse(is.na(st$letter), st$file, st$letter))
  dup = duplicated(key) | duplicated(key, fromLast = TRUE)
  st$include = !is.na(st$conc_M) & !(dup & !st$main & ave(st$main, key, FUN = any))
  st$include = ifelse(is.na(inc_sheet), st$include, inc_sheet & !is.na(st$conc_M))
  st$fit = fit_flags(st, fit_letters)
  st$note = ifelse(dup & !st$include, "repeat replaced by MAIN", ifelse(dup & st$main, "MAIN reading used", ifelse(dup, "duplicate standard: check", "")))
  st$note = ifelse(!is.na(in_sheet) & !in_sheet, "not in the sample sheet", st$note)
  st$note = ifelse(no_day & length(known) == 1, "no date in the file: day taken from the other files", st$note)
  st$note = ifelse(st$day == "run", "day unknown: type it in the day cell", st$note)
  st$note = ifelse(suggested & !(!is.na(in_sheet) & in_sheet), "suggested from the recording order: check", st$note)
  st$note = ifelse(is.na(st$conc_M), "concentration unknown: type it in the conc_M cell", st$note)
  st[order(st$day, st$conc_M), c("file", "day", "letter", "conc_M", "include", "fit", "fit_sheet", "note", "datapath")]
}

# Parsed files are cached by path, so changing the reading time does not reread them
.emf_cache = new.env()
# Files that are not EMF16 exports are read as plain text: the first two numbers on
# each row are taken as time (s) and potential; potentials below 2 in absolute value
# are taken as volts and converted to mV.
read_trace_any = function(path) {
  emf = tryCatch(read_emf16_grf(path)$long, error = function(e) NULL)
  if (!is.null(emf)) return(emf)
  x = sub("\r$", "", readLines(path, warn = FALSE))
  num = lapply(strsplit(trimws(x), "[,;\t ]+"), function(t) suppressWarnings(as.numeric(t)))
  ok = vapply(num, function(v) length(v) >= 2 && !anyNA(v[1:2]), TRUE)
  if (sum(ok) < 5) stop("Could not read ", basename(path), ": it is neither an EMF16 file nor a text file ",
                        "with two numeric columns (time in s, potential).")
  m = do.call(rbind, lapply(num[ok], `[`, 1:2))
  mv = if (max(abs(m[, 2])) < 2) m[, 2] * 1000 else m[, 2]
  data.frame(time_s = m[, 1], channel = "1", emf_mV = mv)
}
emf_trace = function(path) {
  if (is.null(.emf_cache[[path]])) assign(path, read_trace_any(path), envir = .emf_cache)
  .emf_cache[[path]]
}

emf_at = function(path, read_at_s) {
  z = emf_trace(path)
  hit = abs(z$time_s - read_at_s) < 0.05
  if (!any(hit)) hit = which.min(abs(z$time_s - read_at_s))
  median(z$emf_mV[hit])
}

add_emf = function(st, read_at_s, background_I = 0, z = 1) {
  st = st[st$include, ]
  st$E_mV = vapply(st$datapath, emf_at, 0, read_at_s = read_at_s)
  g = davies_gamma(st$conc_M + background_I, z)
  st$neg_log_c = -log10(st$conc_M)
  st$neg_log_a = -log10(g * st$conc_M)
  st
}

# Per-day calibration over the standards marked fit, plus a summary over days
calibrate = function(pts, x = "neg_log_a", T_C = 25, z = 1) {
  S = nernst_S(T_C, z)
  per_day = do.call(rbind, lapply(split(pts, pts$day), function(d) {
    f = d[d$fit, ]
    if (nrow(f) < 2) return(NULL)
    m = lm(f$E_mV ~ f[[x]])
    data.frame(day = d$day[1], n = nrow(f), from_M = min(f$conc_M), to_M = max(f$conc_M),
               slope_mV_dec = unname(coef(m)[2]), intercept_mV = unname(coef(m)[1]),
               r2 = summary(m)$r.squared, pct_nernst = 100 * abs(coef(m)[2]) / S)
  }))
  if (is.null(per_day)) return(NULL)
  summary = data.frame(days = nrow(per_day), slope_mean = mean(per_day$slope_mV_dec),
                       slope_sd = if (nrow(per_day) > 1) sd(per_day$slope_mV_dec) else NA,
                       pct_nernst = 100 * abs(mean(per_day$slope_mV_dec)) / S,
                       r2_min = min(per_day$r2), r2_max = max(per_day$r2),
                       range = sprintf("%g to %g M", min(per_day$from_M), max(per_day$to_M)))
  list(per_day = per_day, summary = summary, S = S)
}

# Change between neighbouring standards against the Nernstian expectation
steps_table = function(pts, x = "neg_log_a", T_C = 25, z = 1) {
  S = nernst_S(T_C, z)
  pts = pts[!pts$letter %in% c("X", "Y", "Z"), ]   # check-standard repeats are not calibration steps
  do.call(rbind, lapply(split(pts, pts$day), function(d) {
    d = d[order(d$conc_M), ]
    if (nrow(d) < 2) return(NULL)
    dx = diff(d[[x]])
    f = d[d$fit, ]
    sgn = if (nrow(f) >= 2) sign(coef(lm(f$E_mV ~ f[[x]]))[2]) else 1
    data.frame(day = d$day[1], from = d$letter[-nrow(d)], to = d$letter[-1],
               decades = -dx, measured_mV = diff(d$E_mV), expected_mV = sgn * S * dx,
               pct_nernst = 100 * diff(d$E_mV) / (sgn * S * dx))
  }))
}

# Fixed interference method (IUPAC; Umezawa et al. 2000, eq. 2)
fim = function(pts, flat_letters, nernst_letters, interf_M, zA = 1, zB = 1, T_C = 25) {
  per_day = do.call(rbind, lapply(split(pts, pts$day), function(d) {
    fl = d$E_mV[d$letter %in% flat_letters]
    nf = d[d$letter %in% nernst_letters, ]
    if (length(fl) < 1 || nrow(nf) < 2) return(NULL)
    e_flat = mean(fl)
    m = lm(E_mV ~ neg_log_a, data = nf)
    aA = 10^(-(e_flat - coef(m)[1]) / coef(m)[2])
    aB = davies_gamma(interf_M + aA, zB) * interf_M
    K = unname(aA / aB^(zA / zB))
    data.frame(day = d$day[1], flat_mV = e_flat, nernst_slope = unname(coef(m)[2]),
               a_primary_cross = unname(aA), a_interferent = aB, K = K, lgK = log10(K))
  }))
  if (is.null(per_day)) return(NULL)
  list(per_day = per_day,
       summary = data.frame(days = nrow(per_day), K_mean = mean(per_day$K),
                            K_sd = if (nrow(per_day) > 1) sd(per_day$K) else NA,
                            lgK_mean = mean(per_day$lgK),
                            lgK_sd = if (nrow(per_day) > 1) sd(per_day$lgK) else NA))
}

# Separate solution method, aA = aB (IUPAC eq. 6). The primary-ion potentials are
# averaged over the calibration days at each concentration. The sign follows the
# primary calibration, so the result does not depend on the instrument's sign
# convention.
ssm = function(primary, interferent, concs, zA = 1, zB = 1, T_C = 25) {
  S = nernst_S(T_C, zA)
  f = primary[primary$fit, ]
  s_sign = sign(coef(lm(f$E_mV ~ I(-f$neg_log_a)))[2])   # slope sign against log a
  out = do.call(rbind, lapply(concs, function(cc) {
    eA = primary$E_mV[abs(log10(primary$conc_M) - log10(cc)) < 1e-6]
    if (!length(eA)) return(NULL)
    ib = interferent[abs(log10(interferent$conc_M) - log10(cc)) < 1e-6, ]
    if (!nrow(ib)) return(NULL)
    aA = davies_gamma(cc, zA) * cc
    do.call(rbind, lapply(seq_len(nrow(ib)), function(i) {
      lgk = s_sign * (ib$E_mV[i] - mean(eA)) / S + (1 - zA / zB) * log10(aA)
      data.frame(conc_M = cc, interferent_day = ib$day[i], E_primary_mean = mean(eA),
                 E_interferent = ib$E_mV[i], lgK = unname(lgk), K = unname(10^lgk))
    }))
  }))
  out
}

# Background-corrected Nernstian fit (Dagar et al. 2026, Sens. Actuators B 468,
# 140650, eq. 1): E = E0 + S log10(x + b_bkg), fitted to every included standard,
# so no linear range has to be chosen by eye. b_bkg is the effective background
# level that sets the lower limit of response. x is the activity or concentration,
# matching the chosen axis; the slope is reported against -log x, as the other fits.
fit_background = function(pts, x = "neg_log_a", T_C = 25, z = 1) {
  S = nernst_S(T_C, z)
  one_day = function(d) {
    v = 10^(-d[[x]])
    top = d[order(v, decreasing = TRUE)[seq_len(min(3, nrow(d)))], ]
    s0 = coef(lm(top$E_mV ~ I(-top[[x]])))[2]
    e0 = mean(top$E_mV - s0 * log10(10^(-top[[x]])))
    sse = function(p) sum((d$E_mV - (p[1] + p[2] * log10(v + 10^p[3])))^2)
    starts = lapply(-8:-3, function(lb) c(e0, s0, lb))
    best = NULL
    for (st in starts) {
      o = tryCatch(optim(st, sse, control = list(maxit = 5000, reltol = 1e-12)), error = function(e) NULL)
      if (!is.null(o) && (is.null(best) || o$value < best$value)) best = o
    }
    p = best$par
    r2 = 1 - best$value / sum((d$E_mV - mean(d$E_mV))^2)
    data.frame(day = d$day[1], n = nrow(d), E0_mV = p[1], slope_mV_dec = -p[2], b_bkg_M = 10^p[3],
               r2 = r2, pct_nernst = 100 * abs(p[2]) / S)
  }
  per_day = do.call(rbind, lapply(split(pts, pts$day), function(d) if (nrow(d) >= 4) one_day(d)))
  if (is.null(per_day)) return(NULL)
  summary = data.frame(days = nrow(per_day), slope_mean = mean(per_day$slope_mV_dec),
                       slope_sd = if (nrow(per_day) > 1) sd(per_day$slope_mV_dec) else NA,
                       pct_nernst = 100 * abs(mean(per_day$slope_mV_dec)) / S,
                       b_bkg_mean = mean(per_day$b_bkg_M), b_bkg_sd = if (nrow(per_day) > 1) sd(per_day$b_bkg_M) else NA,
                       r2_min = min(per_day$r2), r2_max = max(per_day$r2))
  list(per_day = per_day, summary = summary, S = S, x = x)
}

# Curve of the background-corrected fit for plotting
background_curve = function(bg, x_range, n = 200) {
  do.call(rbind, lapply(seq_len(nrow(bg$per_day)), function(i) {
    r = bg$per_day[i, ]
    xs = seq(x_range[1], x_range[2], length.out = n)
    data.frame(day = r$day, xval = xs, E_mV = r$E0_mV - r$slope_mV_dec * log10(10^(-xs) + r$b_bkg_M))
  }))
}

# Precision from repeated readings of one check standard, with day as a random
# effect: linear mixed model fitted by REML with lme4 (Carnoli et al. 2024, Anal.
# Chim. Acta 1304, 342444). Within-day SD = repeatability; within and between
# combined = intermediate precision (ICH Q2(R2) terms). SDs in mV are also given as
# the equivalent concentration error at the theoretical Nernstian slope.
precision_check = function(pts, conc, T_C = 25, z = 1, min_days = 5) {
  d = pts[abs(log10(pts$conc_M) - log10(conc)) < 1e-6, ]
  S = nernst_S(T_C, z)
  pct = function(s) 100 * (10^(s / S) - 1)
  per_day = table(d$day)
  notes = character()
  if (length(per_day) < 2) return(list(notes = "Readings of this standard are needed on at least 2 days."))
  if (length(per_day) < min_days)
    notes = c(notes, sprintf("Only %d days: below the rule-of-thumb minimum of %d levels for a random effect (Carnoli et al. 2024). Treat the between-day SD as indicative only.",
                             length(per_day), min_days))
  if (all(per_day < 2)) {
    notes = c(notes, "Only one reading per day: repeatability (within-day SD) cannot be separated. The SD across days is the intermediate precision only.")
    s_ip = sd(d$E_mV)
    metrics = data.frame(metric = c("Repeatability (within-day SD)", "Between-day SD", "Intermediate precision"),
                         SD_mV = c(NA, NA, s_ip), conc_error_pct = c(NA, NA, pct(s_ip)))
    return(list(metrics = metrics, blup = NULL, notes = notes, n_days = length(per_day), n_readings = nrow(d)))
  }
  m = lme4::lmer(E_mV ~ 1 + (1 | day), data = d, REML = TRUE)
  vc = as.data.frame(lme4::VarCorr(m))
  s_b = vc$sdcor[vc$grp == "day"]
  s_w = vc$sdcor[vc$grp == "Residual"]
  if (lme4::isSingular(m))
    notes = c(notes, "Singular fit: the between-day variance could not be separated from zero with these data.")
  s_ip = sqrt(s_b^2 + s_w^2)
  metrics = data.frame(metric = c("Repeatability (within-day SD)", "Between-day SD", "Intermediate precision"),
                       SD_mV = c(s_w, s_b, s_ip), conc_error_pct = pct(c(s_w, s_b, s_ip)))
  re = lme4::ranef(m, condVar = TRUE)$day
  se = sqrt(attr(re, "postVar")[1, 1, ])
  blup = data.frame(day = rownames(re), effect_mV = re[, 1], lower = re[, 1] - 1.96 * se,
                    upper = re[, 1] + 1.96 * se, readings = as.integer(per_day[rownames(re)]))
  list(metrics = metrics, blup = blup, grand_mean = unname(lme4::fixef(m)[1]), notes = notes,
       n_days = length(per_day), n_readings = nrow(d))
}

# Calibration across days with a common slope and a random day baseline (and a
# random day slope when the data support it)
precision_calibration = function(pts, x = "neg_log_a", min_days = 5) {
  f = pts[pts$fit, ]
  f$xv = f[[x]]
  n_days = length(unique(f$day))
  notes = character()
  if (n_days < 2) return(list(notes = "Calibrations on at least 2 days are needed."))
  if (n_days < min_days)
    notes = c(notes, sprintf("Only %d days: below the rule-of-thumb minimum of %d (Carnoli et al. 2024). Between-day estimates are indicative only.", n_days, min_days))
  m1 = lme4::lmer(E_mV ~ xv + (1 | day), data = f, REML = TRUE)
  vc1 = as.data.frame(lme4::VarCorr(m1))
  if (lme4::isSingular(m1)) notes = c(notes, "Singular fit for the day baseline: its between-day variance could not be separated from zero.")
  m2 = suppressWarnings(suppressMessages(tryCatch(lme4::lmer(E_mV ~ xv + (1 + xv | day), data = f, REML = TRUE), error = function(e) NULL)))
  slope_sd = if (!is.null(m2) && !lme4::isSingular(m2)) {
    v2 = as.data.frame(lme4::VarCorr(m2)); v2$sdcor[v2$grp == "day" & v2$var1 == "xv" & is.na(v2$var2)]
  } else NA
  if (is.na(slope_sd)) notes = c(notes, "Day-to-day slope variation could not be estimated from these data (random slope singular or not fitted).")
  list(table = data.frame(days = n_days, pooled_slope = unname(lme4::fixef(m1)[2]),
                          pooled_slope_se = sqrt(vcov(m1)[2, 2]),
                          between_day_baseline_sd = vc1$sdcor[vc1$grp == "day"],
                          between_day_slope_sd = slope_sd),
       notes = notes)
}

# Settling of each standard: the first time after which the drift over a sliding
# window stays below drift_mV_min (mV per minute) to the end of the recording.
# NA means the standard never settled within the file.
settling_table = function(pts, drift_mV_min = 1, window_s = 60, step_s = 5) {
  do.call(rbind, lapply(seq_len(nrow(pts)), function(i) {
    z = emf_trace(pts$datapath[i]); z = z[!duplicated(z$time_s), ]
    end = max(z$time_s); grid = seq(window_s, end, by = step_s)
    e1 = approx(z$time_s, z$emf_mV, xout = grid)$y; e0 = approx(z$time_s, z$emf_mV, xout = grid - window_s)$y
    bad = which(abs(e1 - e0) >= drift_mV_min * window_s / 60)
    settled = if (!length(bad)) window_s else if (max(bad) == length(grid)) NA else grid[max(bad) + 1]
    data.frame(day = pts$day[i], letter = pts$letter[i], conc_M = pts$conc_M[i], fit = pts$fit[i],
               length_s = end, settled_s = settled)
  }))
}

# Calibration slope per day when the EMF is read at each of several times. A day is
# left out at a time longer than any of its fitted recordings. spread = range of the
# slope over the times from 'from_s' on; flagged when it exceeds flag_mV.
slope_by_time = function(pts, times = c(30, 60, 120, 180, 300, 420, 600), x = "neg_log_a", T_C = 25, z = 1,
                         from_s = 120, flag_mV = 2) {
  len = vapply(pts$datapath, function(p) max(emf_trace(p)$time_s), 0)
  long = do.call(rbind, lapply(times, function(t) {
    short_days = unique(pts$day[pts$fit & len < t])
    p = pts[!pts$day %in% short_days, ]
    if (!nrow(p)) return(NULL)
    p$E_mV = vapply(p$datapath, emf_at, 0, read_at_s = t)
    r = calibrate(p, x, T_C, z)
    if (is.null(r)) NULL else data.frame(read_s = t, r$per_day[, c("day", "slope_mV_dec", "r2", "pct_nernst")])
  }))
  if (is.null(long)) return(NULL)
  later = long[long$read_s >= from_s, ]
  spread = do.call(rbind, lapply(split(later, later$day), function(d) data.frame(
    day = d$day[1], times = nrow(d), slope_min = min(d$slope_mV_dec), slope_max = max(d$slope_mV_dec),
    spread_mV_dec = diff(range(d$slope_mV_dec)), flagged = diff(range(d$slope_mV_dec)) > flag_mV)))
  list(long = long, spread = spread, flag_mV = flag_mV, from_s = from_s)
}

# One reading time for the whole study, chosen by stability rather than by the slope:
# the earliest time by which the given share of the fitted standards, over all days,
# had settled. Choosing the time that gives the best slope separately for each
# dataset would bias the result, so it is not offered.
suggest_read_time = function(settle, share = 0.9, times = c(60, 120, 180, 240, 300, 420, 600)) {
  f = settle[settle$fit, ]
  if (!nrow(f)) return(NULL)
  ok = vapply(times, function(t) mean(!is.na(f$settled_s) & f$settled_s <= t), 0)
  hit = which(ok >= share)
  data.frame(share_wanted = share, read_s = if (length(hit)) times[min(hit)] else NA,
             share_settled = if (length(hit)) ok[min(hit)] else max(ok), fitted_standards = nrow(f),
             never_settled = sum(is.na(f$settled_s)))
}

# A small, faint credit on every figure the app or the scripts make: the name only, in light grey, so it does not
# draw the eye in a slide or a paper figure
tool_name = "ElectroFast"
# In the app it follows the "Name on figures" box of the user's own session (unticked: no text at all). Reading the
# input inside a plot makes the plot redraw when the box changes; outside a reactive context (a download) it is
# read without that link. Scripts and tests have no session, so they keep the credit.
watermark = function(text = tool_name) {
  s = shiny::getDefaultReactiveDomain()
  on = if (is.null(s)) TRUE else tryCatch(s$input$fig_credit, error = function(e) shiny::isolate(s$input$fig_credit))
  if (isFALSE(on)) return(list())
  list(ggplot2::labs(caption = text), ggplot2::theme(plot.caption = ggplot2::element_text(size = 6, colour = "grey80", hjust = 1)))
}

# "Label plots" box: when ticked, CV plots get sweep-direction arrows and their peak values, and Multi-Potential
# Steps traces get the potential of each step, the cycle numbers and the reading windows. Read like watermark();
# scripts and tests have no session and stay unlabelled.
labels_on = function() {
  s = shiny::getDefaultReactiveDomain()
  if (is.null(s)) return(FALSE)
  isTRUE(tryCatch(s$input$fig_labels, error = function(e) shiny::isolate(s$input$fig_labels)))
}

# Every figure the app saves for download goes through fig_save: ggsave, then an invisible note of the software in
# the file itself, so a figure made with ElectroFast can be recognised even with "Name on figures" unticked.
# PNG: tEXt chunks "Software" and "Creation Time" (PNG specification, section 11.3.4), shown by file properties
# and exiftool. SVG: a <metadata> element. The note is lost if the figure is converted or screenshotted.
fig_credit_text = function() paste(tool_name, if (exists("electrofast_version")) electrofast_version)
png_crc32 = local({
  tab = vapply(0:255, function(n) { c = n; for (k in 1:8) c = if (bitwAnd(c, 1) == 1) bitwXor(0xEDB88320 - 2^32, bitwShiftR(c, 1)) else bitwShiftR(c, 1); c }, 0)
  function(bytes) { c = -1L
    for (b in as.integer(bytes)) c = bitwXor(tab[bitwAnd(bitwXor(c, b), 255) + 1], bitwShiftR(c, 8))
    bitwXor(c, -1L) }
})
png_text_chunk = function(key, value) {
  data = c(charToRaw("tEXt"), charToRaw(enc2native(key)), as.raw(0), charToRaw(iconv(value, "UTF-8", "latin1", sub = "?")))
  be32 = function(x) { x = as.numeric(x); if (x < 0) x = x + 2^32; as.raw(c(x %/% 2^24, (x %/% 2^16) %% 256, (x %/% 2^8) %% 256, x %% 256)) }
  c(be32(length(data) - 4), data, be32(png_crc32(data)))
}
fig_stamp = function(file) {
  ext = tolower(tools::file_ext(file)); note = fig_credit_text()
  if (ext == "png") {
    b = readBin(file, "raw", file.size(file))
    if (length(b) < 33 || !identical(b[13:16], charToRaw("IHDR"))) return(invisible(file))
    add = c(png_text_chunk("Software", note), png_text_chunk("Creation Time", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
    writeBin(c(b[1:33], add, b[-(1:33)]), file)   # after the 8-byte signature and the 25-byte IHDR chunk
  } else if (ext == "svg") {
    x = readLines(file, warn = FALSE, encoding = "UTF-8"); i = grep("<svg", x, fixed = TRUE)[1]
    if (!is.na(i)) { j = i; while (j < length(x) && !grepl(">", x[j], fixed = TRUE)) j = j + 1
      x = append(x, sprintf("<metadata>Made with %s</metadata>", note), after = j); writeLines(x, file, useBytes = TRUE) }
  }
  invisible(file)
}
fig_save = function(filename, ...) { ggplot2::ggsave(filename, ...); fig_stamp(filename) }

# Error notification that skips errors without a message. req() stops with an empty "silent" error; when a tryCatch
# catches one, showing it gave an empty red box.
notify_error = function(e, duration = 10) {
  m = conditionMessage(e)
  # the same message gets the same id, so Shiny replaces the pop-up instead of stacking a copy for every table and
  # plot that read the same file
  if (nzchar(trimws(m))) shiny::showNotification(m, type = "error", duration = duration, id = paste0("err_", gsub("[^A-Za-z0-9]", "", substr(m, 1, 80))))
  invisible(NULL)
}
