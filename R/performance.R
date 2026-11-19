# Sensor performance: the "S-criteria" of chemical sensors (sensitivity, selectivity, speed and stability), with
# specificity and precision, from result values of any technique.
# Each function takes a small table that the user prepares (one row per measurement) and returns plain results.
# ISE-specific methods (Nernstian slope, IUPAC selectivity methods, the check standard) stay under Potentiometry.

# Column names are matched loosely: case, spaces and units in brackets do not matter
perf_cols = function(d) { n = tolower(gsub("\\s*\\(.*\\)\\s*$", "", names(d))); names(d) = gsub("[^a-z0-9]+", "_", n); d }
perf_need = function(d, cols, what) {
  miss = setdiff(cols, names(d))
  if (length(miss)) stop(sprintf("%s needs the column(s) %s (download the template to see the layout).", what, paste(miss, collapse = ", ")))
  d
}
perf_num = function(v) suppressWarnings(as.numeric(v))

# Sensitivity: the calibration slope, its linear range, LOD and LOQ. The noise s is the SD of repeated blanks
# (concentration 0) when there are at least two, otherwise the residual SD of the line (ICH Q2: LOD = 3.3 s/slope,
# LOQ = 10 s/slope). With an electrode area, the slope per unit area is given as well.
perf_sensitivity = function(d, from = -Inf, to = Inf, area_cm2 = NA) {
  d = perf_need(perf_cols(d), c("concentration", "signal"), "Sensitivity")
  conc = perf_num(d$concentration); y = perf_num(d$signal)
  blank = y[!is.na(conc) & conc == 0 & !is.na(y)]
  blank_sd = if (length(blank) >= 2) stats::sd(blank) else NA_real_
  f = calibration_line(conc, y, from, to, blank_sd)
  if (is.null(f)) stop("Sensitivity needs at least three standards above zero in the chosen range.")
  k = !is.na(conc) & !is.na(y) & conc > 0 & conc >= from & conc <= to
  s = if (!is.na(blank_sd)) blank_sd else summary(stats::lm(y[k] ~ conc[k]))$sigma
  f$loq = 10 * s / abs(f$slope)
  f$slope_per_area = if (is.finite(area_cm2) && area_cm2 > 0) f$slope / area_cm2 else NA_real_
  f$blanks = length(blank)
  list(fit = f, points = data.frame(concentration = conc, signal = y, used = k))
}

# Selectivity: the change of the analyte signal caused by each interferent (analyte alone against analyte plus
# interferent), in percent, against a tolerance the user sets. Replicates of one interferent are averaged. This is the
# second selectivity procedure of Thévenot et al., Pure Appl. Chem. 1999, 71, 2333, section 5.2: the interferent at its
# expected concentration, the analyte at the middle of its expected range.
perf_selectivity = function(d, tolerance_pct = 5) {
  d = perf_need(perf_cols(d), c("interferent", "signal_analyte", "signal_mixture"), "Selectivity")
  d$a = perf_num(d$signal_analyte); d$m = perf_num(d$signal_mixture)
  conc = if ("interferent_concentration" %in% names(d)) as.character(d$interferent_concentration) else ""
  g = paste(d$interferent, conc)
  out = do.call(rbind, lapply(split(seq_len(nrow(d)), factor(g, levels = unique(g))), function(i) {
    a = mean(d$a[i], na.rm = TRUE); m = mean(d$m[i], na.rm = TRUE)
    data.frame(interferent = d$interferent[i[1]], interferent_concentration = conc[i[1]], replicates = length(i),
               signal_analyte = a, signal_mixture = m, change_pct = 100 * (m - a) / a)
  }))
  out$within_tolerance = abs(out$change_pct) <= tolerance_pct
  out
}

# Specificity (cross-reactivity): the response to each other species alone, at the same concentration as the target, as a
# percentage of the target response, and as the ratio of the target signal to that species' signal. The ratio is the first
# selectivity procedure of Thévenot et al., Pure Appl. Chem. 1999, 71, 2333, section 5.2 (no ratio when a species gives no
# positive signal). IUPAC recommends the graded term selectivity; specificity is absolute (Vessman et al., Pure Appl. Chem. 2001, 73, 1381).
perf_specificity = function(d) {
  d = perf_need(perf_cols(d), c("species", "signal", "target"), "Specificity")
  d$y = perf_num(d$signal); tg = toupper(trimws(as.character(d$target))) %in% c("TRUE", "T", "YES", "Y", "1")
  if (!any(tg)) stop("Specificity: mark the target species with yes in the target column.")
  ref = mean(d$y[tg], na.rm = TRUE)
  out = do.call(rbind, lapply(split(seq_len(nrow(d)), factor(d$species, levels = unique(d$species))), function(i)
    data.frame(species = d$species[i[1]], target = any(tg[i]), replicates = length(i), signal = mean(d$y[i], na.rm = TRUE))))
  out$cross_reactivity_pct = 100 * out$signal / ref
  out$selectivity_ratio = ifelse(out$signal > 0, ref / out$signal, NA_real_)
  out
}

# Stability: the signal over time (days, hours or cycles) as a percentage of the first value, the straight-line
# change per unit time, and when the signal first falls below the retention limit (90% by default), interpolated.
perf_stability = function(d, limit_pct = 90) {
  d = perf_need(perf_cols(d), c("time", "signal"), "Stability")
  tt = perf_num(d$time); y = perf_num(d$signal); k = !is.na(tt) & !is.na(y); tt = tt[k]; y = y[k]
  m = tapply(y, tt, mean); tm = as.numeric(names(m)); m = as.numeric(m)
  if (length(m) < 2) stop("Stability needs signals at two or more times.")
  ret = 100 * m / m[1]
  fit = stats::lm(ret ~ tm)
  below = which(ret < limit_pct)
  t_limit = if (length(below)) { j = below[1]; tm[j - 1] + (limit_pct - ret[j - 1]) * (tm[j] - tm[j - 1]) / (ret[j] - ret[j - 1]) } else NA_real_
  list(points = data.frame(time = tm, signal = m, retention_pct = ret),
       summary = data.frame(times = length(m), retained_at_last_pct = ret[length(ret)], change_pct_per_unit = unname(stats::coef(fit)[2]),
                            first_time_below_limit = t_limit, limit_pct = limit_pct, last_time = max(tm)))
}

# Speed: the response time after an addition at t_add. Baseline = median over pre_s before the addition, final level
# = median over the last post_s of the trace; the trace is smoothed with a running median over smooth_s, and t90 (t95)
# is the first time after the addition at which it has covered 90% (95%) of the change.
perf_speed = function(t, y, t_add, pre_s = 2, post_s = 10, smooth_s = 0.2) {
  o = order(t); t = t[o]; y = y[o]
  dt = stats::median(diff(t[seq_len(min(length(t), 10000))])); w = max(1, 2 * floor(smooth_s / dt / 2) + 1)
  ys = if (w > 1 && w < length(y)) as.numeric(stats::runmed(y, w, endrule = "median")) else y
  base = stats::median(y[t >= t_add - pre_s & t < t_add]); fin = stats::median(y[t >= max(t) - post_s])
  if (!is.finite(base) || !is.finite(fin)) stop("Speed: no points in the baseline or final window. Check the addition time.")
  ch = fin - base; after = t > t_add
  cross = function(f) { lev = base + f * ch; i = which(after & (if (ch > 0) ys >= lev else ys <= lev))[1]; if (is.na(i)) NA_real_ else t[i] - t_add }
  data.frame(t_add_s = t_add, baseline = base, final = fin, change = ch, t90_s = cross(0.9), t95_s = cross(0.95))
}

# A trace for the speed tab: EMF16 recordings as potential (mV) of one channel; other files through the instrument
# readers, as current when the file has a current, otherwise potential (V)
perf_trace = function(path, name, instrument = "auto", channel = NULL) {
  emf = tryCatch(read_emf16_grf(path)$long, error = function(e) NULL)
  if (!is.null(emf) && nrow(emf)) {
    chs = unique(emf$channel); ch = if (isTRUE(channel %in% chs)) channel else chs[1]; e = emf[emf$channel == ch, ]
    return(list(file = name, t = e$time_s, y = e$emf_mV, label = "Potential (mV)", scale = 1, channels = chs, channel = ch))
  }
  d = read_echem(path, name, instrument)$data
  if (is.null(d$t_s)) stop(name, ": the file has no time column")
  if (!is.null(d$I_A)) { u = current_unit(d$I_A)
    list(file = name, t = d$t_s, y = d$I_A, label = paste0("Current (", u$label, ")"), scale = u$scale, channels = NULL, channel = NULL) }
  else list(file = name, t = d$t_s, y = d$E_V, label = "Potential (V)", scale = 1, channels = NULL, channel = NULL)
}

# Precision in the terms of ICH Q2(R2): repeatability (same conditions), intermediate precision (within one
# laboratory, here across days and electrodes) and reproducibility (between laboratories). Variance components from
# a linear mixed model with the grouping columns present (day, electrode, laboratory) as random effects, fitted by
# REML (lme4). Electrode-to-electrode variation, often called reproducibility in sensor papers, is part of the
# intermediate precision here; reproducibility needs results from two or more laboratories.
perf_precision = function(d, min_levels = 5) {
  d = perf_need(perf_cols(d), c("value"), "Precision")
  d$value = perf_num(d$value); d = d[!is.na(d$value), , drop = FALSE]
  groups = intersect(c("day", "electrode", "laboratory"), names(d))
  groups = groups[vapply(groups, function(g) length(unique(stats::na.omit(d[[g]]))) >= 2, TRUE)]
  for (g in groups) d[[g]] = factor(d[[g]])
  notes = character(0); sd_g = setNames(rep(NA_real_, 3), c("day", "electrode", "laboratory"))
  if (length(groups)) {
    f = stats::as.formula(paste("value ~ 1 +", paste(sprintf("(1 | %s)", groups), collapse = " + ")))
    m = suppressMessages(lme4::lmer(f, data = d, REML = TRUE))
    vc = as.data.frame(lme4::VarCorr(m)); rep_sd = vc$sdcor[vc$grp == "Residual"]
    for (g in groups) sd_g[[g]] = vc$sdcor[vc$grp == g]
    if (lme4::isSingular(m)) notes = c(notes, "The model is singular: at least one variance component is estimated as zero, often because there are few groups.")
    for (g in groups) if (nlevels(d[[g]]) < min_levels) notes = c(notes, sprintf("Only %d %s levels: the between-%s SD is indicative only (five or more are usually needed).", nlevels(d[[g]]), g, g))
  } else { rep_sd = stats::sd(d$value); notes = c(notes, "No day, electrode or laboratory column with two or more levels: only the repeatability is given.") }
  z = function(x) ifelse(is.na(x), 0, x)
  inter = sqrt(rep_sd^2 + z(sd_g[["day"]])^2 + z(sd_g[["electrode"]])^2)
  repro = if (!is.na(sd_g[["laboratory"]])) sqrt(inter^2 + sd_g[["laboratory"]]^2) else NA_real_
  mu = mean(d$value)
  comp = data.frame(component = c("Repeatability (within-run SD)", "Between-day SD", "Between-electrode SD", "Intermediate precision (within the laboratory)",
                                  "Between-laboratory SD", "Reproducibility (between laboratories)"),
                    SD = c(rep_sd, sd_g[["day"]], sd_g[["electrode"]], if (any(c("day", "electrode") %in% groups)) inter else NA, sd_g[["laboratory"]], repro))
  comp$RSD_pct = 100 * comp$SD / abs(mu)
  list(components = comp[!is.na(comp$SD), ], mean = mu, n = nrow(d), groups = groups, notes = notes)
}

# Templates: the columns each table needs (and the optional ones), with no example values
perf_templates = list(
  sensitivity = data.frame(concentration = numeric(0), signal = numeric(0)),
  selectivity = data.frame(interferent = character(0), interferent_concentration = character(0), signal_analyte = numeric(0), signal_mixture = numeric(0)),
  specificity = data.frame(species = character(0), signal = numeric(0), target = character(0)),
  stability = data.frame(time = numeric(0), signal = numeric(0)),
  precision = data.frame(value = numeric(0), day = character(0), electrode = character(0), laboratory = character(0)))
