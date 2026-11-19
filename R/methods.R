# Analyses for voltammetry, amperometry, open-circuit potential, chronopotentiometry
# and rotating ring-disk electrodes. Every function takes a file read by read_echem()
# (or read_chi_txt / read_pine_txt) and returns plain data frames, so the app, the
# report and the tests all use the same numbers.
# Currents are in A, potentials in V, times in s, charges in C, concentrations in the
# user's unit (whatever the file names or the table say).

F_const = 96485.33212

# Concentration from file names
# Suggests a concentration for each file of a series; the user confirms it in a table.
# "+ 50mM EtOH" gives 50 (the analyte added to the electrolyte). Without a "+", a name
# with a single concentration gives that one. A file with no added concentration in a
# series where the others have one is taken as the blank (0). Anything else is NA.
conc_from_name = function(names, unit = "mM") {
  to = c(M = 1, mM = 1e-3, uM = 1e-6, nM = 1e-9)
  # units are matched in any case, because instruments often save names in lower case
  # ("5mm etoh" is 5 mM)
  canon = c(m = "M", mm = "mM", um = "uM", nm = "nM")
  one = function(n) {
    n = gsub("µ|μ", "u", n)
    pat = "[0-9]*\\.?[0-9]+(e-?[0-9]+)?\\s*([mun]?[mM])\\b"
    hit = regmatches(n, gregexpr(pat, n, ignore.case = TRUE))[[1]]
    plus = regmatches(n, regexpr(paste0("\\+\\s*", pat), n, ignore.case = TRUE))
    pick = if (length(plus)) sub("^\\+\\s*", "", plus) else if (length(hit) == 1) hit else NA
    if (is.na(pick)) return(c(value = NA, plus = length(plus) > 0))
    v = as.numeric(sub("\\s*[A-Za-z]+$", "", pick)); u = canon[[tolower(sub("^.*?([A-Za-z]+)$", "\\1", pick))]]
    c(value = v * to[[u]] / to[[unit]], plus = length(plus) > 0)
  }
  r = t(vapply(names, one, c(value = 0, plus = 0)))
  v = r[, "value"]
  # in a series where the analyte is written as "+ x", a file without "+" is the blank
  if (any(r[, "plus"] == 1)) v[r[, "plus"] == 0] = 0
  unname(v)
}
rpm_from_name = function(names) as.numeric(sub(".*?([0-9]+)\\s*rpm.*", "\\1", ifelse(grepl("[0-9]\\s*rpm", names, ignore.case = TRUE), names, NA), ignore.case = TRUE))

# Straight-line calibration with the metrics journals report: slope, intercept, R2,
# n, range, and a detection limit from the blank SD (3.3 s/slope, ICH Q2) when blank
# replicates exist, otherwise from the residual SD of the line.
calibration_line = function(conc, y, from = -Inf, to = Inf, blank_sd = NA) {
  k = !is.na(conc) & !is.na(y) & conc > 0 & conc >= from & conc <= to
  if (sum(k) < 3) return(NULL)
  m = stats::lm(y[k] ~ conc[k])
  s = unname(stats::coef(m)[2])
  s_yx = summary(m)$sigma
  data.frame(n = sum(k), from = min(conc[k]), to = max(conc[k]), slope = s,
             slope_se = summary(m)$coefficients[2, 2], intercept = unname(stats::coef(m)[1]),
             r2 = summary(m)$r.squared,
             lod = if (!is.na(blank_sd)) 3.3 * blank_sd / abs(s) else 3.3 * s_yx / abs(s),
             lod_basis = if (!is.na(blank_sd)) "3.3 x blank SD / slope" else "3.3 x residual SD / slope")
}

# Cyclic voltammetry
# Monotonic sweeps of a CV: start and end rows and direction
cv_sweeps = function(E) {
  d = sign(diff(E)); d[d == 0] = 1
  cuts = which(diff(d) != 0) + 1
  b = c(1, cuts, length(E))
  data.frame(from = b[-length(b)], to = b[-1], dir = ifelse(E[b[-1]] < E[b[-length(b)]], "cathodic", "anodic"))
}
# The last full sweep in one direction (spanning more than min_span V), which is the
# most stable cycle of a repeated CV
last_sweep = function(x, dir = "cathodic", min_span = 0.5) {
  E = x$data$E_V; I = x$data$I_A
  s = cv_sweeps(E)
  s = s[s$dir == dir & abs(E[s$to] - E[s$from]) > min_span, ]
  if (!nrow(s)) stop("No ", dir, " sweep longer than ", min_span, " V in ", x$file)
  r = s[nrow(s), ]
  data.frame(E_V = E[r$from:r$to], I_A = I[r$from:r$to])
}
# Charge of a peak between two potentials on one sweep, above a straight baseline
# joining the two window ends: Q = |integral of (I - baseline) dE| / scan rate
peak_charge = function(sweep, window, v) {
  w = sweep[sweep$E_V >= min(window) & sweep$E_V <= max(window), ]
  base = stats::approx(c(w$E_V[1], w$E_V[nrow(w)]), c(w$I_A[1], w$I_A[nrow(w)]), xout = w$E_V)$y
  net = w$I_A - base
  q = abs(sum(diff(w$E_V) * (net[-1] + net[-length(net)]) / 2)) / v
  list(Q_C = q, window = cbind(w, base_A = base))
}
# Surface area from a reduction peak with a known charge density: gold oxide 390,
# PdO 405 (uC/cm2), or any other value the user enters
ecsa_window = function(x, window, q_ref_uC_cm2, geometric_area_cm2 = NA, dir = "cathodic") {
  v = x$params$scan_rate_V_s
  s = last_sweep(x, dir)
  pc = peak_charge(s, window, v)
  ecsa = pc$Q_C / (q_ref_uC_cm2 * 1e-6)
  list(result = data.frame(file = x$file, scan_rate_V_s = v, window = sprintf("%.3f to %.3f V", min(window), max(window)),
                           Q_uC = pc$Q_C * 1e6, q_ref_uC_cm2 = q_ref_uC_cm2, ECSA_cm2 = ecsa,
                           roughness_factor = ecsa / geometric_area_cm2),
       sweep = s, window = pc$window)
}
# Peak current against scan rate for both peaks, before the log-log slope: straight-line fits of |ip| against v
# (linear for a surface-confined species) and against v^1/2 (linear through the origin for a diffusing species),
# each with slope, intercept, their standard errors and R2, then the log-log slope (0.5 for diffusion, 1 for an
# adsorbed species). Needs at least 3 scan rates for each peak.
scan_rate_fits = function(v, ipa_A, ipc_A) {
  one = function(peak, y) {
    ok = is.finite(v) & is.finite(y) & v > 0 & y != 0; vv = v[ok]; yy = abs(y[ok])
    if (length(unique(vv)) < 3) return(NULL)
    lin = function(x, against) { m = stats::lm(yy ~ x); cf = summary(m)$coefficients
      data.frame(peak = peak, against = against, n = length(yy), slope = cf[2, 1], slope_se = cf[2, 2],
                 slope_unit = if (against == "v") "A s/V" else "A (s/V)^1/2",
                 intercept_A = cf[1, 1], intercept_se_A = cf[1, 2], intercept_pct_of_smallest_ip = 100 * cf[1, 1] / min(yy), r2 = summary(m)$r.squared) }
    ll = stats::lm(log10(yy) ~ log10(vv)); cl = summary(ll)$coefficients
    list(fits = rbind(lin(vv, "v"), lin(sqrt(vv), "√v")),
         loglog = data.frame(peak = peak, n = length(yy), loglog_slope = cl[2, 1], loglog_slope_se = cl[2, 2], r2 = summary(ll)$r.squared))
  }
  r = Filter(Negate(is.null), list(one("anodic", ipa_A), one("cathodic", ipc_A)))
  if (!length(r)) return(NULL)
  list(fits = do.call(rbind, lapply(r, `[[`, "fits")), loglog = do.call(rbind, lapply(r, `[[`, "loglog")))
}
# What the scan-rate fits point to, one line per peak: which value the log-log slope is closer to, and the intercept of
# the line against v^1/2 as a share of the smallest peak current. It is flagged when it differs from zero by more than
# 2 standard errors and is larger than 5% of the smallest peak (a default of this app, not a published criterion);
# with noise-free data the standard error alone would flag intercepts of a few parts per million.
scan_rate_reading = function(sf) {
  if (is.null(sf)) return("At least 3 scan rates are needed for each peak.")
  vapply(seq_len(nrow(sf$loglog)), function(i) { L = sf$loglog[i, ]; q = sf$fits[sf$fits$peak == L$peak & sf$fits$against == "√v", ]
    sprintf("%s peak: log-log slope %.3f ± %.3f, closer to %s. Against √v the intercept is %s.",
            tools::toTitleCase(L$peak), L$loglog_slope, L$loglog_slope_se,
            if (L$loglog_slope < 0.75) "0.5 (diffusion-controlled)" else "1 (surface-confined species)",
            if (abs(q$intercept_A) <= 2 * q$intercept_se_A || abs(q$intercept_pct_of_smallest_ip) <= 5)
              sprintf("%.2g%% of the smallest peak, close to zero as expected for diffusion", q$intercept_pct_of_smallest_ip)
            else sprintf("%.3g µA (%.2g%% of the smallest peak): charging current, background or slow kinetics may contribute", q$intercept_A * 1e6, q$intercept_pct_of_smallest_ip))
  }, "")
}
# Tafel slope over one potential window of a slow LSV. Because the answer depends on
# the window, tafel_scan() fits every window of a few widths so the choice is visible.
tafel_fit = function(x, window) {
  d = x$data; k = d$E_V >= min(window) & d$E_V <= max(window) & d$I_A > 0
  if (sum(k) < 3) return(NULL)
  m = stats::lm(log10(d$I_A[k]) ~ d$E_V[k])
  s = unname(stats::coef(m)[2])
  data.frame(from_V = min(window), to_V = max(window), n = sum(k), tafel_mV_dec = 1000 / s, r2 = summary(m)$r.squared)
}
tafel_scan = function(x, from = min(x$data$E_V), to = max(x$data$E_V), widths = c(0.006, 0.020, 0.050), step = 0.002, min_r2 = 0.99) {
  out = do.call(rbind, lapply(widths, function(w) do.call(rbind, lapply(seq(from, to - w, by = step), function(lo) {
    t = tafel_fit(x, c(lo, lo + w)); if (!is.null(t)) cbind(width_mV = w * 1000, t) }))))
  if (is.null(out)) return(NULL)
  out[out$r2 > min_r2 & out$tafel_mV_dec > 0 & out$tafel_mV_dec < 400, ]
}

# Pulse voltammetry (DPV, SWV)
# The response of one voltammogram: the current at a set potential (nearest point), or
# the largest current within a potential window
pulse_response = function(x, mode = c("at_potential", "peak"), E_read = NA, window = NULL) {
  mode = match.arg(mode); d = x$data
  if (mode == "at_potential") { k = which.min(abs(d$E_V - E_read)); return(data.frame(file = x$file, E_V = d$E_V[k], I_A = d$I_A[k])) }
  w = d[d$E_V >= min(window) & d$E_V <= max(window), ]
  k = which.max(abs(w$I_A))
  data.frame(file = x$file, E_V = w$E_V[k], I_A = w$I_A[k])
}

# Open-circuit potential
ocp_summary = function(x, last_s = 60) {
  d = x$data; tail_ = d[d$t_s >= max(d$t_s) - last_s, ]
  drift = unname(stats::coef(stats::lm(tail_$E_V ~ tail_$t_s))[2]) * 1000 * 60
  data.frame(file = x$file, run_s = max(d$t_s), E_end_V = d$E_V[nrow(d)], E_min_V = min(d$E_V), E_max_V = max(d$E_V),
             E_mean_last_V = mean(tail_$E_V), drift_last_mV_min = drift)
}

# Chronoamperometry
# Charge passed in an i-t record (for example a deposition step), and the same charge
# expressed as monolayers when an area and a monolayer charge density are given
charge_passed = function(x, area_cm2 = NA, q_ml_uC_cm2 = NA) {
  d = x$data
  q = abs(sum(diff(d$t_s) * (d$I_A[-1] + d$I_A[-nrow(d)]) / 2))
  data.frame(file = x$file, E_V = chi_setting(x, "Init E (V)"), time_s = max(d$t_s), Q_uC = q * 1e6,
             I_start_uA = d$I_A[1] * 1e6, I_end_uA = mean(utils::tail(d$I_A, 100)) * 1e6,
             monolayers = q / (q_ml_uC_cm2 * 1e-6 * area_cm2))
}
# Mean current over the last avg_s seconds before each time in t_ends
current_at_ends = function(d, t_ends, avg_s) vapply(t_ends, function(te) {
  k = d$t_s > te - avg_s & d$t_s <= te
  if (any(k)) mean(d$I_A[k]) else NA_real_
}, 0)
# Step chronoamperometry (CHI "Chronoamperometry", one or two steps of equal width):
# the current at the end of each step
step_currents = function(x, avg_s = 2) {
  w = chi_setting(x, "Pulse Width (sec)"); n = chi_setting(x, "Step", 1)
  i = current_at_ends(x$data, w * seq_len(n), avg_s)
  stats::setNames(data.frame(x$file, t(i)), c("file", paste0("I_step", seq_len(n), "_A")))
}
# Multi-potential steps (CHI "Multi-Potential Steps"): steps E1, T1, E2, T2, ... repeated
# for a number of cycles. The response is the current at the end of the sensing step
# (by default the last step) in each cycle. The first cycle is left out of the mean
# because the surface is still changing; the spread over the others gives the RSD.
multistep_currents = function(x, sense_step = NULL, avg_s = 1, skip_cycles = 1) {
  st = x$params$settings
  k = sort(as.integer(sub("^E([0-9]+) \\(V\\)$", "\\1", grep("^E[0-9]+ \\(V\\)$", names(st), value = TRUE))))
  E = vapply(k, function(j) st[[sprintf("E%d (V)", j)]], 0)
  Tt = vapply(k, function(j) st[[sprintf("T%d (s)", j)]], 0)
  n = chi_setting(x, "Cycle", 1)
  if (is.null(sense_step)) sense_step = length(k)
  per = sum(Tt)
  ends = (seq_len(n) - 1) * per + sum(Tt[seq_len(sense_step)])
  i = current_at_ends(x$data, ends, avg_s)
  use = i[-seq_len(min(skip_cycles, n - 1))]
  data.frame(file = x$file, sense_E_V = E[sense_step], steps = paste(sprintf("%g V %g s", E, Tt), collapse = ", "),
             cycles = n, I_mean_A = mean(use), I_sd_A = if (length(use) > 1) stats::sd(use) else NA,
             rsd_pct = if (length(use) > 1) 100 * stats::sd(use) / abs(mean(use)) else NA,
             I_cycles = paste(signif(i, 4), collapse = "; "))
}

# Chronopotentiometry
# Transition time. At the start the potential moves fast while the double layer charges
# and the surface reaches the redox couple; it then slows on the plateau of the couple,
# and moves fast again at the transition, when the couple is used up at the surface.
# The transition is therefore the steepest point that comes after the slope has dipped:
# a point counts only if its slope is at least rise times the lowest slope before it.
# When the plateau is too short to separate from the start, no transition is reported.
# The run must also have settled afterwards (end slope below a tenth of the transition
# slope). E at tau/4 is the quarter-wave potential, close to E1/2 for a reversible couple.
#
# The default (steepest point) matches the usual reading of tau as the point of greatest
# slope (Pine Research, chronopotentiometry support article). method = "tangent" is a
# tangent-intersection construction written for comparison; no published procedure for it
# has been checked yet, so treat it as an alternative reading, not a reference method.
# It draws a tangent on the plateau (at its
# flattest point before the inflection) and a tangent at the inflection; tau is the time
# where they cross, which marks the start of the transition and so comes earlier than
# the inflection. The time where the inflection tangent crosses a tangent on the flattest
# part after the transition (within three inflection times) is given as tau_end.
cp_transition = function(x, skip_s = 0.3, smooth_s = 0.5, rise = 1.1, min_step_V = 0.05, method = c("inflection", "tangent")) {
  method = match.arg(method)
  d = x$data[!duplicated(x$data$t_s), ]
  dt = stats::median(diff(d$t_s)); k = max(3, round(smooth_s / dt))
  Es = stats::filter(d$E_V, rep(1 / k, k), sides = 2)
  s = abs(c(NA, diff(Es)) / dt)
  use = which(d$t_s > skip_s & !is.na(s))
  low_before = cummin(s[use])
  cand = use[s[use] >= rise * low_before]
  i_app = abs(chi_setting(x, "Cathodic Current (A)", NA))
  j = if (length(cand)) cand[which.max(s[cand])] else NA
  end_slope = mean(utils::tail(s[!is.na(s)], k))
  # a real transition moves the potential by tens of mV between tau/2 and 2 tau; a bump
  # in a slow drift does not
  step = if (is.na(j)) 0 else abs(diff(stats::approx(d$t_s, d$E_V, c(d$t_s[j] / 2, min(2 * d$t_s[j], max(d$t_s))))$y))
  ok = !is.na(j) && end_slope < 0.1 * s[j] && d$t_s[j] < 0.9 * max(d$t_s) && step >= min_step_V
  tau = if (ok) d$t_s[j] else NA
  tau_end = NA
  if (ok && method == "tangent") {
    sg = c(NA, diff(Es)) / dt                         # signed slope of the smoothed curve
    line = function(i) c(a = Es[i] - sg[i] * d$t_s[i], b = sg[i])
    cross = function(L1, L2) unname((L2["a"] - L1["a"]) / (L1["b"] - L2["b"]))
    pre = which(d$t_s > skip_s & d$t_s < d$t_s[j] & !is.na(sg))
    post = which(d$t_s > d$t_s[j] & d$t_s <= 3 * d$t_s[j] & !is.na(sg))
    Lj = line(j)
    tau = cross(line(pre[which.min(abs(sg[pre]))]), Lj)
    if (length(post)) tau_end = cross(Lj, line(post[which.min(abs(sg[post]))]))
  }
  data.frame(file = x$file, method = method, current_A = i_app, run_s = max(d$t_s), tau_s = tau, tau_end_s = tau_end,
             E_at_tau_V = if (ok) stats::approx(d$t_s, d$E_V, tau)$y else NA,
             E_quarter_V = if (ok) stats::approx(d$t_s, d$E_V, tau / 4)$y else NA,
             step_mV = 1000 * step,
             note = if (ok) "" else if (is.na(j)) "no separate transition: the plateau is too short to tell it from the start of the run"
                    else if (step < min_step_V) "no clear transition: no sharp potential step after the plateau"
                    else "no clear transition: the potential had not settled again by the end of the run")
}
# Sand equation: i tau^1/2 / C is constant for a diffusion-controlled process. With the
# electrode area, D follows from i tau^1/2 = n F A C (pi D)^1/2 / 2.
sand_table = function(tr, conc_mM, n = 1, area_cm2 = NA) {
  C = conc_mM * 1e-6                                  # mol/cm3
  it = tr$current_A * sqrt(tr$tau_s)
  cbind(tr[, c("file", "current_A", "tau_s")], conc_mM = conc_mM, i_sqrt_tau = it, i_sqrt_tau_per_C = it / conc_mM,
        D_cm2_s = (2 * it / (n * F_const * area_cm2 * C))^2 / pi)
}
# Sand line across concentrations: i tau^1/2 against C. A non-zero intercept points to
# a contribution that does not scale with concentration (double-layer charging,
# adsorption), which is why D is taken from the slope rather than from single files.
sand_fit = function(tab, n = 1, area_cm2 = NA) {
  k = !is.na(tab$tau_s) & !is.na(tab$conc_mM)
  if (sum(k) < 2) return(NULL)
  m = stats::lm(tab$i_sqrt_tau[k] ~ tab$conc_mM[k])
  s = unname(stats::coef(m)[2])                         # A s^1/2 per mM
  data.frame(n_files = sum(k), slope_A_s05_per_mM = s, intercept_A_s05 = unname(stats::coef(m)[1]),
             r2 = summary(m)$r.squared,
             D_from_slope_cm2_s = (2 * s / (n * F_const * area_cm2 * 1e-6))^2 / pi)
}

# Rotating ring-disk electrode
# Collection efficiency N = -I_ring / I_disk over a disk-potential window on the
# limiting-current plateau. The two currents have opposite signs whatever the
# instrument's sign convention, so the ratio is positive.
rrde_collection = function(x, disk_window, rpm = NA) {
  d = x$data
  if (is.null(d$I2_A)) stop(x$file, " has no ring current (second channel)")
  w = d[d$E_V >= min(disk_window) & d$E_V <= max(disk_window), ]
  n = -w$I2_A / w$I_A
  data.frame(file = x$file, rpm = rpm, N = mean(n), N_sd = stats::sd(n), points = nrow(w),
             I_disk_A = mean(w$I_A), I_ring_A = mean(w$I2_A),
             ring_E_V = if (!is.null(d$E2_V)) stats::median(d$E2_V) else chi_setting(x, "2nd E (V)", NA))
}
# Levich check across rotation rates: limiting disk current against omega^1/2
levich = function(cal) {
  cal = cal[!is.na(cal$rpm), ]
  if (nrow(cal) < 2) return(NULL)
  w = sqrt(2 * pi * cal$rpm / 60)
  m = stats::lm(abs(cal$I_disk_A) ~ w)
  data.frame(n = nrow(cal), slope_A_per_rad_s = unname(stats::coef(m)[2]), intercept_A = unname(stats::coef(m)[1]),
             r2 = summary(m)$r.squared)
}
# Faradaic efficiency of a product made at the disk and collected at the ring:
# FE = |I_ring| / (N |I_disk|), point by point and averaged over a disk-potential window.
# An optional ring baseline (the ring current where the disk makes no product) is
# subtracted first.
rrde_fe = function(x, N, disk_window, ring_baseline_window = NULL) {
  d = x$data
  if (is.null(d$I2_A)) stop(x$file, " has no ring current (second channel)")
  base = if (is.null(ring_baseline_window)) 0 else
    stats::median(d$I2_A[d$E_V >= min(ring_baseline_window) & d$E_V <= max(ring_baseline_window)])
  d$FE_pct = 100 * abs(d$I2_A - base) / (N * abs(d$I_A))
  w = d[d$E_V >= min(disk_window) & d$E_V <= max(disk_window), ]
  list(summary = data.frame(file = x$file, N = N, window = sprintf("%.3f to %.3f V", min(disk_window), max(disk_window)),
                            ring_baseline_A = base, FE_mean_pct = mean(w$FE_pct), FE_sd_pct = stats::sd(w$FE_pct),
                            points = nrow(w)),
       curve = d)
}
