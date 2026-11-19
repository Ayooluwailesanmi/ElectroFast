# Electrochemical impedance: reading, checking, plotting data, equivalent circuits written as text,
# fitting, ranking, Kramers-Kronig validation and cautious interpretation.
#
# Circuits are written in a short notation: elements joined by "-" are in series and p(a, b, ...)
# holds elements in parallel, for example "Rs-p(Rct,Q)" or "Rs-p(Q,Rct-W)". Elements:
#   R  resistor (ohm)                 C  capacitor (F)
#   Q  constant phase element, Z = 1 / (Q (j w)^n); parameters Q (S s^n) and n
#   W  semi-infinite Warburg, Z = sigma w^-1/2 (1 - j); parameter sigma (ohm s^-1/2)
#   Wo finite Warburg, reflective (open): R coth(s)/s, s = (j w tau)^1/2; parameters R and tau
#   Ws finite Warburg, transmissive (short): R tanh(s)/s; parameters R and tau
#   L  inductor (H)
# An element can carry a label after its letter (Rs, Rct, Q1, Wo2, ...), which names its parameters:
# Q1 has Q1 and n1, W has sigma, W2 has sigma2, Wo1 has R_Wo1 and tau_Wo1.
# Zim is the imaginary part as measured (negative for a capacitor); plots show -Z''.

eis_kinds = c("Wo", "Ws", "R", "C", "Q", "W", "L")

# Parse the notation into a tree of series and parallel nodes
eis_parse = function(s) {
  s = gsub("[[:space:]]", "", s)
  split_top = function(x, sep) {             # split at sep where not inside brackets
    out = character(); depth = 0; start = 1
    for (i in seq_len(nchar(x))) { ch = substr(x, i, i)
      if (ch == "(") depth = depth + 1 else if (ch == ")") depth = depth - 1
      else if (ch == sep && depth == 0) { out = c(out, substr(x, start, i - 1)); start = i + 1 } }
    c(out, substr(x, start, nchar(x)))
  }
  node = function(x) {
    parts = split_top(x, "-")
    if (length(parts) > 1) return(list(type = "series", kids = lapply(parts, node)))
    if (grepl("^p\\(.*\\)$", x)) return(list(type = "parallel", kids = lapply(split_top(substr(x, 3, nchar(x) - 1), ","), node)))
    if (grepl("^\\(.*\\)$", x)) return(node(substr(x, 2, nchar(x) - 1)))
    k = eis_kinds[startsWith(x, eis_kinds)][1]
    if (is.na(k)) stop("unknown element '", x, "' (use R, C, Q, W, Wo, Ws or L, with an optional label)")
    lab = substr(x, nchar(k) + 1, nchar(x))
    if (grepl("[^A-Za-z0-9_]", lab)) stop("element '", x, "': a label may only hold letters, digits and _")
    list(type = "element", kind = k, name = x, label = lab)
  }
  node(s)
}
eis_elements = function(t) if (t$type == "element") list(t) else do.call(c, lapply(t$kids, eis_elements))
eis_param_names = function(el) switch(el$kind, R = el$name, C = el$name, L = el$name,
  Q = c(el$name, paste0("n", el$label)), W = paste0("sigma", el$label),
  Wo = c(paste0("R_", el$name), paste0("tau_", el$name)), Ws = c(paste0("R_", el$name), paste0("tau_", el$name)))

# Compile a circuit: parameter names, which are CPE exponents, the impedance function and a description
eis_compile = function(s) {
  tree = eis_parse(s); els = eis_elements(tree)
  if (anyDuplicated(vapply(els, `[[`, "", "name"))) stop("circuit '", s, "': each element needs its own label")
  pn = unlist(lapply(els, eis_param_names)); idx = list(); k = 0
  for (e in els) { m = length(eis_param_names(e)); idx[[e$name]] = k + seq_len(m); k = k + m }
  z = function(t, w, p) switch(t$type,
    series = Reduce(`+`, lapply(t$kids, z, w = w, p = p)),
    parallel = 1 / Reduce(`+`, lapply(t$kids, function(k) 1 / z(k, w, p))),
    element = { q = p[idx[[t$name]]]
      switch(t$kind, R = complex(real = q[1], imaginary = 0 * w), C = 1 / (1i * w * q[1]), L = 1i * w * q[1],
             Q = 1 / (q[1] * (1i * w)^q[2]), W = q[1] * (1 - 1i) / sqrt(w),
             Wo = { u = sqrt(1i * w * q[2]); q[1] / (tanh(u) * u) }, Ws = { u = sqrt(1i * w * q[2]); q[1] * tanh(u) / u }) })
  list(string = s, tree = tree, elements = els, p = pn, exponent = grepl("^n", pn) & pn %in% unlist(lapply(els, function(e) if (e$kind == "Q") paste0("n", e$label))),
       f = function(w, p) z(tree, w, p), notation = eis_notation(tree))
}
eis_notation = function(t) switch(t$type, element = t$name,
  series = paste(vapply(t$kids, eis_notation, ""), collapse = " + "),
  parallel = paste0("(", paste(vapply(t$kids, eis_notation, ""), collapse = " || "), ")"))   # plain || so every font and PDF can show it

# Candidate circuits tried by default (the first four are the ones listed in the EIS requirements)
eis_library = c("Rs-p(Rct,C)", "Rs-p(Rct,Q)", "Rs-p(C,Rct-W)", "Rs-Q", "Rs-C", "Rs-p(Q,Rct-W)", "Rs-p(Q,Rct-Wo)",
                "Rs-p(Q,Rct-Ws)", "Rs-p(R1,Q1)-p(R2,Q2)", "Rs-p(Q1,R1-p(R2,Q2))", "Rs-p(Rct,Q)-W")

# Starting values from the shape of the spectrum: series R from the high-frequency real part, other
# resistances from the arc width, capacitances from the frequency at the top of the arc, sigma from the
# low-frequency end, tau from the lowest frequency
eis_start = function(d, cc) {
  w = 2 * pi * d$freq_Hz; Zr = Re(d$Z); Zi = Im(d$Z)
  Rs = max(min(Zr), 1e-6); Rp = max(max(Zr) - Rs, Rs); Cg = 1 / (w[which.max(-Zi)] * Rp)
  sg = max(abs(Zi[which.min(w)]) * sqrt(min(w)), 1e-6)
  top = if (cc$tree$type == "series") vapply(cc$tree$kids, function(k) if (k$type == "element") k$name else "", "") else ""
  nR = sum(vapply(cc$elements, function(e) e$kind == "R" && !e$name %in% top[1], TRUE))
  v = unlist(lapply(cc$elements, function(e) switch(e$kind,
    R = if (e$name %in% top) Rs else Rp / max(nR, 1), C = Cg, L = 1e-7, Q = c(Cg, 0.9), W = sg, Wo = c(Rp, 1 / min(w)), Ws = c(Rp, 1 / min(w)))))
  stats::setNames(v, cc$p)
}
# "R1=1000, n=0.9" or "R1=1e3:1e8" into a named vector or a two-column matrix
eis_named = function(txt) { if (is.null(txt) || !nzchar(trimws(txt))) return(NULL)
  kv = strsplit(trimws(strsplit(txt, "[,;]")[[1]]), "=")
  stats::setNames(vapply(kv, function(x) trimws(x[2]), ""), vapply(kv, function(x) trimws(x[1]), "")) }

# Weighted complex nonlinear least squares. Weighting: "modulus" (1/|Z| on both parts, the default),
# "proportional" (1/|Z'| and 1/|Z''| separately) or "unit" (none). Parameters are fitted on a log scale
# (CPE exponents on a logistic scale between 0 and 1), from several starts. Optional starting values and
# bounds as text, for example start = "Rct=2e6" and bounds = "Rct=1e5:1e8, n=0.5:1".
eis_fit = function(d, circuit, weighting = c("modulus", "proportional", "unit"), start = NULL, bounds = NULL, starts = 12, seed = 1) {
  weighting = match.arg(weighting); cc = eis_compile(circuit); pn = cc$p; ex = cc$exponent
  w = 2 * pi * d$freq_Hz; Zm = d$Z
  wr = switch(weighting, modulus = 1 / Mod(Zm), proportional = 1 / pmax(abs(Re(Zm)), 1e-12), unit = rep(1, length(w)))
  wi = switch(weighting, modulus = 1 / Mod(Zm), proportional = 1 / pmax(abs(Im(Zm)), 1e-12), unit = rep(1, length(w)))
  lo = ifelse(ex, 0, 0); hi = ifelse(ex, 1, Inf)
  b = eis_named(bounds)
  for (k in intersect(names(b), pn)) { r = as.numeric(strsplit(b[[k]], ":")[[1]]); lo[pn == k] = r[1]; hi[pn == k] = r[2] }
  to_p = function(u) ifelse(ex, lo + (pmin(hi, 1) - lo) / (1 + exp(-u)), exp(u))
  to_u = function(p) ifelse(ex, -log((pmin(hi, 1) - lo) / pmin(pmax(p - lo, 1e-6), pmin(hi, 1) - lo - 1e-6) - 1), log(p))
  ssr = function(u) { p = to_p(u)
    if (any(!ex & (p < lo | p > hi))) return(1e30)
    z = cc$f(w, p); s = sum(((Re(z) - Re(Zm)) * wr)^2 + ((Im(z) - Im(Zm)) * wi)^2); if (is.finite(s)) s else 1e30 }
  p0 = eis_start(d, cc); s0 = eis_named(start)
  for (k in intersect(names(s0), pn)) p0[[k]] = as.numeric(s0[[k]])
  p0 = ifelse(!ex & is.finite(hi), pmin(pmax(p0, lo), hi), pmax(p0, ifelse(ex, 0, lo)))
  set.seed(seed); best = NULL
  for (k in seq_len(starts)) {
    u0 = to_u(p0) + if (k == 1) 0 else stats::rnorm(length(pn), 0, ifelse(ex, 0.5, 1.5))
    o = tryCatch(stats::optim(u0, ssr, method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-12)), error = function(e) NULL)
    if (!is.null(o)) o = tryCatch(stats::optim(o$par, ssr, method = "BFGS", control = list(maxit = 2000, reltol = 1e-14)), error = function(e) o)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best = o
  }
  N = 2 * length(w); np = length(pn); p = to_p(best$par); s2 = best$value / max(N - np, 1)
  # standard errors from the Hessian of the weighted sum of squares, on the fitting scale, carried to the
  # parameter scale by the delta method
  H = tryCatch(stats::optimHess(best$par, ssr), error = function(e) NULL)
  se_u = if (!is.null(H)) tryCatch(sqrt(pmax(diag(solve(H / 2)) * s2, 0)), error = function(e) rep(NA, np)) else rep(NA, np)
  dpdu = ifelse(ex, (p - lo) * (1 - (p - lo) / (pmin(hi, 1) - lo)), p)
  se = abs(dpdu) * se_u
  fz = cc$f(w, p); res = fz - Zm
  at_bound = (ex & (p > pmin(hi, 1) - 0.005 | p < lo + 0.005)) | (!ex & ((is.finite(hi) & p > 0.99 * hi) | (lo > 0 & p < 1.01 * lo)))
  list(circuit = circuit, notation = cc$notation, weighting = weighting,
       params = data.frame(parameter = pn, value = p, se = se, rel_se_pct = 100 * se / abs(p), at_bound = at_bound,
                           unit = eis_param_unit(pn, cc)),
       start_values = p0, chi2 = best$value, chi2_red = s2, n_points = N, n_params = np,
       rmse_ohm = sqrt(mean(c(Re(res), Im(res))^2)), nrmse_pct = 100 * sqrt(mean(c(Re(res), Im(res))^2 / Mod(Zm)^2)),
       aic = N * log(best$value / N) + 2 * np, aicc = N * log(best$value / N) + 2 * np + 2 * np * (np + 1) / max(N - np - 1, 1),
       bic = N * log(best$value / N) + np * log(N), converged = best$convergence == 0, fit_Z = fz,
       residuals = data.frame(freq_Hz = d$freq_Hz, res_re_pct = 100 * Re(res) / Mod(Zm), res_im_pct = 100 * Im(res) / Mod(Zm)))
}
eis_param_unit = function(pn, cc) vapply(pn, function(k) {
  e = cc$elements[[which(vapply(cc$elements, function(e) k %in% eis_param_names(e), TRUE))[1]]]
  switch(e$kind, R = "ohm", C = "F", L = "H", W = "ohm s^-1/2", Q = if (startsWith(k, "n")) "" else "S s^n",
         Wo = , Ws = if (startsWith(k, "tau")) "s" else "ohm")
}, "")

# Linear Kramers-Kronig test (Schonleber et al., Electrochim. Acta 2014, 131, 20): a series resistance,
# M RC elements with time constants spread evenly on a log scale over the measured range, a series
# capacitance and an inductance, all linear in the unknowns. M grows until mu = 1 - sum|R_k<0| / sum|R_k>=0|
# falls below 0.85. Small, structureless residuals mean the data behave as a linear, causal, stable system.
lin_kk = function(d, c_mu = 0.85, M_max = 60) {
  w = 2 * pi * d$freq_Hz; Z = d$Z; wt = 1 / Mod(Z)
  solve_M = function(M) {
    tau = 10^seq(log10(1 / max(w)), log10(1 / min(w)), length.out = M)
    Xr = cbind(1, sapply(tau, function(t) 1 / (1 + (w * t)^2)), 0, 0)
    Xi = cbind(0, sapply(tau, function(t) -w * t / (1 + (w * t)^2)), -1 / w, w)
    b = qr.solve(rbind(Xr * wt, Xi * wt), c(Re(Z) * wt, Im(Z) * wt))
    Rk = b[2:(M + 1)]
    list(M = M, mu = 1 - sum(abs(Rk[Rk < 0])) / sum(abs(Rk[Rk >= 0])), fit = complex(real = drop(Xr %*% b), imaginary = drop(Xi %*% b)))
  }
  for (M in seq_len(M_max)) { s = solve_M(M); if (M > 1 && s$mu <= c_mu) break }
  rr = (Re(Z) - Re(s$fit)) / Mod(Z); ri = (Im(Z) - Im(s$fit)) / Mod(Z)
  list(M = s$M, mu = s$mu, residuals = data.frame(freq_Hz = d$freq_Hz, res_re_pct = 100 * rr, res_im_pct = 100 * ri),
       max_abs_res_pct = 100 * max(abs(c(rr, ri))), rms_res_pct = 100 * sqrt(mean(c(rr, ri)^2)))
}

# Checks on the data as read: missing or repeated frequencies, sign convention, and agreement of
# |Z| and phase recomputed from Z' and Z'' with the columns in the file when it has them
eis_validate = function(e) {
  d = e$data; out = list()
  out$points = nrow(d); out$freq_range = sprintf("%g to %g Hz", min(d$freq_Hz), max(d$freq_Hz))
  out$repeated_frequencies = sum(duplicated(d$freq_Hz))
  out$missing_or_nonfinite = if (!is.null(e$dropped)) e$dropped else 0
  out$share_capacitive = mean(Im(d$Z) < 0)
  if (!is.null(e$file_mod)) out$max_mod_diff_pct = 100 * max(abs(Mod(d$Z) - e$file_mod) / e$file_mod)
  if (!is.null(e$file_phase)) out$max_phase_diff_deg = max(abs(abs(Arg(d$Z) * 180 / pi) - abs(e$file_phase)))
  flags = character()
  if (out$repeated_frequencies) flags = c(flags, sprintf("%d repeated frequencies", out$repeated_frequencies))
  if (out$missing_or_nonfinite) flags = c(flags, sprintf("%d rows without usable numbers were left out", out$missing_or_nonfinite))
  if (out$share_capacitive < 0.5) flags = c(flags, "most points have Z'' > 0: check the sign convention of the imaginary column")
  if (!is.null(out$max_mod_diff_pct) && out$max_mod_diff_pct > 1) flags = c(flags, "|Z| in the file does not match Z' and Z''")
  if (!is.null(out$max_phase_diff_deg) && out$max_phase_diff_deg > 1) flags = c(flags, "the phase in the file does not match Z' and Z''")
  out$flags = if (length(flags)) paste(flags, collapse = "; ") else "none"
  out$note = e$note
  out
}

# Features of the spectrum that point to circuit families, read before any fit
eis_features = function(d) {
  f = d$freq_Hz; Z = d$Z; ph = Arg(Z) * 180 / pi; o = order(f); f = f[o]; Z = Z[o]; ph = ph[o]
  low = f <= 10 * min(f); high = f >= max(f) / 10
  slope_low = unname(stats::coef(stats::lm(log10(Mod(Z[low])) ~ log10(f[low])))[2])
  mz = -Im(Z); peaks = which(diff(sign(diff(stats::filter(mz, rep(1 / 3, 3), sides = 2)))) == -2) + 1
  out = data.frame(high_f_real_ohm = Re(Z[which.max(f)]), low_f_phase_deg = ph[1], low_f_logZ_slope = slope_low,
                   min_phase_deg = min(ph), f_at_min_phase_Hz = f[which.min(ph)], local_maxima_of_minus_Zim = length(peaks))
  hints = c(sprintf("Series resistance about %.3g ohm (real part at the highest frequency).", out$high_f_real_ohm))
  hints = c(hints, if (slope_low < -0.85 && ph[1] < -75) "Low-frequency |Z| slope near -1 and phase near -90 degrees: a blocking (capacitive) interface. Try Rs-C or Rs-Q. An Rct too large to resolve looks the same."
    else if (slope_low > -0.65 && slope_low < -0.35 && ph[1] < -30 && ph[1] > -60) "Low-frequency |Z| slope near -0.5 and phase near -45 degrees: diffusion (Warburg) or a strongly non-ideal interface (CPE n near 0.5)."
    else if (slope_low > -0.2) "Low-frequency |Z| levels off: a finite resistance in parallel with the interface (charge transfer or a transmissive finite diffusion)."
    else "Low-frequency behaviour between the ideal cases: compare circuits with a CPE and with diffusion elements.")
  if (out$local_maxima_of_minus_Zim >= 2) hints = c(hints, "-Z'' has more than one maximum: more than one time constant (two RC/RQ elements or nested circuits).")
  list(features = out, hints = hints)
}

# Effective capacitance of a CPE, only where a resistance is in parallel with it:
# Hsu and Mansfeld, C = (Q R^(1-n))^(1/n) with R the parallel resistance (the form in the EIS
# requirements; Zahner documentation eq. 7), and Brug et al. (J. Electroanal. Chem. 1984, 176, 275),
# C = Q^(1/n) (Rs Rp / (Rs + Rp))^((1-n)/n) (Zahner eq. 8)
hsu_mansfeld_capacitance = function(Q, n, R) (Q * R^(1 - n))^(1 / n)
brug_capacitance = function(Q, n, Rs, Rp) Q^(1 / n) * (Rs * Rp / (Rs + Rp))^((1 - n) / n)
# k0 = R T / (n^2 F^2 A Rct C), equal concentrations of both forms of the probe, small amplitude at the
# formal potential (Microchim. Acta 2023, 190, 251, eq. 17)
k0_from_rct = function(Rct_ohm, area_cm2, C_mM, n = 1, T_K = 298.15) R_const * T_K / (n^2 * F_const^2 * area_cm2 * Rct_ohm * mM_to_mol_cm3(C_mM))
# D from the Warburg coefficient with equal D and C of both forms:
# sigma = R T / (n^2 F^2 A 2^1/2) (1/(C_O D_O^1/2) + 1/(C_R D_R^1/2)) (Bard and Faulkner)
D_from_warburg = function(sigma, area_cm2, C_mM, n = 1, T_K = 298.15)
  (2 * R_const * T_K / (n^2 * F_const^2 * area_cm2 * sqrt(2) * mM_to_mol_cm3(C_mM) * sigma))^2

# Quantities that follow from a fit, each computed only when its circuit supports it
eis_derived = function(fit, area_cm2 = NA, C_mM = NA, n_e = 1, T_K = 298.15) {
  cc = eis_compile(fit$circuit); p = stats::setNames(fit$params$value, fit$params$parameter); out = list()
  series_R = if (cc$tree$type == "series") Filter(function(k) k$type == "element" && k$kind == "R", cc$tree$kids) else list()
  Rs = if (length(series_R)) p[[series_R[[1]]$name]] else NA
  # every parallel pair of one R and one Q gives an effective capacitance
  walk = function(t) if (t$type == "element") list() else c(
    if (t$type == "parallel" && length(t$kids) == 2 && all(vapply(t$kids, function(k) k$type == "element", TRUE))) {
      kinds = vapply(t$kids, `[[`, "", "kind")
      if (setequal(kinds, c("R", "Q"))) list(list(R = t$kids[[which(kinds == "R")]], Q = t$kids[[which(kinds == "Q")]])) }, do.call(c, lapply(t$kids, walk)))
  for (pr in walk(cc$tree)) {
    Q = p[[pr$Q$name]]; nn = p[[paste0("n", pr$Q$label)]]; R = p[[pr$R$name]]
    out[[paste0("C_eff_HsuMansfeld_", pr$Q$name, "_F")]] = hsu_mansfeld_capacitance(Q, nn, R)
    if (!is.na(Rs)) out[[paste0("C_eff_Brug_", pr$Q$name, "_F")]] = brug_capacitance(Q, nn, Rs, R)
  }
  for (k in names(p)) if (grepl("^Rct", k) || (k %in% c("R1", "R2") && !k %in% if (!is.na(Rs)) series_R[[1]]$name)) {
    if (!is.na(area_cm2)) out[[paste0(k, "_area_ohm_cm2")]] = p[[k]] * area_cm2
    if (grepl("^Rct", k) && !is.na(area_cm2) && !is.na(C_mM)) out[[paste0("k0_from_", k, "_cm_s")]] = k0_from_rct(p[[k]], area_cm2, C_mM, n_e, T_K)
  }
  for (k in grep("^sigma", names(p), value = TRUE)) if (!is.na(area_cm2) && !is.na(C_mM)) out[[paste0("D_from_", k, "_cm2_s")]] = D_from_warburg(p[[k]], area_cm2, C_mM, n_e, T_K)
  out
}

# A cautious reading of one fit: what is supported by the numbers and what is not
eis_interpret = function(fit, kk = NULL, feats = NULL) {
  p = fit$params; v = stats::setNames(p$value, p$parameter); s = stats::setNames(p$rel_se_pct, p$parameter); txt = character()
  if (!is.null(kk)) txt = c(txt, sprintf("Kramers-Kronig check: largest residual %.1f%% of |Z|, RMS %.2f%%. %s", kk$max_abs_res_pct, kk$rms_res_pct,
    if (kk$max_abs_res_pct < 2) "The data are consistent with a linear, stable system." else if (kk$max_abs_res_pct < 5) "Small departures. Look at where the residuals sit in frequency." else "Large departures: drift or nonlinearity is likely, and fitted values should not be trusted."))
  if ("Rs" %in% names(v)) txt = c(txt, sprintf("Series (electrolyte and contact) resistance Rs = %.3g ohm (error %.1f%%).", v[["Rs"]], s[["Rs"]]))
  for (k in grep("^n", names(v), value = TRUE)) txt = c(txt, sprintf("CPE exponent %s = %.2f: %s", k, v[[k]],
    if (v[[k]] > 0.9) "close to an ideal capacitor." else if (v[[k]] > 0.7) "a non-ideal capacitive interface (roughness, inhomogeneity)." else if (v[[k]] > 0.4) "far from a capacitor, as for diffusion or a distributed process." else "resistive behaviour. The element may not be a capacitance at all."))
  for (k in grep("^R(ct|[0-9])", names(v), value = TRUE)) txt = c(txt, if (is.finite(s[[k]]) && s[[k]] < 30) sprintf("%s = %.3g ohm is determined by the data (error %.0f%%).", k, v[[k]], s[[k]])
    else sprintf("%s = %.3g ohm is not well determined (error %s): do not report it as a charge-transfer resistance.", k, v[[k]], if (is.finite(s[[k]])) sprintf("%.0f%%", s[[k]]) else "not estimable"))
  for (k in grep("^sigma", names(v), value = TRUE)) txt = c(txt, sprintf("Warburg coefficient %s = %.3g ohm s^-1/2 (error %.0f%%): %s", k, v[[k]], s[[k]],
    if (s[[k]] < 30) "a diffusion contribution is supported by this fit." else "the diffusion term is weakly determined."))
  for (k in grep("^R_W[os]", names(v), value = TRUE)) { tk = sub("^R_", "tau_", k)
    txt = c(txt, sprintf("Finite-length diffusion %s: R = %.3g ohm (error %.0f%%), tau = %.3g s (error %.0f%%). %s", sub("^R_", "", k), v[[k]], s[[k]], v[[tk]], s[[tk]],
      if (startsWith(sub("^R_", "", k), "Wo")) "A reflective (open) boundary: the low-frequency branch turns capacitive, as for a layer the species cannot leave."
      else "A transmissive (short) boundary: the low-frequency branch turns resistive, as for diffusion through a layer of finite thickness.")) }
  if (any(p$at_bound)) txt = c(txt, sprintf("At a bound: %s. The optimum lies outside the allowed range or the parameter is not identifiable.", paste(p$parameter[p$at_bound], collapse = ", ")))
  if (!fit$converged) txt = c(txt, "The optimizer did not report convergence: treat the values with caution.")
  txt = c(txt, "A better fit statistic does not prove a circuit is the physically correct one. Keep the circuit whose elements make sense for the system.")
  txt
}

# Layout of a circuit diagram: boxes for elements, lines for wires (drawn by the app with ggplot2)
eis_diagram = function(circuit) {
  tree = eis_parse(circuit)
  size = function(t) switch(t$type, element = c(1.4, 1),
    series = { s = sapply(t$kids, size); c(sum(s[1, ]) + 0.4 * (ncol(s) - 1), max(s[2, ])) },
    parallel = { s = sapply(t$kids, size); c(max(s[1, ]) + 0.8, sum(s[2, ]) + 0.3 * (ncol(s) - 1)) })
  boxes = list(); wires = list()
  draw = function(t, x0, yc) { sz = size(t)
    if (t$type == "element") { boxes[[length(boxes) + 1]] <<- data.frame(xmin = x0 + 0.1, xmax = x0 + sz[1] - 0.1, ymin = yc - 0.3, ymax = yc + 0.3, label = t$name)
      wires[[length(wires) + 1]] <<- data.frame(x = c(x0, x0 + sz[1] - 0.1), xend = c(x0 + 0.1, x0 + sz[1]), y = yc, yend = yc); return(invisible()) }
    if (t$type == "series") { x = x0
      for (i in seq_along(t$kids)) { k = t$kids[[i]]; s = size(k); draw(k, x, yc)
        if (i < length(t$kids)) wires[[length(wires) + 1]] <<- data.frame(x = x + s[1], xend = x + s[1] + 0.4, y = yc, yend = yc)
        x = x + s[1] + 0.4 }
      return(invisible()) }
    s = sapply(t$kids, size); tot = sum(s[2, ]) + 0.3 * (ncol(s) - 1); y = yc + tot / 2; ys = numeric()
    for (i in seq_along(t$kids)) { yk = y - s[2, i] / 2; ys = c(ys, yk)
      draw(t$kids[[i]], x0 + 0.4 + (sz[1] - 0.8 - s[1, i]) / 2, yk)
      wires[[length(wires) + 1]] <<- data.frame(x = c(x0 + 0.4, x0 + 0.4 + (sz[1] - 0.8 + s[1, i]) / 2), xend = c(x0 + 0.4 + (sz[1] - 0.8 - s[1, i]) / 2, x0 + sz[1] - 0.4), y = yk, yend = yk)
      y = y - s[2, i] - 0.3 }
    wires[[length(wires) + 1]] <<- data.frame(x = c(x0, x0 + 0.4, x0 + sz[1] - 0.4, x0 + sz[1] - 0.4), xend = c(x0 + 0.4, x0 + 0.4, x0 + sz[1] - 0.4, x0 + sz[1]),
                                              y = c(yc, max(ys), max(ys), yc), yend = c(yc, min(ys), min(ys), yc))
  }
  sz = size(tree); draw(tree, 0.6, 0)
  wires[[length(wires) + 1]] = data.frame(x = c(0, sz[1] + 0.6), xend = c(0.6, sz[1] + 1.2), y = 0, yend = 0)
  list(boxes = do.call(rbind, boxes), wires = do.call(rbind, wires), width = sz[1] + 1.2, height = sz[2])
}

# Display units: values are divided only for plots and tables; fitting always uses ohm, Hz and F
unit_factor = function(u) c(ohm = 1, kohm = 1e3, Mohm = 1e6, Gohm = 1e9, Hz = 1, kHz = 1e3, MHz = 1e6,
                            F = 1, mF = 1e-3, uF = 1e-6, nF = 1e-9, pF = 1e-12)[[u]]
unit_label = function(u) c(ohm = "Ω", kohm = "kΩ", Mohm = "MΩ", Gohm = "GΩ", Hz = "Hz", kHz = "kHz", MHz = "MHz",
                           F = "F", mF = "mF", uF = "µF", nF = "nF", pF = "pF")[[u]]
auto_ohm = function(z) { m = max(abs(c(Re(z), Im(z)))); if (m >= 1e9) "Gohm" else if (m >= 1e6) "Mohm" else if (m >= 1e3) "kohm" else "ohm" }

# Impedance files: a PalmSens session, or any text or CSV export with a frequency column and the real and
# imaginary parts (CHI "Freq/Hz, Z'/ohm, Z"/ohm, Z/ohm, Phase/deg", or names such as Zre, Re(Z), Zim,
# -Im(Z), -Z''). A leading minus in the name means the column holds -Z''. When no column says so and most
# values are positive, the column is taken as -Z'' and a note is returned. |Z| and phase columns, when
# present, are kept for the validation check.
read_eis = function(path, name = basename(path), instrument = "auto") {
  fmt = instrument_format(path, name, instrument)
  if (fmt == "pssession") {
    r = read_pssession(path, name)
    if (r$kind != "eis") stop(name, " is not an impedance measurement")
    return(list(file = name, data = r$data, note = "", dropped = 0))
  }
  if (fmt %in% c("biologic", "gamry", "mpr")) {
    r = switch(fmt, biologic = read_biologic_mpt(path, name), gamry = read_gamry_dta(path, name), mpr = read_biologic_mpr(path, name))
    if (r$kind != "eis") stop(name, " is not an impedance measurement (", r$params$technique, ")")
    d = r$data[order(-r$data$freq_Hz), ]
    return(list(file = name, data = d, note = "", dropped = 0))
  }
  L = readLines(path, warn = FALSE, encoding = "UTF-8")
  # CHI "data only" export: frequency, Z', Z'', |Z|, phase with no header
  if (fmt == "chi" && !any(grepl("freq", L[1:min(60, length(L))], ignore.case = TRUE))) {
    r = read_chi_txt(path, name)
    if (!identical(r$kind, "eis")) stop(name, " is not an impedance measurement (", r$params$technique, ")")
    d = r$data[is.finite(r$data$freq_Hz) & r$data$freq_Hz > 0, ]
    return(list(file = name, data = d[order(-d$freq_Hz), ], note = "", dropped = 0))
  }
  hdr = which(grepl("freq", L, ignore.case = TRUE) & grepl("[,;\t]", L))
  if (!length(hdr)) stop(name, ": no header line with a frequency column was found")
  h = hdr[length(hdr)]
  sep = if (grepl("\t", L[h])) "\t" else if (grepl(";", L[h])) ";" else ","
  tab = utils::read.table(text = L[h:length(L)], sep = sep, header = TRUE, check.names = FALSE, fill = TRUE,
                          strip.white = TRUE, comment.char = "", quote = "", stringsAsFactors = FALSE)
  # names written by CSV writers come in quotes; strip them so a quote is not read as Z"
  nm = gsub('^"|"$', "", trimws(names(tab))); low = tolower(gsub("[[:space:]]", "", nm))
  fcol = which(grepl("freq", low))[1]
  icol = which(grepl("z''|z\"|zim|im\\(z\\)|zimag|imag", low))[1]
  rcol = which(grepl("z'|zre|re\\(z\\)|zreal|real", low) & !seq_along(low) %in% icol)[1]
  mcol = which(grepl("^\\|?z\\|?(/|$|\\()|^mod|^\\|z\\|", low) & !seq_along(low) %in% c(rcol, icol))[1]
  pcol = which(grepl("phase", low))[1]
  if (any(is.na(c(fcol, rcol, icol)))) stop(name, ": columns for frequency, Z' and Z'' were not all found (columns: ", paste(nm, collapse = ", "), ")")
  num = function(v) suppressWarnings(as.numeric(v))
  f = num(tab[[fcol]]); zr = num(tab[[rcol]]); zi = num(tab[[icol]])
  zm = if (!is.na(mcol)) num(tab[[mcol]]) else NULL; ph = if (!is.na(pcol)) num(tab[[pcol]]) else NULL
  k = is.finite(f) & is.finite(zr) & is.finite(zi) & f > 0
  dropped = sum(!k & rowSums(!is.na(tab)) > 0)
  f = f[k]; zr = zr[k]; zi = zi[k]; zm = zm[k]; ph = ph[k]
  note = ""
  if (grepl("^-", low[icol])) zi = -zi
  else if (mean(zi > 0) > 0.5) { zi = -zi; note = "imaginary column read as -Z'' (most values were positive)" }
  o = order(-f)
  list(file = name, data = data.frame(freq_Hz = f[o], Z = complex(real = zr[o], imaginary = zi[o])), note = note, dropped = dropped,
       file_mod = if (length(zm)) zm[o] else NULL, file_phase = if (length(ph)) ph[o] else NULL)
}
