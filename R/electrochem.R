# Electrode geometry, surface-area references and the standard equations of each
# technique (voltammetry, chronoamperometry, chronocoulometry, hydrodynamic electrodes,
# impedance). Every function returns plain numbers or data frames in SI-style units
# (A, V, s, C, cm, cm2, mol/cm3 inside the equations), so the app, the report and the
# tests use the same numbers. Sources are given next to each equation and listed in
# equation_sources() for the Equations tab.

R_const = 8.314462618
if (!exists("F_const")) F_const = 96485.33212
mM_to_mol_cm3 = function(c_mM) c_mM * 1e-6

# Electrode geometry
# Disk sizes sold by the main makers, checked on the makers' pages (October 2026).
# The list only helps to fill the diameter; any other size can be typed.
electrode_sizes = data.frame(
  label = c("2 mm disk (CH Instruments CHI101 Au, CHI102 Pt, CHI103 Ag)",
            "3 mm disk (CH Instruments CHI104 glassy carbon, BASi MF-2012 glassy carbon)",
            "1.6 mm disk (BASi MF-2013 Pt, MF-2014 Au)",
            "3 mm printed carbon (CH Instruments SE101)",
            "4 mm screen-printed carbon (Metrohm DropSens 110, maker gives 0.11 cm²)",
            "5 mm disk (BioLogic glassy carbon, Pine RRDE disk)",
            "4.57 mm disk (Pine RRDE)", "5.5 mm disk (Pine RRDE)", "5.61 mm disk (Pine RRDE)",
            "25 um microdisk (CH Instruments CHI106 Au, CHI108 Pt)",
            "12.5 um microdisk (CH Instruments CHI105 Au)",
            "10 um microdisk (CH Instruments CHI107 Pt)"),
  d_mm = c(2, 3, 1.6, 3, 4, 5, 4.57, 5.5, 5.61, 0.025, 0.0125, 0.010))

# Geometric area (cm2) of common electrode shapes. Lengths in mm.
# rectangle: a flat plate or strip; faces = 1 for one exposed face (the rest masked),
#   2 for both faces; edges = TRUE adds the immersed edges of a strip of the given
#   thickness (two long edges and the bottom edge: 2 L t + W t)
# wire: a cylinder of diameter d immersed to length L, with or without its flat end
# band: a band microelectrode of length L and width W (one face)
electrode_area = function(shape = c("disk", "rectangle", "wire", "band", "hemisphere", "sphere", "area"),
                          d_mm = NA, length_mm = NA, width_mm = NA, faces = 1, edges = FALSE,
                          thickness_mm = 0, wire_end = FALSE, area_cm2 = NA) {
  shape = match.arg(shape)
  mm2 = switch(shape,
    disk = pi * (d_mm / 2)^2,
    rectangle = faces * length_mm * width_mm + if (isTRUE(edges)) 2 * length_mm * thickness_mm + width_mm * thickness_mm else 0,
    wire = pi * d_mm * length_mm + if (isTRUE(wire_end)) pi * (d_mm / 2)^2 else 0,
    band = length_mm * width_mm,
    hemisphere = 2 * pi * (d_mm / 2)^2,
    sphere = 4 * pi * (d_mm / 2)^2,
    area = area_cm2 * 100)
  mm2 / 100
}

# Surface-area reference charges
# Charge per cm2 of real surface for one monolayer, from Lukaszewski, Soszko and
# Czerwinski, Int. J. Electrochem. Sci. 11 (2016) 4442 (review, values as listed there;
# hydrogen and oxide values follow Trasatti and Petrii, Pure Appl. Chem. 63 (1991) 711).
# For PdO the review lists 424. Any other value (for example 405, which has no published source here) is typed in
# by the user, and the value used is reported.
ecsa_references = data.frame(
  key = c("pt_h", "rh_h", "ir_h", "pd_ox424", "au_ox", "pt_ox", "rh_ox", "pt_co", "cu_upd", "ag_upd", "pb_upd"),
  label = c("Pt, hydrogen desorption, 210 µC/cm²", "Rh, hydrogen desorption, 221 µC/cm²", "Ir, hydrogen desorption, 218 µC/cm²",
            "Pd, PdO reduction, 424 µC/cm²", "Au, oxide reduction, 390 µC/cm²",
            "Pt, oxide reduction, 420 µC/cm²", "Rh, oxide reduction, 442 µC/cm²", "Pt, CO stripping, 420 µC/cm²",
            "Cu underpotential deposition (stripping), 420 µC/cm²", "Ag underpotential deposition (stripping), 210 µC/cm²",
            "Pb underpotential deposition (stripping), 300 µC/cm²"),
  q_uC_cm2 = c(210, 221, 218, 424, 390, 420, 442, 420, 420, 210, 300),
  sweep = c("anodic", "anodic", "anodic", "cathodic", "cathodic", "cathodic", "cathodic", "anodic", "anodic", "anodic", "anodic"),
  source = c(rep("Lukaszewski et al. 2016 (Trasatti and Petrii 1991)", 3), "Lukaszewski et al. 2016 (Woods 1976)",
             "Lukaszewski et al. 2016 (390 ± 10)",
             rep("Lukaszewski et al. 2016", 6)))

# Double-layer capacitance from CVs in a potential range without faradaic current.
# At each scan rate the charging current is half the gap between the anodic and the
# cathodic sweep at E_at; its slope against scan rate is C_dl. ECSA = C_dl / C_s, where
# C_s (the capacitance of 1 cm2 of the same surface) must be assumed: reported values
# for metals lie mostly between 20 and 40 uF/cm2 and the error can reach 100%
# (Lukaszewski et al. 2016). Best used to compare electrodes of the same material.
cdl_from_cvs = function(xs, E_at, Cs_uF_cm2 = 40) {
  pts = do.call(rbind, lapply(xs, function(x) {
    a = last_sweep(x, "anodic", 0.02); c_ = last_sweep(x, "cathodic", 0.02)
    ia = stats::approx(a$E_V, a$I_A, E_at, ties = mean)$y; ic = stats::approx(c_$E_V, c_$I_A, E_at, ties = mean)$y
    data.frame(file = x$file, v_V_s = x$params$scan_rate_V_s, I_anodic_A = ia, I_cathodic_A = ic, I_charging_A = (ia - ic) / 2)
  }))
  pts = pts[order(pts$v_V_s), ]
  fit = if (nrow(pts) >= 2) { m = stats::lm(I_charging_A ~ v_V_s, pts)
    data.frame(n = nrow(pts), E_at_V = E_at, Cdl_uF = unname(stats::coef(m)[2]) * 1e6, intercept_A = unname(stats::coef(m)[1]),
               r2 = if (nrow(pts) > 2) summary(m)$r.squared else NA, Cs_uF_cm2 = Cs_uF_cm2,
               ECSA_cm2 = unname(stats::coef(m)[2]) * 1e6 / Cs_uF_cm2) } else NULL
  list(points = pts, fit = fit)
}

# Cyclic voltammetry: peaks and diagnostics
# Peak currents are measured from zero current (no baseline). For a second peak on a
# decaying background, or a small reverse peak, the reverse peak current read this way
# is only approximate.
sweep_peak = function(s, dir) {
  i = if (dir == "anodic") which.max(s$I_A) else which.min(s$I_A)
  ip = s$I_A[i]
  # half-peak potential: on the rising side, where the current first reaches half the peak
  before = seq_len(i)
  h = before[if (dir == "anodic") s$I_A[before] >= ip / 2 else s$I_A[before] <= ip / 2][1]
  data.frame(Ep_V = s$E_V[i], ip_A = ip, Ep_half_V = s$E_V[h])
}
cv_peaks = function(x, window = NULL, T_K = 298.15) {
  one = function(dir) {
    s = tryCatch(last_sweep(x, dir), error = function(e) NULL)
    if (is.null(s)) return(data.frame(Ep_V = NA, ip_A = NA, Ep_half_V = NA))
    if (!is.null(window)) s = s[s$E_V >= min(window) & s$E_V <= max(window), ]
    if (nrow(s) < 3) return(data.frame(Ep_V = NA, ip_A = NA, Ep_half_V = NA))
    sweep_peak(s, dir)
  }
  a = one("anodic"); c_ = one("cathodic")
  dEp = a$Ep_V - c_$Ep_V
  data.frame(file = x$file, v_V_s = x$params$scan_rate_V_s, Epa_V = a$Ep_V, ipa_A = a$ip_A, Epc_V = c_$Ep_V, ipc_A = c_$ip_A,
             dEp_mV = 1000 * dEp, E_half_V = (a$Ep_V + c_$Ep_V) / 2, ratio_ipa_ipc = abs(a$ip_A / c_$ip_A),
             Epa_minus_Ep2_mV = 1000 * abs(a$Ep_V - a$Ep_half_V), Epc_minus_Ep2_mV = 1000 * abs(c_$Ep_V - c_$Ep_half_V))
}
# Transfer coefficient from the peak shape of an irreversible wave:
# |Ep - Ep/2| = 1.857 RT / (alpha n F), 47.7/alpha mV at 25 C (Microchim. Acta 2023, 190, 251, eq. 9)
alpha_from_peak_width = function(width_mV, n = 1, T_K = 298.15) 1.857 * R_const * T_K / (n * F_const * width_mV / 1000)

# Randles-Sevcik. ip = 0.4463 n F A C (n F v D / R T)^1/2 for a reversible wave, which is
# 2.69e5 n^3/2 A D^1/2 C v^1/2 at 25 C; ip = 0.4958 n F A C (alpha n F v D / R T)^1/2 for an
# irreversible wave, 2.99e5 n alpha^1/2 A D^1/2 C v^1/2 at 25 C with n' = n = 1
# (Bard and Faulkner; Microchim. Acta 2023, 190, 251, eqs. 3 and 7).
# Units: A in cm2, D in cm2/s, C in mol/cm3, v in V/s, ip in A.
rs_coef = function(n = 1, T_K = 298.15, reversible = TRUE, alpha = 0.5)
  if (reversible) 0.4463 * n * F_const * sqrt(n * F_const / (R_const * T_K)) else 0.4958 * n * F_const * sqrt(alpha * n * F_const / (R_const * T_K))
# slope of |ip| against v^1/2 gives the area when D is known, or D when the area is known
rs_solve = function(slope_A_per_sqrtV_s, n = 1, C_mM, D_cm2_s = NA, area_cm2 = NA, T_K = 298.15, reversible = TRUE, alpha = 0.5) {
  k = rs_coef(n, T_K, reversible, alpha) * mM_to_mol_cm3(C_mM)
  data.frame(area_cm2 = if (!is.na(D_cm2_s)) slope_A_per_sqrtV_s / (k * sqrt(D_cm2_s)) else NA,
             D_cm2_s = if (!is.na(area_cm2)) (slope_A_per_sqrtV_s / (k * area_cm2))^2 else NA)
}

# Standard rate constant k0 (cm/s) from the peak separation of one CV.
# Nicholson, with the Lavagnini et al. (Electroanalysis 2004, 16, 505) fit of the working
# curve: psi = (-0.6288 + 0.0021 X) / (1 - 0.017 X), X = n dEp in mV, valid up to about
# 200 mV; k0 = psi (pi D n v F / R T)^1/2.
# Klingler and Kochi (J. Phys. Chem. 1981, 85, 1731) for n dEp above 150 mV:
# k0 = 2.18 (alpha D n v F / R T)^1/2 exp(-(alpha^2 n F / R T) dEp).
# Both as given in Microchim. Acta 2023, 190, 251 (eqs. 11 to 13, Table 1).
k0_nicholson = function(dEp_mV, v_V_s, D_cm2_s, n = 1, T_K = 298.15) {
  X = n * dEp_mV
  psi = (-0.6288 + 0.0021 * X) / (1 - 0.017 * X)
  k0 = psi * sqrt(pi * D_cm2_s * n * v_V_s * F_const / (R_const * T_K))
  ifelse(X > 61 & X <= 200, k0, NA)
}
k0_klingler_kochi = function(dEp_mV, v_V_s, D_cm2_s, alpha = 0.5, n = 1, T_K = 298.15) {
  f = n * F_const / (R_const * T_K)
  k0 = 2.18 * sqrt(alpha * D_cm2_s * n * v_V_s * F_const / (R_const * T_K)) * exp(-alpha^2 * f * dEp_mV / 1000)
  ifelse(n * dEp_mV > 150, k0, NA)
}
# Which regime the peak separation points to (Microchim. Acta 2023: reversible about
# 57 to 63 mV, quasi-reversible 63 to 200 mV, irreversible above 200 mV for n dEp)
cv_regime = function(dEp_mV, n = 1) {
  X = n * dEp_mV
  ifelse(is.na(X), NA, ifelse(X <= 63, "reversible", ifelse(X <= 200, "quasi-reversible", "irreversible")))
}

# Quasi-reversible Randles-Sevcik: ip = ip(reversible) x K(Lambda, alpha) (Matsuda and Ayabe; used for the
# electroactive area by Trachioti et al., Microchim. Acta 2023, 190, 251, eq. 4). The paper reads K off published
# graphs. Here K, the peak width Delta = (Ep/2 - Ep) F/RT and the peak separation come from a table made with a
# Butler-Volmer finite-difference simulation (one form in the bulk, equal D, first sweep, peaks from zero), checked
# against the reversible and irreversible Randles-Sevcik peaks (0.4% low, cancelled by taking K against the
# reversible run on the same grid) and Nicholson's dEp-psi values (within 2 mV). See rs_quasi_table in
# R/rs_quasi_table.R and scripts/12_rs_quasi_table.R. Lambda = k0 / sqrt(D f v) = psi sqrt(pi).
# Steps as in the paper: Lambda from n dEp, alpha from the width of the forward peak, then K(Lambda, alpha).
rs_quasi_K = function(dEp_mV, width_mV, n = 1, T_K = 298.15, tab = rs_quasi_table) {
  u = 1000 * R_const * T_K / F_const                        # RT/F in mV
  X = n * dEp_mV / u; W = n * width_mV / u
  alphas = sort(unique(tab$alpha))
  one = function(X, W) {
    if (!is.finite(X) || !is.finite(W)) return(c(NA, NA, NA, NA))
    # for each alpha in the table, the log Lambda that gives this dEp, and the width and K there
    at = t(vapply(alphas, function(a) { s = tab[tab$alpha == a, ]; s = s[order(s$dEp_RT_F), ]
      if (X < min(s$dEp_RT_F) || X > max(s$dEp_RT_F)) return(c(NA, NA, NA))
      lL = stats::approx(s$dEp_RT_F, s$log10_Lambda, X)$y
      s = s[order(s$log10_Lambda), ]
      c(lL, stats::approx(s$log10_Lambda, s$Delta, lL)$y, stats::approx(s$log10_Lambda, s$K, lL)$y) }, numeric(3)))
    ok = is.finite(at[, 1]); if (sum(ok) < 2) return(c(NA, NA, NA, NA))
    a = alphas[ok]; at = at[ok, , drop = FALSE]
    # width falls as alpha rises, so alpha comes from the measured width; outside the table it is held at the edge
    o = order(at[, 2]); a_hat = stats::approx(at[o, 2], a[o], W, rule = 2)$y
    c(a_hat, stats::approx(a, at[, 1], a_hat)$y, stats::approx(a, at[, 3], a_hat)$y,
      as.numeric(W < min(at[, 2]) || W > max(at[, 2])))
  }
  r = t(mapply(one, X, W))
  data.frame(alpha = r[, 1], log10_Lambda = r[, 2], psi = 10^r[, 2] / sqrt(pi), K = r[, 3], alpha_at_table_edge = r[, 4] == 1)
}

# Surface-confined (adsorbed) species: surface coverage from the peak charge,
# Gamma = Q / (n F A) (mol/cm2), and the peak current of an ideal reversible surface wave,
# ip = n^2 F^2 v A Gamma / (4 R T) (Bard and Faulkner), used as a check.
surface_coverage = function(Q_C, n = 1, area_cm2) Q_C / (n * F_const * area_cm2)
ip_surface = function(Gamma_mol_cm2, v_V_s, n = 1, area_cm2, T_K = 298.15) n^2 * F_const^2 * v_V_s * area_cm2 * Gamma_mol_cm2 / (4 * R_const * T_K)
# Laviron plot at high scan rates: Epc against log10 v has slope -2.303 RT/(alpha n F) and
# Epa slope 2.303 RT/((1 - alpha) n F) (Laviron, J. Electroanal. Chem. 1979, 101, 19, eqs. 20 and 21). ks from
# eq. 22 of the same paper, valid for n dEp above 200 mV (checked against the paper; it writes 2.3 for ln 10):
# log ks = alpha log(1 - alpha) + (1 - alpha) log alpha - log(RT / n F v) - alpha (1 - alpha) n F dEp / (2.303 RT)
laviron = function(v_V_s, Epa_V, Epc_V, n = 1, T_K = 298.15, min_dEp_mV = 200 / n) {
  k = abs(Epa_V - Epc_V) * 1000 >= min_dEp_mV
  if (sum(k) < 2) return(NULL)
  lv = log10(v_V_s[k]); sc = unname(stats::coef(stats::lm(Epc_V[k] ~ lv))[2]); sa = unname(stats::coef(stats::lm(Epa_V[k] ~ lv))[2])
  g = 2.303 * R_const * T_K / (n * F_const)
  alpha = -g / sc                                     # from the cathodic slope
  alpha_a = 1 - g / sa                                # from the anodic slope, for comparison
  vmax = max(v_V_s[k]); dE = abs(Epa_V - Epc_V)[k][which.max(v_V_s[k])]
  logks = alpha * log10(1 - alpha) + (1 - alpha) * log10(alpha) - log10(R_const * T_K / (n * F_const * vmax)) -
          alpha * (1 - alpha) * n * F_const * dE / (2.303 * R_const * T_K)
  data.frame(points = sum(k), slope_cathodic_V = sc, slope_anodic_V = sa, alpha_cathodic = alpha, alpha_anodic = alpha_a,
             ks_per_s = if (alpha > 0 && alpha < 1) 10^logks else NA)
}

# Chronoamperometry
# Cottrell: i = n F A C D^1/2 / (pi t)^1/2 (Bard and Faulkner). The slope of i against
# t^-1/2 over the chosen time window gives D (area known) or the area (D known).
cottrell_fit = function(x, t_from, t_to, n = 1, C_mM, D_cm2_s = NA, area_cm2 = NA) {
  d = x$data; k = d$t_s >= t_from & d$t_s <= t_to & d$t_s > 0
  if (sum(k) < 5) return(NULL)
  m = stats::lm(abs(d$I_A[k]) ~ I(1 / sqrt(d$t_s[k])))
  s = unname(stats::coef(m)[2]); C = mM_to_mol_cm3(C_mM)
  data.frame(file = x$file, from_s = t_from, to_s = t_to, points = sum(k), slope_A_s05 = s, intercept_A = unname(stats::coef(m)[1]),
             r2 = summary(m)$r.squared,
             D_cm2_s = if (!is.na(area_cm2)) pi * (s / (n * F_const * area_cm2 * C))^2 else NA,
             area_cm2 = if (!is.na(D_cm2_s)) s * sqrt(pi) / (n * F_const * C * sqrt(D_cm2_s)) else NA)
}
# Steady-state current at a microelectrode: disk 4 n F D C r (Saito), hemisphere and
# sphere 2 pi n F D C r and 4 pi n F D C r (Bard and Faulkner). Solves for whichever of
# D, C or r is NA.
ume_factor = function(geometry = c("disk", "hemisphere", "sphere")) switch(match.arg(geometry), disk = 4, hemisphere = 2 * pi, sphere = 4 * pi)
ume_steady = function(i_ss_A, n = 1, D_cm2_s = NA, C_mM = NA, r_cm = NA, geometry = "disk") {
  g = ume_factor(geometry); i = abs(i_ss_A); C = mM_to_mol_cm3(C_mM)
  data.frame(i_ss_A = i,
             D_cm2_s = if (is.na(D_cm2_s)) i / (g * n * F_const * C * r_cm) else D_cm2_s,
             C_mM = if (is.na(C_mM)) i / (g * n * F_const * D_cm2_s * r_cm) * 1e6 else C_mM,
             r_um = if (is.na(r_cm)) i / (g * n * F_const * D_cm2_s * C) * 1e4 else r_cm * 1e4)
}

# Chronocoulometry (Anson plot)
# Q = 2 n F A C D^1/2 t^1/2 / pi^1/2 + Q_dl + n F A Gamma (Anson; Microchim. Acta 2023,
# eq. 1). The slope gives D or A; the intercept minus the blank intercept (Q_dl, from the
# same step in electrolyte alone) gives the adsorbed amount Gamma.
anson_fit = function(x, t_from, t_to, n = 1, C_mM, D_cm2_s = NA, area_cm2 = NA, Q_dl_C = NA) {
  d = x$data
  q = if (!is.null(d$Q_C)) d$Q_C else c(0, cumsum(diff(d$t_s) * (d$I_A[-1] + d$I_A[-nrow(d)]) / 2))
  k = d$t_s >= t_from & d$t_s <= t_to
  if (sum(k) < 5) return(NULL)
  m = stats::lm(abs(q[k]) ~ sqrt(d$t_s[k]))
  s = unname(stats::coef(m)[2]); b = unname(stats::coef(m)[1]); C = mM_to_mol_cm3(C_mM)
  A = if (!is.na(area_cm2)) area_cm2 else if (!is.na(D_cm2_s)) s * sqrt(pi) / (2 * n * F_const * C * sqrt(D_cm2_s)) else NA
  data.frame(file = x$file, from_s = t_from, to_s = t_to, points = sum(k), slope_C_per_s05 = s, intercept_C = b, r2 = summary(m)$r.squared,
             D_cm2_s = if (!is.na(area_cm2)) pi * (s / (2 * n * F_const * area_cm2 * C))^2 else NA,
             area_cm2 = A,
             Gamma_mol_cm2 = if (!is.na(Q_dl_C) && !is.na(A)) (b - Q_dl_C) / (n * F_const * A) else NA)
}

# Rotating disk electrode
rpm_to_omega = function(rpm) 2 * pi * rpm / 60
# limiting current against omega^1/2: the slope gives D (or n, when D is given)
levich_fit = function(rpm, i_L_A, n = 1, area_cm2, C_mM, nu_cm2_s = 0.01, D_cm2_s = NA) {
  w = sqrt(rpm_to_omega(rpm)); m = stats::lm(abs(i_L_A) ~ w); s = unname(stats::coef(m)[2])
  k = 0.620 * F_const * area_cm2 * nu_cm2_s^(-1 / 6) * mM_to_mol_cm3(C_mM)
  data.frame(rates = length(rpm), slope_A_per_rad05 = s, intercept_A = unname(stats::coef(m)[1]), r2 = summary(m)$r.squared,
             D_cm2_s = (s / (n * k))^(3 / 2), n_from_slope = if (!is.na(D_cm2_s)) s / (k * D_cm2_s^(2 / 3)) else NA)
}
# Koutecky-Levich: 1/i = 1/i_k + 1/(B omega^1/2). Intercept gives the kinetic current
# i_k, the slope 1/B gives n when D, C, A and nu are known (or D when n is known).
koutecky_levich = function(rpm, i_A, n = NA, area_cm2, C_mM, nu_cm2_s = 0.01, D_cm2_s = NA) {
  x = 1 / sqrt(rpm_to_omega(rpm)); y = 1 / abs(i_A); m = stats::lm(y ~ x)
  s = unname(stats::coef(m)[2]); b = unname(stats::coef(m)[1])
  k = 0.620 * F_const * area_cm2 * nu_cm2_s^(-1 / 6) * mM_to_mol_cm3(C_mM)
  data.frame(rates = length(rpm), slope = s, intercept = b, r2 = summary(m)$r.squared, i_k_A = 1 / b,
             n = if (!is.na(D_cm2_s)) 1 / (s * k * D_cm2_s^(2 / 3)) else NA,
             D_cm2_s = if (!is.na(n)) (1 / (s * k * n))^(3 / 2) else NA)
}
# RRDE, oxygen reduction: electron number and peroxide yield from the disk and ring
# currents and the collection efficiency N:
# n = 4 |I_D| / (|I_D| + |I_R| / N), H2O2 % = 200 (|I_R| / N) / (|I_D| + |I_R| / N)
rrde_n_peroxide = function(I_disk_A, I_ring_A, N) {
  d = abs(I_disk_A); r = abs(I_ring_A) / N
  data.frame(n = 4 * d / (d + r), H2O2_pct = 200 * r / (d + r))
}

# Tafel analysis
# Tafel slope b (V/decade) = 2.303 R T / (alpha n F), so alpha n = 2.303 R T / (b F).
# The exchange current is the Tafel line extrapolated to the equilibrium potential.
tafel_kinetics = function(slope_dec_per_V, intercept_log10_A, E_eq_V = NA, area_cm2 = NA, T_K = 298.15) {
  b = 1 / slope_dec_per_V                                   # the fit is log10 I against E
  i0 = if (!is.na(E_eq_V)) 10^(intercept_log10_A + slope_dec_per_V * E_eq_V) else NA
  data.frame(tafel_mV_dec = 1000 * abs(b), alpha_n = 2.303 * R_const * T_K / (abs(b) * F_const),
             i0_A = i0, j0_A_cm2 = i0 / area_cm2)
}

# Probe diffusion coefficients used in the literature for electroactive-area work
redox_probes = data.frame(
  key = c("ferri", "ruhex", "other"),
  label = c("Ferricyanide, 7.6e-6 cm2/s", "Hexaammineruthenium(III), 9.1e-6 cm2/s", "Other (enter D)"),
  D_cm2_s = c(7.6e-6, 9.1e-6, NA),
  source = c("Microchim. Acta 2023, 190, 251", "Microchim. Acta 2023, 190, 251", ""))

# Electrochemical impedance: see R/eis.R

# Equations used in the app, with their sources, for the Equations tab
equation_sources = function() data.frame(
  Technique = c("Surface area", "Surface area", "Surface area", "Electrode area",
                "CV", "CV", "CV", "CV", "CV", "CV", "CV", "CV", "CV",
                "Chronoamperometry", "Chronoamperometry", "Chronoamperometry", "Chronocoulometry",
                "Chronopotentiometry", "RDE", "RDE", "RRDE", "RRDE", "Tafel", "Tafel",
                "EIS", "EIS", "EIS", "EIS", "EIS", "EIS", "Calibration"),
  Quantity = c("ECSA from a monolayer charge", "ECSA from double-layer capacitance", "Electroactive area from a redox probe", "Geometric area of each shape",
               "Peak current, reversible (Randles–Ševčík)", "Peak current, irreversible", "Peak current, quasi-reversible", "Mixed solutions (both forms present)", "Transfer coefficient from peak width",
               "k⁰, Nicholson (Lavagnini fit)", "k⁰, Klingler-Kochi", "Surface coverage of an adsorbed species", "Laviron plot (adsorbed species)",
               "Cottrell", "Microelectrode steady state", "Microdisk transient", "Anson",
               "Sand", "Levich", "Koutecky-Levich", "Collection efficiency", "Electron number and peroxide yield", "Tafel slope", "Exchange current",
               "CPE and Warburg elements", "k⁰ from charge-transfer resistance", "Effective capacitance of a CPE (Hsu-Mansfeld)", "Effective capacitance of a CPE (Brug)", "Validity of the spectrum", "Model ranking", "Detection limit"),
  Equation = c("ECSA = Q / q_ref", "C_dl = slope of (i_a - i_c)/2 against v, then ECSA = C_dl / C_s",
               "slope of ip against v^1/2 with the Randles–Ševčík equation and a known D",
               "disk pi d^2/4, plate L W per face, strip edges 2 L t + W t, wire pi d L (+ pi d^2/4), hemisphere 2 pi r^2",
               "ip = 0.4463 n F A C (n F v D / R T)^1/2 = 2.69e5 n^3/2 A D^1/2 C v^1/2 at 25 C",
               "ip = 0.4958 n F A C (alpha n F v D / R T)^1/2 = 2.99e5 n alpha^1/2 A D^1/2 C v^1/2 at 25 C",
               "ip = ip(reversible) x K(Lambda, alpha), Lambda = k⁰ / (D F v / R T)^1/2 from n dEp, alpha from Ep - Ep/2. Area from the slope of ip/K against v^1/2",
               "Randles–Ševčík assumes one form in the bulk. With both forms, use Cottrell or Anson after a step from the rest potential",
               "|Ep - Ep/2| = 1.857 R T / (alpha n F) = 47.7/(alpha n) mV at 25 C",
               "psi = (-0.6288 + 0.0021 n dEp)/(1 - 0.017 n dEp), k0 = psi (pi D n v F / R T)^1/2, n dEp 61 to 200 mV",
               "k0 = 2.18 (alpha D n v F / R T)^1/2 exp(-alpha^2 n F dEp / R T), n dEp above 150 mV",
               "Gamma = Q / (n F A) and ip = n^2 F^2 v A Gamma / (4 R T)",
               "Epc slope against log v = -2.303 R T/(alpha n F). log ks = alpha log(1 - alpha) + (1 - alpha) log alpha - log(R T/(n F v)) - alpha (1 - alpha) n F dEp/(2.303 R T), n dEp above 200 mV",
               "i = n F A C D^1/2 / (pi t)^1/2",
               "disk i = 4 n F D C r, hemisphere 2 pi n F D C r, sphere 4 pi n F D C r",
               "i = 4 n F D C r f(tau), tau = 4 D t / r^2, f = 0.7854 + 0.8862 tau^-1/2 + 0.2146 exp(-0.7823 tau^-1/2)",
               "Q = 2 n F A C D^1/2 t^1/2 / pi^1/2 + Q_dl + n F A Gamma",
               "i tau^1/2 = n F A C (pi D)^1/2 / 2",
               "i_L = 0.620 n F A D^2/3 omega^1/2 nu^-1/6 C (omega in rad/s)",
               "1/i = 1/i_k + 1/(B omega^1/2), B = 0.620 n F A D^2/3 nu^-1/6 C",
               "N = -I_ring / I_disk on the limiting plateau",
               "n = 4 |I_D| / (|I_D| + |I_R|/N) and H2O2 % = 200 (|I_R|/N) / (|I_D| + |I_R|/N)",
               "b = 2.303 R T / (alpha n F)", "Tafel line extrapolated to the equilibrium potential",
               "Z_CPE = 1/(Q (j w)^n), Z_W = sigma w^-1/2 (1 - j), finite Warburg R tanh(s)/s or R coth(s)/s, s = (j w tau)^1/2",
               "k0 = R T / (n^2 F^2 A R_ct C)", "C = (Q R^(1-n))^(1/n), R in parallel with the CPE", "C = Q^(1/n) (Rs Rp/(Rs + Rp))^((1-n)/n)",
               "linear Kramers-Kronig fit with M RC elements, M set by mu = 0.85", "AICc = N ln(chi2/N) + 2p + 2p(p+1)/(N-p-1), N = 2 x frequencies", "3.3 x SD of the blank (or residual SD) / slope"),
  Source = c("Lukaszewski et al., Int. J. Electrochem. Sci. 2016, 11, 4442, and Trasatti and Petrii, Pure Appl. Chem. 1991, 63, 711",
             "Lukaszewski et al. 2016 (C_s mostly 20 to 40 µF/cm², error up to 100%)", "Microchim. Acta 2023, 190, 251",
             "geometry", "Bard and Faulkner, and Microchim. Acta 2023, 190, 251", "Bard and Faulkner, and Microchim. Acta 2023, 190, 251",
             "Matsuda and Ayabe, and Microchim. Acta 2023, 190, 251, eq. 4 and Fig. 6. K, the peak width and ΔEp from a Butler–Volmer simulation in this app, checked against Randles–Ševčík and Nicholson 1965",
             "Bard and Faulkner (derivation with one form in the bulk). 1 + 1 mM ferri/ferro read as one form gives about 1.6x the area in this app's simulation",
             "Microchim. Acta 2023, 190, 251, eq. 9", "Nicholson, Anal. Chem. 1965, 37, 1351, and Lavagnini et al., Electroanalysis 2004, 16, 505",
             "Klingler and Kochi, J. Phys. Chem. 1981, 85, 1731", "Bard and Faulkner", "Laviron, J. Electroanal. Chem. 1979, 101, 19, eqs. 20 to 22",
             "Bard and Faulkner", "Bard and Faulkner (disk: Saito 1968)", "Shoup and Szabo, J. Electroanal. Chem. 1982, 140, 237",
             "Anson, and Microchim. Acta 2023, 190, 251, eq. 1", "Bard and Faulkner", "Bard and Faulkner", "Bard and Faulkner", "Bard and Faulkner",
             "widely used RRDE relations (oxygen reduction)", "Bard and Faulkner", "Bard and Faulkner",
             "Lasia, and Orazem and Tribollet", "Microchim. Acta 2023, 190, 251, eq. 17", "Hsu and Mansfeld (as in the Zahner documentation, eq. 7)", "Brug et al., J. Electroanal. Chem. 1984, 176, 275 (Zahner eq. 8)",
             "Schonleber et al., Electrochim. Acta 2014, 131, 20", "Akaike information criterion with the small-sample correction", "ICH Q2"))

# The same equations, typeset for the Equations tab (LaTeX for MathJax), one row per row of equation_sources() and
# in the same order. formula is the equation alone (empty where the entry is a rule in words); note holds the
# conditions and words that go with it. Together they say exactly what the Equation column says.
equation_display = data.frame(formula = c(
  r"---(\mathrm{ECSA} = \dfrac{Q}{q_\mathrm{ref}})---",
  r"---(C_\mathrm{dl} = \text{slope of } \dfrac{i_a - i_c}{2} \text{ against } v \qquad \mathrm{ECSA} = \dfrac{C_\mathrm{dl}}{C_s})---",
  "",
  r"---(\begin{gathered} \text{disk } \dfrac{\pi d^2}{4} \qquad \text{plate } LW \text{ per face} \qquad \text{strip edges } 2Lt + Wt \\ \text{wire } \pi d L \ \left(+\,\dfrac{\pi d^2}{4}\right) \qquad \text{hemisphere } 2\pi r^2 \end{gathered})---",
  r"---(i_p = 0.4463\, nFAC \left(\dfrac{nFvD}{RT}\right)^{1/2} = 2.69\times10^{5}\, n^{3/2} A D^{1/2} C v^{1/2})---",
  r"---(i_p = 0.4958\, nFAC \left(\dfrac{\alpha n F v D}{RT}\right)^{1/2} = 2.99\times10^{5}\, n\, \alpha^{1/2} A D^{1/2} C v^{1/2})---",
  r"---(i_p = i_p^\mathrm{rev}\, K(\Lambda, \alpha) \qquad \Lambda = \dfrac{k^0}{(DFv/RT)^{1/2}})---",
  "",
  r"---(|E_p - E_{p/2}| = \dfrac{1.857\, RT}{\alpha n F} = \dfrac{47.7}{\alpha n}\ \mathrm{mV})---",
  r"---(\psi = \dfrac{-0.6288 + 0.0021\, n\Delta E_p}{1 - 0.017\, n\Delta E_p} \qquad k^0 = \psi \left(\dfrac{\pi D n v F}{RT}\right)^{1/2})---",
  r"---(k^0 = 2.18 \left(\dfrac{\alpha D n v F}{RT}\right)^{1/2} \exp\!\left(-\dfrac{\alpha^2 n F \Delta E_p}{RT}\right))---",
  r"---(\Gamma = \dfrac{Q}{nFA} \qquad i_p = \dfrac{n^2 F^2 v A \Gamma}{4RT})---",
  r"---(\begin{gathered} \dfrac{\mathrm{d}E_{pc}}{\mathrm{d}\log v} = -\dfrac{2.303\, RT}{\alpha n F} \\ \log k_s =\alpha \log(1 - \alpha) + (1 - \alpha)\log\alpha - \log\dfrac{RT}{nFv} - \dfrac{\alpha(1 - \alpha)\, nF\Delta E_p}{2.303\, RT} \end{gathered})---",
  r"---(i = \dfrac{nFAC D^{1/2}}{(\pi t)^{1/2}})---",
  r"---(\text{disk } i = 4nFDCr \qquad \text{hemisphere } i = 2\pi nFDCr \qquad \text{sphere } i = 4\pi nFDCr)---",
  r"---(\begin{gathered} i = 4nFDCr\, f(\tau) \qquad \tau = \dfrac{4Dt}{r^2} \\ f = 0.7854 + 0.8862\,\tau^{-1/2} + 0.2146\, \exp\!\left(-0.7823\,\tau^{-1/2}\right) \end{gathered})---",
  r"---(Q = \dfrac{2nFAC D^{1/2} t^{1/2}}{\pi^{1/2}} + Q_\mathrm{dl} + nFA\Gamma)---",
  r"---(i\,\tau^{1/2} = \dfrac{nFAC (\pi D)^{1/2}}{2})---",
  r"---(i_L = 0.620\, nFA D^{2/3} \omega^{1/2} \nu^{-1/6} C)---",
  r"---(\dfrac{1}{i} = \dfrac{1}{i_k} + \dfrac{1}{B\,\omega^{1/2}} \qquad B = 0.620\, nFA D^{2/3} \nu^{-1/6} C)---",
  r"---(N = -\dfrac{I_\mathrm{ring}}{I_\mathrm{disk}})---",
  r"---(n = \dfrac{4|I_D|}{|I_D| + |I_R|/N} \qquad \%\,\mathrm{H_2O_2} = \dfrac{200\, |I_R|/N}{|I_D| + |I_R|/N})---",
  r"---(b = \dfrac{2.303\, RT}{\alpha n F})---",
  "",
  r"---(\begin{gathered} Z_\mathrm{CPE} = \dfrac{1}{Q\,(j\omega)^n} \qquad Z_W = \sigma\,\omega^{-1/2}(1 - j) \\ \text{finite Warburg: } R\,\dfrac{\tanh s}{s} \text{ or } R\,\dfrac{\coth s}{s} \qquad s = (j\omega\tau)^{1/2} \end{gathered})---",
  r"---(k^0 = \dfrac{RT}{n^2 F^2 A R_\mathrm{ct} C})---",
  r"---(C = \left(Q R^{1-n}\right)^{1/n})---",
  r"---(C = Q^{1/n} \left(\dfrac{R_s R_p}{R_s + R_p}\right)^{(1-n)/n})---",
  "",
  r"---(\mathrm{AICc} = N \ln\!\left(\dfrac{\chi^2}{N}\right) + 2p + \dfrac{2p(p+1)}{N - p - 1})---",
  r"---(\mathrm{LOD} = \dfrac{3.3\, s}{\text{slope}})---"),
note = c(
  "", "",
  r"---(Slope of \(i_p\) against \(v^{1/2}\) with the Randles–Ševčík equation and a known \(D\).)---",
  "",
  "At 25 °C.", "At 25 °C.",
  r"---(\(\Lambda\) from \(n\Delta E_p\) and \(\alpha\) from \(E_p - E_{p/2}\). Area from the slope of \(i_p/K\) against \(v^{1/2}\).)---",
  "Randles–Ševčík assumes one form in the bulk. With both forms present, use Cottrell or Anson after a step from the rest potential.",
  "At 25 °C.",
  r"---(For \(n\Delta E_p\) from 61 to 200 mV.)---",
  r"---(For \(n\Delta E_p\) above 150 mV.)---",
  "",
  r"---(\(\alpha\) from the slope of \(E_{pc}\) against \(\log v\). For \(n\Delta E_p\) above 200 mV.)---",
  "", "", "", "", "",
  r"---(\(\omega\) in rad/s.)---",
  "",
  "On the limiting plateau.",
  "", "",
  "Tafel line extrapolated to the equilibrium potential.",
  "", "",
  r"---(\(R\) in parallel with the CPE.)---",
  "",
  r"---(Linear Kramers–Kronig fit with \(M\) RC elements, \(M\) set by \(\mu = 0.85\).)---",
  r"---(\(N = 2 \times\) the number of frequencies.)---",
  r"---(\(s\) is the SD of the blank (or the residual SD).)---"))
