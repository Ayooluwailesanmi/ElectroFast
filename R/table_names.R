# Readable column headers for the tables shown in the app and the reports. The calculations keep short code names
# (Epa_V, k0_nicholson_cm_s), which CSV downloads also keep so that scripts can read them back; round_df() turns them
# into headers such as "Epa (V)" and "k⁰, Nicholson (cm/s)" for display. A name with a space or a bracket is taken as
# readable already and left as it is.

# unit suffixes, longest first so that cm2_s is matched before s
table_units = c(
  "uA_mM_cm2" = "µA mM⁻¹ cm⁻²", "A_per_rad05" = "A rad⁻¹ᐟ²", "uC_cm2" = "µC/cm²", "uF_cm2" = "µF/cm²", "mol_cm2" = "mol/cm²",
  "mV_dec" = "mV/decade", "mV_min" = "mV/min", "cm2_s" = "cm²/s", "cm_s" = "cm/s", "V_s" = "V/s", "uA_mM" = "µA/mM",
  "cm2" = "cm²", "cm" = "cm", "um" = "µm", "mV" = "mV", "V" = "V", "mA" = "mA", "uA" = "µA", "nA" = "nA", "pA" = "pA", "A" = "A",
  "uC" = "µC", "fC" = "fC", "C" = "C", "uF" = "µF", "F" = "F", "ohm" = "Ω", "ms" = "ms", "s" = "s", "min" = "min", "Hz" = "Hz",
  "mM" = "mM", "M" = "M", "pct" = "%", "rpm" = "rpm")

# whole names, and symbols, that the general rule would not write well
table_labels = c(
  file = "File", day = "Day", date = "Date", note = "Note", conc = "Concentration", conc_M = "Concentration (M)",
  n = "n", N = "N", N_sd = "SD of N", K = "K", z = "z", r2 = "R²", se = "SE", sd = "SD", mean = "Mean", AICc = "AICc", dAICc = "ΔAICc", BIC = "BIC",
  chi2_red = "Reduced χ²", nrmse_pct = "NRMSE (%)", rel_se_pct = "Error (%)", worst_rel_se_pct = "Largest parameter error (%)",
  lgK = "log K", pct_nernst = "% of Nernstian", nernst_slope = "Nernstian slope (mV/decade)", decades = "Decades",
  Epa_V = "Epa (V)", Epc_V = "Epc (V)", ipa_A = "ipa (A)", ipc_A = "ipc (A)", dEp_mV = "ΔEp (mV)", E_half_V = "E½ (V)",
  ratio_ipa_ipc = "ipa / ipc", regime = "Regime", v_V_s = "Scan rate (V/s)",
  scan_rate_V_s = "Scan rate (V/s)", k0_nicholson_cm_s = "k⁰, Nicholson (cm/s)", k0_klingler_kochi_cm_s = "k⁰, Klingler-Kochi (cm/s)",
  k0_K_table_cm_s = "k⁰, K table (cm/s)", peak = "Peak", against = "Against", slope_se = "SE of the slope", slope_unit = "Slope unit",
  intercept_se_A = "SE of the intercept (A)", intercept_pct_of_smallest_ip = "Intercept (% of the smallest peak)", loglog_slope = "Log-log slope", loglog_slope_se = "SE of the log-log slope", alpha_K_table = "α from the K table", k0_note = "k⁰ note", alpha_from_width = "α from the peak width (irreversible)",
  E0_mV = "E⁰ (mV)", E_V = "E (V)", I_A = "I (A)", E_quarter_V = "E at τ/4 (V)", tau_s = "τ (s)", ECSA_cm2 = "ECSA (cm²)",
  Q_uC = "Q (µC)", q_ref_uC_cm2 = "Charge density (µC/cm²)", D_cm2_s = "D (cm²/s)", b_bkg_M = "Background b (M)",
  i_k_A = "Kinetic current (A)", I_disk_A = "Disk current (A)", I_ring_A = "Ring current (A)", ring_E_V = "Ring potential (V)",
  FE_mean_pct = "Faradaic efficiency, mean (%)", FE_sd_pct = "Faradaic efficiency, SD (%)", H2O2_mean_pct = "H₂O₂ yield, mean (%)",
  H2O2_sd_pct = "H₂O₂ yield, SD (%)", n_mean = "n, mean", n_sd = "n, SD", n_ref = "Reference sessions", n_from_slope = "n from the slope",
  slope_A_per_rad05 = "Slope (A rad⁻¹ᐟ²)", lod = "LOD", loq = "LOQ", at_bound = "At bound", recommended = "Recommended",
  a_primary_cross = "Primary-ion activity at the crossing", a_interferent = "Interferent activity", E_interferent = "E, interferent (mV)",
  E_primary_mean = "E, primary ion, mean (mV)", cross_reactivity_pct = "Cross-reactivity (%)", within_tolerance = "Within tolerance",
  signal_analyte = "Signal, target alone", signal_mixture = "Signal, with interferent", interferent_concentration = "Interferent concentration",
  median_charge_fC = "Median charge (fC)", median_amplitude_pA = "Median amplitude (pA)", median_duration_ms = "Median duration (ms)",
  median_height_pA = "Median height (pA)", unsettled_share = "Share unsettled", low_end_ratio = "Low-end ratio", background_ratio = "Background ratio",
  slope_min = "Lowest slope (mV/decade)", slope_max = "Highest slope (mV/decade)", times = "Reading times", rates = "Rotation rates", rpm = "Rotation rate (rpm)",
  t_start_s = "Start (s)", t_end_s = "End (s)", t_peak_s = "Peak time (s)", t_step_s = "Step time (s)", largest_pA = "Largest step (pA)",
  flat_mV = "Plateau (mV)", flat_sd = "Plateau SD (mV)", run_s = "Run length (s)", settled_s = "Settled after (s)", bin_ms = "Bin width (ms)",
  interferent_day = "Interferent day", to = "To", to_M = "To (M)", from = "From", from_M = "From (M)",
  pooled_slope = "Pooled slope (mV/decade)", pooled_slope_se = "SE of the pooled slope (mV/decade)",
  between_day_baseline_sd = "Between-day SD of the baseline (mV)", between_day_slope_sd = "Between-day SD of the slope (mV/decade)")

nice_names = function(x) vapply(x, function(nm) {
  if (nm %in% names(table_labels)) return(table_labels[[nm]])
  if (grepl("[ (]", nm) || !nzchar(nm)) return(nm)
  unit = ""
  for (u in names(table_units)) if (endsWith(nm, paste0("_", u)) && nchar(nm) > nchar(u) + 1) { unit = table_units[[u]]; nm = substr(nm, 1, nchar(nm) - nchar(u) - 1); break }
  w = gsub("_", " ", nm)
  w = sub("\\bk0\\b", "k⁰", sub("\\bd([A-Z])", "Δ\\1", w))
  if (grepl("^[a-z]{2,}", w)) w = paste0(toupper(substr(w, 1, 1)), substring(w, 2))
  if (nzchar(unit)) sprintf("%s (%s)", w, unit) else w
}, "", USE.NAMES = FALSE)
