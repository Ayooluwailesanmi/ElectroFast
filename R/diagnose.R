# Possible causes for an ISE calibration session: the signs the app already measures, checked day by day, then combined
# into patterns with what to check. These are prompts for the user, not diagnoses: the same sign can have several causes.
# The cut-offs are this app's defaults, not published standards; the user can change them (th).
#   slope      % of the Nernstian slope over the fit range (default: flagged below 90 or above 110)
#   low end    slope over the three lowest fitted standards / slope over the whole fit (below 0.6)
#   background b of the day / the lowest b of all days, background-corrected fit (above 3)
#   check      end reading minus start reading of the check standard, X then Z (beyond 5 mV)
#   scatter    R2 of the fit (below 0.98)
#   settling   share of fitted standards still drifting at the reading time, Reading time tab rule (above 20%)
diagnose_defaults = list(slope_lo = 90, slope_hi = 110, low_end = 0.6, background = 3, drift_mV = 5, r2 = 0.98, unsettled = 0.2)

diagnose_ise = function(pts, res, bg = NULL, settle = NULL, x = "neg_log_a", read_s = 120, th = diagnose_defaults) {
  days = res$per_day$day
  low_ratio = vapply(days, function(dd) {
    p = pts[pts$day == dd & pts$fit & is.finite(pts$E_mV), ]; p = p[order(p$conc_M), ]
    if (nrow(p) < 5) return(NA_real_)
    # same convention as calibrate(): the plain coefficient of E against the chosen axis, so the ratio is sign-free
    lo = utils::head(p, 3); s_lo = unname(stats::coef(stats::lm(lo$E_mV ~ lo[[x]]))[2])
    s_lo / res$per_day$slope_mV_dec[res$per_day$day == dd]
  }, 0)
  b_ratio = if (!is.null(bg) && !is.null(bg$per_day) && nrow(bg$per_day)) {
    b = bg$per_day$b_bkg_M[match(days, bg$per_day$day)]; b / min(b, na.rm = TRUE) } else rep(NA_real_, length(days))
  check = vapply(days, function(dd) {
    p = pts[pts$day == dd, ]; xv = p$E_mV[toupper(p$letter) == "X"]; zv = p$E_mV[toupper(p$letter) == "Z"]
    if (length(xv) && length(zv)) zv[1] - xv[1] else NA_real_ }, 0)
  unsettled = vapply(days, function(dd) {
    if (is.null(settle)) return(NA_real_)
    s = settle[settle$day == dd & settle$fit, ]; if (!nrow(s)) return(NA_real_)
    mean(is.na(s$settled_s) | s$settled_s > read_s) }, 0)
  signs = data.frame(day = days, pct_nernst = res$per_day$pct_nernst, low_end_ratio = low_ratio, background_ratio = b_ratio,
                     check_drift_mV = check, r2 = res$per_day$r2, unsettled_share = unsettled)
  f = with(signs, data.frame(day,
    slope_low = !is.na(pct_nernst) & pct_nernst < th$slope_lo, slope_high = !is.na(pct_nernst) & pct_nernst > th$slope_hi,
    low_end = !is.na(low_end_ratio) & low_end_ratio < th$low_end, background = !is.na(background_ratio) & background_ratio > th$background,
    drift = !is.na(check_drift_mV) & abs(check_drift_mV) > th$drift_mV, scatter = !is.na(r2) & r2 < th$r2,
    unsettled = !is.na(unsettled_share) & unsettled_share > th$unsettled))
  causes = do.call(rbind, lapply(seq_len(nrow(f)), function(i) {
    g = f[i, ]; out = character(0)
    # the wording names only the signs that were flagged on that day
    if (g$low_end && g$background)
      out = c(out, "Low end flattening with the background up: more Li⁺ than intended in the dilute standards, or Li⁺ released at the membrane surface. Check: fresh dilute standards in clean tubes, rinsing between standards, the membrane surface (clean or soak).")
    else if (g$low_end)
      out = c(out, "Low end flattening: the dilute standards give less response than the rest of the line. Check: fresh dilute standards in clean tubes, rinsing between standards, the membrane surface, and the detection limit expected for this membrane.")
    if (g$slope_low && !g$low_end)
      out = c(out, "Slope low over the whole range. Check the reference electrode (filling solution, junction, and against a second reference). Also check the temperature and the membrane (cracks, seal).")
    if (g$slope_high)
      out = c(out, "Slope above Nernstian: a mixed potential or a changing ionic strength. A constant background electrolyte, or the fit against concentration, tests this.")
    if (g$drift)
      out = c(out, sprintf("Check standard moved by %.1f mV within the session: carry-over from the concentrated standards, or a drifting reference. Check: rinse and re-read the check standard after the highest standards, and compare the reference with a second one.", signs$check_drift_mV[i]))
    if (g$scatter) out = c(out, "Scatter about the line: check the order and labels of the standards, and standards read before the potential settled.")
    if (g$unsettled) out = c(out, "Several standards were still drifting at the reading time: read later or record longer (Reading time tab).")
    if (!length(out)) out = "No sign flagged."
    data.frame(day = g$day, signs = paste(names(g)[-1][unlist(g[-1])], collapse = ", "), check = out)
  }))
  causes$signs[causes$signs == ""] = "none"
  list(signs = signs, flags = f, causes = causes)
}
