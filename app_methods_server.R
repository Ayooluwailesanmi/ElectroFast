# Server side of the technique tabs (UI in R/app_methods_ui.R). Kept outside R/ because Shiny
# runs every file in R/ at start-up, and this one only works inside the app's server function,
# where it is sourced with local = TRUE, so it uses that function's input, output, session
# and its renderPlot (which falls back to a default size on hidden tabs).

g_theme = function() list(ggplot2::theme_classic(base_size = 12), watermark())
pal = function(n) rep(cycle_palette("#CC9900"), length.out = n)

# Files are read once per upload and kept by path
.methods_cache = new.env()
read_up = function(files) if (is.null(files) || !NROW(files)) NULL else lapply(seq_len(NROW(files)), function(i) {
  # the instrument choice is part of the key, so changing it reads the file again with the chosen reader
  inst = input$instrument %||% "auto"; key = paste(files$datapath[i], inst)
  if (is.null(.methods_cache[[key]])) assign(key, read_echem(files$datapath[i], files$name[i], inst), envir = .methods_cache)
  .methods_cache[[key]]
})
# One line per uploaded file saying how it was read: the instrument recognised (or chosen), the technique, the
# number of points, segments and scan rate; EMF16 recordings with their channels and length; impedance files with
# their frequencies. Uses the same cached read as the analysis, so a large file is not read twice.
format_labels = c(chi = "CH Instruments", pssession = "PalmSens session", mpr = "BioLogic EC-Lab binary (.mpr)", biologic = "BioLogic EC-Lab text (.mpt)",
                  gamry = "Gamry", pine = "Pine AfterMath", abf = "Axon pCLAMP (.abf)", text = "text/CSV (columns read by their names)")
readas_line = function(path, name, inst = "auto") {
  fmt = tryCatch(instrument_format(path, name, inst), error = function(e) e)
  if (inherits(fmt, "error")) return(paste("Not read:", conditionMessage(fmt)))
  if (fmt == "text" && any(grepl("^\\s*Time\\s*/\\s*s", readLines(path, n = 40, warn = FALSE)))) {
    z = tryCatch(read_emf16_grf(path), error = function(e) NULL)
    if (!is.null(z)) return(sprintf("Read as: Lawson Labs EMF16, %d channel(s), %d points, %.0f s", length(z$channels), nrow(z$wide), max(z$wide$time_s)))
  }
  key = paste(path, inst)
  x = .methods_cache[[key]] %||% tryCatch({ r = read_echem(path, name, inst); assign(key, r, envir = .methods_cache); r }, error = function(e) e)
  if (inherits(x, "error") || isTRUE(x$kind == "eis")) {
    e = tryCatch(read_eis(path, name, inst), error = function(e2) NULL)
    if (!is.null(e)) return(sprintf("Read as: %s, impedance, %d frequencies", format_labels[[fmt]], nrow(e$data)))
    return(paste("Not read:", if (inherits(x, "error")) conditionMessage(x) else "no data found"))
  }
  d = x$data; p = x$params
  parts = c(format_labels[[fmt]], if (identical(p$technique, "text export")) NA else p$technique %||% NA, sprintf("%d points", nrow(d)),
            if (!is.null(d$segment) && max(d$segment, na.rm = TRUE) > 1) sprintf("%d segments", max(d$segment, na.rm = TRUE)),
            if (isTRUE(is.finite(p$scan_rate_V_s))) sprintf("%g V/s", p$scan_rate_V_s))
  paste("Read as:", paste(parts[!is.na(parts) & nzchar(parts)], collapse = ", "))
}
for (id in c("dl_files", "sr_files", "tf_file", "pv_files", "oc_files", "ca_files", "cc_files", "cc_blank", "cp_files", "rd_files", "rr_cal_files",
             "rr_fe_files", "ei_file", "pq_file", "cal_files", "fim_files", "ssm_files", "ecsa_file", "sa_files")) local({ id = id
  output[[paste0(id, "_readas")]] = renderUI({ f = req(input[[id]]); k = seq_len(min(nrow(f), 20))
    lines = vapply(k, function(i) paste0(f$name[i], ": ", readas_line(f$datapath[i], f$name[i], input$instrument %||% "auto")), "")
    tagList(lapply(lines, tags$div), if (nrow(f) > 20) tags$div(sprintf("and %d more files", nrow(f) - 20))) })
})
safe_read = function(files) tryCatch(read_up(files), error = function(e) {
  notify_error(e, duration = 12); NULL })

# Editable concentration table shared by the calibration tabs
conc_table = function(files_id, unit_id, table_id) {
  tbl = reactiveVal(NULL)
  observeEvent(list(input[[files_id]], input[[unit_id]]), {
    f = req(input[[files_id]])
    v = conc_from_name(f$name, input[[unit_id]])
    tbl(data.frame(file = f$name, conc = v, include = TRUE,
                   note = ifelse(is.na(v), "no concentration found in the name: type it", "suggested from the file name")))
  })
  output[[table_id]] = DT::renderDT(DT::datatable(req(tbl()), rownames = FALSE, selection = "none",
    editable = list(target = "cell", disable = list(columns = c(0, 3))), options = list(pageLength = 50, dom = "t")))
  observeEvent(input[[paste0(table_id, "_cell_edit")]], {
    e = input[[paste0(table_id, "_cell_edit")]]; d = tbl()
    if (e$col == 1) { d$conc[e$row] = suppressWarnings(as.numeric(e$value)); d$note[e$row] = "entered by hand" }
    if (e$col == 2) d$include[e$row] = toupper(trimws(e$value)) %in% c("TRUE", "T", "1", "YES")
    tbl(d)
  })
  tbl
}
# points table with the fitted line appended, for the CSV downloads
add_line = function(p, cl) if (is.null(cl)) p else cbind(p, stats::setNames(cl, paste0("line_", names(cl))))
fmt_line = function(cl, unit, yunit = "A") {
  if (is.null(cl)) return(data.frame(Note = "At least 3 included standards with a concentration above 0 are needed for a line."))
  data.frame(Standards = cl$n, Range = sprintf("%g to %g %s", cl$from, cl$to, unit),
             Slope = sprintf("%.4g %s/%s", cl$slope, yunit, unit), Intercept = sprintf("%.4g %s", cl$intercept, yunit),
             `R²` = sprintf("%.4f", cl$r2), LOD = sprintf("%.3g %s", cl$lod, unit), `LOD basis` = cl$lod_basis, check.names = FALSE)
}
cal_plot_obj = function(pts, cl, unit, ylab) {
  g = ggplot2::ggplot(pts, ggplot2::aes(conc, y * 1e6)) + ggplot2::geom_point(ggplot2::aes(shape = used), size = 2.6)
  if (!is.null(cl)) g = g + ggplot2::geom_abline(intercept = cl$intercept * 1e6, slope = cl$slope * 1e6, colour = "#1A6FDF")
  g + ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "In the fit", `FALSE` = "Not in the fit"), name = NULL) +
    ggplot2::labs(x = sprintf("Concentration (%s)", unit), y = ylab,
                  subtitle = if (is.null(cl)) NULL else sprintf("Slope %.4g µA/%s, R² %.4f, LOD %.3g %s", cl$slope * 1e6, unit, cl$r2, cl$lod, unit)) +
    g_theme()
}
traces_plot = function(xs, x, y, xlab, ylab, yscale = 1e6) {
  d = do.call(rbind, lapply(xs, function(o) data.frame(file = o$file, xv = o$data[[x]], yv = o$data[[y]] * yscale)))
  g = ggplot2::ggplot(d, ggplot2::aes(xv, yv, colour = file)) + ggplot2::geom_path(linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = pal(length(xs)), name = NULL) + ggplot2::labs(x = xlab, y = ylab) +
    g_theme() + ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
  if (labels_on() && x == "E_V") g = g + cv_labels(d, sub("^.*\\((.*)\\).*$", "\\1", ylab))
  g
}

# "Label plots" on a voltammogram: on the last two sweeps of the first file, an arrow beside each sweep in its
# direction (4% of the current range off the curve, above an anodic sweep and below a cathodic one, so it does not
# hide the line) and the peak current with its potential. A maximum at the very end of a sweep is a vertex, not a
# peak, and is not labelled.
cv_labels = function(d, unit) {
  z = d[d$file == d$file[1], ]; s = cv_sweeps(z$xv); n = length(s$from)
  if (!n) return(list())
  yr = diff(range(z$yv, na.rm = TRUE)); out = list()
  for (i in utils::tail(seq_len(n), 2)) {
    w = z[s$from[i]:s$to[i], ]; m = nrow(w)
    if (m < 20) next
    an = w$xv[m] > w$xv[1]; dy = if (an) 0.04 * yr else -0.04 * yr
    a = w[round(0.35 * m), ]; b = w[round(0.45 * m), ]
    out = c(out, list(ggplot2::annotate("segment", x = a$xv, xend = b$xv, y = a$yv + dy, yend = b$yv + dy, linewidth = 0.7,
                                        arrow = grid::arrow(length = grid::unit(0.1, "in"), type = "closed"))))
    # the most prominent feature of the sweep, of either sign: a reduction peak on a cathodic sweep, or an oxidation
    # peak on it (ethanol on Pd oxidises on the return sweep too)
    dev = w$yv - stats::median(w$yv); k = which.max(abs(dev))
    if (k > 0.03 * m && k < 0.97 * m)
      out = c(out, list(ggplot2::annotate("point", x = w$xv[k], y = w$yv[k], size = 2),
                        # forward-sweep label to the right of its peak, return-sweep label to the left, so two
                        # close peaks do not print over each other
                        ggplot2::annotate("text", x = w$xv[k] + (if (an) 0.015 else -0.015) * diff(range(z$xv)),
                                          y = w$yv[k] + sign(dev[k]) * 0.05 * yr, size = 3.3, hjust = if (an) 0 else 1,
                                          label = sub("-", "−", sprintf("%.3g %s at %.3f V", w$yv[k], unit, w$xv[k])))))
  }
  out
}

# "Label plots" on Multi-Potential Steps traces: from the first file's settings, the potential of each step along
# the top, the cycle numbers, dashed lines between cycles, and the reading windows of the cycles that are used (the
# last avg_s of the last step of each cycle, after the cycles left out), as multistep_currents() reads them.
step_labels = function(x, ytop, yr, avg_s, skip) {
  st = x$params$settings
  k = sort(as.integer(sub("^E([0-9]+) \\(V\\)$", "\\1", grep("^E[0-9]+ \\(V\\)$", names(st), value = TRUE))))
  if (!length(k)) return(list())
  E = vapply(k, function(j) st[[sprintf("E%d (V)", j)]], 0); Tt = vapply(k, function(j) st[[sprintf("T%d (s)", j)]], 0)
  n = chi_setting(x, "Cycle", 1); per = sum(Tt)
  mid = rep((seq_len(n) - 1) * per, each = length(k)) + rep(cumsum(Tt) - Tt / 2, n)
  ends = seq_len(n) * per; used = seq_len(n) > min(skip, n - 1)
  list(ggplot2::annotate("rect", xmin = ends[used] - avg_s, xmax = ends[used], ymin = -Inf, ymax = Inf, fill = "#D55E00", alpha = 0.3),
       ggplot2::annotate("text", x = mid, y = ytop, label = sub("-", "−", sprintf("%+.2f V", rep(E, n))), size = 2.7),
       if (n > 1) ggplot2::annotate("text", x = (seq_len(n) - 0.5) * per, y = ytop + 0.07 * yr, label = paste("Cycle", seq_len(n)), size = 3.3, fontface = "bold"),
       if (n > 1) ggplot2::geom_vline(xintercept = seq_len(n - 1) * per, linetype = "dashed", colour = "grey60", linewidth = 0.4))
}

# Electrode panels (geometric area, radius) of the tabs that need an area
el = lapply(c(sa = "sa_el", dl = "dl_el", sr = "sr_el", tf = "tf_el", ca = "ca_el", cc = "cc_el", cp = "cp_el",
              rd = "rd_el", ei = "ei_el"), electrode_server)
area_of = function(k) el[[k]]()$area_cm2
T_K = function(p) (input[[paste0(p, "_T")]] %||% 25) + 273.15
# picking a literature probe fills in its D; "Other" leaves the box for the user
for (p in c("sr", "ca", "cc", "rd", "ei")) local({ p = p
  observeEvent(input[[paste0(p, "_probe")]], {
    D = redox_probes$D_cm2_s[redox_probes$key == input[[paste0(p, "_probe")]]]
    if (length(D) && !is.na(D)) updateNumericInput(session, paste0(p, "_D"), value = D)
  })
})
nz = function(v) if (is.null(v) || length(v) == 0) NA else v

# Surface area from a CV peak
observeEvent(input$sa_preset, {
  r = ecsa_references[ecsa_references$key == input$sa_preset, ]
  if (!nrow(r)) return()
  updateNumericInput(session, "sa_q", value = r$q_uC_cm2)
  updateRadioButtons(session, "sa_dir", selected = r$sweep)
  # default windows for Pd and Au CVs vs Ag/AgCl (1 M KCl); for the other
  # methods only the charge and sweep are set, because the regions depend on the reference
  # electrode and electrolyte used
  w = list(pd_ox424 = c(0.60, 0.10), au_ox = c(1.30, 0.55))[[r$key]]
  if (!is.null(w)) { updateNumericInput(session, "sa_from", value = w[1]); updateNumericInput(session, "sa_to", value = w[2]) }
  else showNotification("Set the potential window to the peak of your own CV.", duration = 8)
})
sa_fits = reactive({
  xs = req(safe_read(req(input$sa_files)))
  lapply(xs, function(x) ecsa_window(x, c(input$sa_from, input$sa_to), input$sa_q, area_of("sa"), input$sa_dir))
})
sa_res = reactive(do.call(rbind, lapply(sa_fits(), `[[`, "result")))
output$sa_table = renderTable(round_df(sa_res()))
sa_plot_obj = reactive({
  f = sa_fits()
  sw = do.call(rbind, lapply(f, function(z) cbind(file = z$result$file, z$sweep)))
  wn = do.call(rbind, lapply(f, function(z) cbind(file = z$result$file, z$window)))
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = wn, ggplot2::aes(E_V, ymin = I_A * 1e6, ymax = base_A * 1e6, fill = file), alpha = 0.25) +
    ggplot2::geom_path(data = sw, ggplot2::aes(E_V, I_A * 1e6, colour = file), linewidth = 0.6) +
    ggplot2::geom_line(data = wn, ggplot2::aes(E_V, base_A * 1e6, colour = file), linetype = "dashed") +
    ggplot2::scale_colour_manual(values = pal(length(f)), name = NULL) + ggplot2::scale_fill_manual(values = pal(length(f)), guide = "none") +
    ggplot2::labs(x = "Potential (V)", y = "Current (µA)", subtitle = sprintf("Last %s sweep, shaded: charge integrated (%g µC/cm²)", input$sa_dir, input$sa_q)) +
    g_theme()
})
output$sa_plot = renderPlot(sa_plot_obj())
output$sa_csv = downloadHandler("surface_area.csv", function(f) write.csv(sa_res(), f, row.names = FALSE))
output$sa_png = downloadHandler("surface_area.png", function(f) fig_save(f, sa_plot_obj(), width = 7, height = 5, dpi = 300))

# Scan rate
output$sr_points = renderTable({ p = sr_pk(); round_df(p[, c("file", "v_V_s", "Epa_V", "ipa_A", "Epc_V", "ipc_A")]) })
# both peaks against v and v^1/2, then the log-log slopes (same peak values as the Peaks and k0 table)
sr_rate = reactive({ p = sr_pk(); scan_rate_fits(p$v_V_s, p$ipa_A, p$ipc_A) })
output$sr_fits = renderTable({ f = sr_rate(); if (is.null(f)) data.frame(Note = "At least 3 scan rates are needed.") else
  round_df(f$fits) })
output$sr_loglog = renderTable({ f = sr_rate(); if (is.null(f)) NULL else round_df(f$loglog) })
output$sr_reading = renderUI(helpText(lapply(scan_rate_reading(sr_rate()), tags$p)))
sr_lin_plot_obj = reactive({ p = sr_pk()
  d = rbind(data.frame(peak = "Anodic", x = p$v_V_s, against = "italic(v)~(V/s)", ip = abs(p$ipa_A)), data.frame(peak = "Cathodic", x = p$v_V_s, against = "italic(v)~(V/s)", ip = abs(p$ipc_A)),
            data.frame(peak = "Anodic", x = sqrt(p$v_V_s), against = "italic(v)^{1/2}~((V/s)^{1/2})", ip = abs(p$ipa_A)), data.frame(peak = "Cathodic", x = sqrt(p$v_V_s), against = "italic(v)^{1/2}~((V/s)^{1/2})", ip = abs(p$ipc_A)))
  d = d[is.finite(d$ip), ]
  d$against = factor(d$against, levels = c("italic(v)~(V/s)", "italic(v)^{1/2}~((V/s)^{1/2})"))   # v first, as in the tables
  ggplot2::ggplot(d, ggplot2::aes(x, ip * 1e6, colour = peak)) + ggplot2::geom_point(size = 2.4) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6, fullrange = TRUE) +
    ggplot2::facet_grid(peak ~ against, scales = "free", labeller = ggplot2::label_parsed) + ggplot2::expand_limits(x = 0, y = 0) +
    ggplot2::scale_colour_manual(values = c(Anodic = "#F14040", Cathodic = "#1A6FDF"), guide = "none") +
    ggplot2::labs(x = NULL, y = expression("|"*italic(i)[p]*"| (µA)")) + g_theme() + ggplot2::theme(strip.background = ggplot2::element_blank()) })
output$sr_lin_plot = renderPlot(sr_lin_plot_obj())
sr_plot_obj = reactive({ p = sr_pk()
  d = rbind(data.frame(peak = "Anodic", v = p$v_V_s, ip = abs(p$ipa_A)), data.frame(peak = "Cathodic", v = p$v_V_s, ip = abs(p$ipc_A)))
  d = d[is.finite(d$ip) & d$ip > 0, ]; f = sr_rate()
  ggplot2::ggplot(d, ggplot2::aes(v, ip * 1e6, colour = peak)) + ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() + ggplot2::scale_colour_manual(values = c(Anodic = "#F14040", Cathodic = "#1A6FDF"), name = NULL) +
    ggplot2::labs(x = "Scan rate (V/s, log scale)", y = expression("|"*italic(i)[p]*"| (µA, log scale)"),
                  subtitle = if (is.null(f)) NULL else paste(sprintf("%s %.3f", tools::toTitleCase(f$loglog$peak), f$loglog$loglog_slope), collapse = ", ")) + g_theme()
})
output$sr_plot = renderPlot(sr_plot_obj())
sr_win = reactive(if (isTRUE(input$sr_use_window)) c(input$sr_from, input$sr_to))
# Peaks of every CV with the regime and k0 from the peak separation
sr_pk = reactive({
  xs = req(safe_read(req(input$sr_files)))
  p = do.call(rbind, lapply(xs, cv_peaks, window = sr_win(), T_K = T_K("sr")))
  p = p[order(p$v_V_s), ]
  n = input$sr_n %||% 1; w = if (identical(input$sr_dir, "cathodic")) p$Epc_minus_Ep2_mV else p$Epa_minus_Ep2_mV
  p$regime = cv_regime(p$dEp_mV, n)
  # |Ep - Ep/2| = 1.857 RT/(alpha n F) holds for an irreversible wave only; a reversible wave is 56.5 mV wide whatever alpha is
  p$alpha_from_width = ifelse(p$regime == "irreversible", alpha_from_peak_width(w, n, T_K("sr")), NA)
  p$k0_nicholson_cm_s = k0_nicholson(p$dEp_mV, p$v_V_s, input$sr_D, n, T_K("sr"))
  p$k0_klingler_kochi_cm_s = k0_klingler_kochi(p$dEp_mV, p$v_V_s, input$sr_D, input$sr_alpha, n, T_K("sr"))
  # alpha and k0 from the simulated K(Lambda, alpha) table (peak separation and width together), quasi-reversible CVs
  kq = rs_quasi_K(p$dEp_mV, w, n, T_K("sr")); q = p$regime == "quasi-reversible"
  p$alpha_K_table = ifelse(q, kq$alpha, NA)
  p$k0_K_table_cm_s = ifelse(q, kq$psi * sqrt(pi * input$sr_D * n * p$v_V_s * F_const / (R_const * T_K("sr"))), NA)
  # near the reversible limit 1 to 2 mV of peak separation (sampling, switching potential) changes k0 by tens of percent
  p$k0_note = ifelse(n * p$dEp_mV < 70 & (q | is.finite(p$k0_nicholson_cm_s)), "nΔEp below 70 mV: k⁰ is very sensitive to ΔEp, so treat it as an estimate", "")
  p
})
output$sr_peaks = renderTable(round_df(sr_pk()[, c("file", "v_V_s", "Epa_V", "ipa_A", "Epc_V", "ipc_A", "dEp_mV", "E_half_V", "ratio_ipa_ipc",
                                                   "regime", "alpha_from_width", "alpha_K_table", "k0_nicholson_cm_s", "k0_klingler_kochi_cm_s", "k0_K_table_cm_s", "k0_note")]))
sr_cv_plot_obj = reactive(traces_plot(req(safe_read(req(input$sr_files))), "E_V", "I_A", "Potential (V)", "Current (µA)"))
output$sr_cv_plot = renderPlot(sr_cv_plot_obj())
# Electroactive area or D from the forward peak against v^1/2 (diffusing species), or the
# surface coverage from the peak charge (adsorbed species)
sr_ip = reactive({ p = sr_pk(); data.frame(v = p$v_V_s, ip = abs(if (input$sr_dir == "anodic") p$ipa_A else p$ipc_A)) })
sr_area_tab = reactive({
  if (input$sr_species == "ads") {
    xs = req(safe_read(req(input$sr_files))); A = area_of("sr")
    return(do.call(rbind, lapply(xs, function(x) {
      q = peak_charge(last_sweep(x, input$sr_dir), c(input$sr_qfrom, input$sr_qto), x$params$scan_rate_V_s)$Q_C
      G = surface_coverage(q, input$sr_n, A)
      data.frame(file = x$file, v_V_s = x$params$scan_rate_V_s, Q_uC = q * 1e6, Gamma_mol_cm2 = G,
                 ip_ideal_reversible_A = ip_surface(G, x$params$scan_rate_V_s, input$sr_n, A, T_K("sr")))
    })))
  }
  # Randles-Sevcik assumes only one form in the bulk. With both forms the CV starts with the other form
  # already being converted at the electrode, and the peaks read from zero come out larger (about 1.6x
  # for 1 + 1 mM ferri/ferro cycled from -0.2 V, by simulation), so no area or D is given from them.
  validate(need(!isTRUE(input$sr_mix), paste(
    "Randles–Ševčík needs a solution with only one form of the couple (for example 1 mM ferricyanide alone).",
    "With both forms present the peaks are larger than the equation assumes and the area would be too high.",
    "Use the Chronoamperometry tab (Cottrell) or the Chronocoulometry tab (Anson) instead. After a step from the rest potential",
    "to a potential where one form is fully converted, the current depends only on that form's concentration.")))
  d = sr_ip(); quasi = input$sr_rev == "quasi"; rev = input$sr_rev != "irr"
  # quasi-reversible: each peak is divided by its K(Lambda, alpha), then the reversible equation is used (Trachioti 2023, Fig. 6)
  kq = NULL
  if (quasi) {
    p = sr_pk(); w = if (input$sr_dir == "anodic") p$Epa_minus_Ep2_mV else p$Epc_minus_Ep2_mV
    kq = rs_quasi_K(p$dEp_mV, w, input$sr_n, T_K("sr"))
    kq$K[p$regime == "reversible"] = 1
    validate(need(!any(p$regime == "irreversible", na.rm = TRUE), "Some CVs are irreversible (n dEp above 200 mV): use the irreversible form for those."),
             need(all(is.finite(kq$K)), "K could not be read for every CV (peak separation outside 59 to 200 mV): check the peaks table."))
    d$ip = d$ip / kq$K
  }
  s = if (nrow(d) >= 2) unname(stats::coef(stats::lm(ip ~ sqrt(v), d))[2]) else d$ip / sqrt(d$v)
  r2 = if (nrow(d) >= 3) summary(stats::lm(ip ~ sqrt(d$v), d))$r.squared else NA
  out = rs_solve(s, input$sr_n, input$sr_C, if (input$sr_solve == "area") input$sr_D else NA,
                 if (input$sr_solve == "D") area_of("sr") else NA, T_K("sr"), rev, input$sr_alpha)
  note = if (quasi && any(kq$alpha_at_table_edge, na.rm = TRUE)) "peak width outside the table (alpha 0.3 to 0.7): alpha held at the edge"
    else if (!quasi && rev && any(sr_pk()$regime == "quasi-reversible", na.rm = TRUE))
      "some CVs are quasi-reversible. Their peaks are lower than reversible ones, so the reversible form underestimates the area. Use the quasi-reversible form"
    else ""
  parts = list(data.frame(scan_rates = nrow(d), slope_A_per_sqrtV_s = s, r2 = r2, form = c(rev = "reversible", quasi = "quasi-reversible (K)", irr = "irreversible")[[input$sr_rev]],
                           geometric_area_cm2 = area_of("sr")), out,
               if (quasi) data.frame(K_per_scan_rate = paste(signif(kq$K, 3), collapse = ", "), alpha_from_width = paste(signif(kq$alpha, 2), collapse = ", ")),
               data.frame(note = note))
  do.call(cbind, Filter(Negate(is.null), parts))
})
output$sr_area = renderTable(round_df(sr_area_tab()))
sr_rs_plot_obj = reactive({ d = sr_ip()
  ggplot2::ggplot(d, ggplot2::aes(sqrt(v), ip * 1e6)) + ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#1A6FDF") +
    ggplot2::labs(x = expression(italic(v)^{1/2}*" ((V/s)"^{1/2}*")"), y = expression("|"*italic(i)[p]*"| (µA)")) + g_theme() })
output$sr_rs_plot = renderPlot(sr_rs_plot_obj())
sr_lav = reactive({ p = sr_pk(); laviron(p$v_V_s, p$Epa_V, p$Epc_V, input$sr_n, T_K("sr")) })
output$sr_lav = renderTable({ l = sr_lav(); if (is.null(l)) data.frame(Note = "At least 2 scan rates with dEp of 200/n mV or more are needed.") else round_df(l) })
sr_lav_plot_obj = reactive({ p = sr_pk()
  d = rbind(data.frame(v = p$v_V_s, E = p$Epa_V, peak = "anodic"), data.frame(v = p$v_V_s, E = p$Epc_V, peak = "cathodic"))
  ggplot2::ggplot(d, ggplot2::aes(v, E, colour = peak)) + ggplot2::geom_point(size = 2.6) + ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = c(anodic = "#F14040", cathodic = "#1A6FDF"), name = NULL) +
    ggplot2::labs(x = "Scan rate (V/s, log scale)", y = "Peak potential (V)") + g_theme() })
output$sr_lav_plot = renderPlot(sr_lav_plot_obj())
output$sr_csv = downloadHandler("cv_analysis.csv", function(f) {
  p = sr_pk(); a = tryCatch(sr_area_tab(), error = function(e) NULL)
  write.csv(cbind(p, if (!is.null(a) && nrow(a) == 1) a[rep(1, nrow(p)), ] else NULL), f, row.names = FALSE) })

# Tafel
tf_x = reactive(req(safe_read(req(input$tf_file)))[[1]])
tf_one = reactive(tafel_fit(tf_x(), c(input$tf_from, input$tf_to)))
tf_all = reactive(tafel_scan(tf_x()))
output$tf_table = renderTable({ t = tf_one()
  if (is.null(t)) return(data.frame(Note = "Fewer than 3 points with positive current in the window."))
  data.frame(Window = sprintf("%.3f to %.3f V", t$from_V, t$to_V), Points = t$n, `Tafel slope (mV/dec)` = sprintf("%.2f", t$tafel_mV_dec),
             `R²` = sprintf("%.4f", t$r2), check.names = FALSE) })
output$tf_range = renderTable({ a = tf_all(); if (is.null(a) || !nrow(a)) return(NULL)
  do.call(rbind, lapply(split(a, a$width_mV), function(d) data.frame(`Window width (mV)` = d$width_mV[1], Windows = nrow(d),
    `Lowest slope` = sprintf("%.1f", min(d$tafel_mV_dec)), `Highest slope` = sprintf("%.1f", max(d$tafel_mV_dec)), check.names = FALSE))) })
tf_plot_obj = reactive({
  d = tf_x()$data; d = d[d$I_A > 0, ]
  ggplot2::ggplot(d, ggplot2::aes(E_V, log10(I_A))) + ggplot2::geom_path(linewidth = 0.5) +
    ggplot2::annotate("rect", xmin = input$tf_from, xmax = input$tf_to, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#F14040") +
    ggplot2::labs(x = "Potential (V)", y = expression(log[10]*"("*italic(I)*" / A)"), subtitle = "Shaded: fitted window") + g_theme()
})
output$tf_plot = renderPlot(tf_plot_obj())
# alpha n from the slope, and the exchange current when an equilibrium potential is given
output$tf_kin = renderTable({
  d = tf_x()$data; k = d$E_V >= min(input$tf_from, input$tf_to) & d$E_V <= max(input$tf_from, input$tf_to) & d$I_A > 0
  if (sum(k) < 3) return(NULL)
  m = stats::lm(log10(d$I_A[k]) ~ d$E_V[k]); cf = unname(stats::coef(m))
  t = tafel_kinetics(cf[2], cf[1], nz(input$tf_Eeq), area_of("tf"), T_K("tf"))
  data.frame(`alpha n` = sprintf("%.3f", t$alpha_n), `Exchange current (A)` = if (is.na(t$i0_A)) "enter the equilibrium potential" else sprintf("%.3g", t$i0_A),
             `Exchange current density (A/cm²)` = if (is.na(t$j0_A_cm2)) "" else sprintf("%.3g", t$j0_A_cm2), check.names = FALSE)
})
output$tf_csv = downloadHandler("tafel_windows.csv", function(f) write.csv(tf_all(), f, row.names = FALSE))

# Pulse voltammetry calibration
pv_tbl = conc_table("pv_files", "pv_unit", "pv_conc")
pv_x = reactive(req(safe_read(req(input$pv_files))))
pv_pts = reactive({
  r = do.call(rbind, lapply(pv_x(), pulse_response, mode = input$pv_mode, E_read = input$pv_E, window = c(input$pv_from, input$pv_to)))
  t = req(pv_tbl()); r$conc = t$conc[match(r$file, t$file)]; r$include = t$include[match(r$file, t$file)]
  r$y = r$I_A; r$used = r$include & !is.na(r$conc) & r$conc > 0 & r$conc >= input$pv_cfrom & r$conc <= input$pv_cto
  r[order(r$conc), ]
})
pv_line = reactive({ p = pv_pts(); p = p[p$include, ]; calibration_line(p$conc, p$y, input$pv_cfrom, input$pv_cto) })
output$pv_fit = renderTable(fmt_line(pv_line(), input$pv_unit))
output$pv_points = renderTable(round_df(pv_pts()[, c("file", "conc", "E_V", "I_A", "used")]))
pv_plot_obj = reactive(cal_plot_obj(pv_pts(), pv_line(), input$pv_unit, "Current (µA)"))
output$pv_plot = renderPlot(pv_plot_obj())
pv_traces_obj = reactive(traces_plot(pv_x(), "E_V", "I_A", "Potential (V)", "Current (µA)"))
output$pv_traces = renderPlot(pv_traces_obj())
output$pv_csv = downloadHandler("pulse_calibration.csv", function(f) write.csv(add_line(pv_pts(), pv_line()), f, row.names = FALSE))
output$pv_png = downloadHandler("pulse_calibration.png", function(f) fig_save(f, pv_plot_obj(), width = 6.5, height = 4.5, dpi = 300))

# Open-circuit potential
oc_x = reactive(req(safe_read(req(input$oc_files))))
oc_res = reactive(do.call(rbind, lapply(oc_x(), ocp_summary, last_s = input$oc_last)))
output$oc_table = renderTable(round_df(oc_res()))
oc_plot_obj = reactive(traces_plot(oc_x(), "t_s", "E_V", "Time (s)", "Potential (V)", 1))
output$oc_plot = renderPlot(oc_plot_obj())
output$oc_csv = downloadHandler("ocp.csv", function(f) write.csv(oc_res(), f, row.names = FALSE))

# Chronoamperometry
ca_tbl = conc_table("ca_files", "ca_unit", "ca_conc")
ca_x = reactive(req(safe_read(req(input$ca_files))))
ca_eq = function() input$ca_mode %in% c("cottrell", "ume", "spikes", "esteps")
ca_res = reactive({
  xs = ca_x(); A = area_of("ca"); g = el$ca()
  if (input$ca_mode == "spikes") return(do.call(rbind, lapply(ca_spk(), function(z) spike_summary(z$file, z$t, z$res, input$ca_k))))
  if (input$ca_mode == "esteps") return(do.call(rbind, lapply(ca_est(), function(z) step_summary(z$file, z$t, z$res, input$ca_ek))))
  r = switch(input$ca_mode,
    multistep = do.call(rbind, lapply(xs, multistep_currents, avg_s = input$ca_avg, skip_cycles = input$ca_skip)),
    steps = { s = lapply(xs, step_currents, avg_s = input$ca_avg); k = Reduce(union, lapply(s, names))
              do.call(rbind, lapply(s, function(z) { z[setdiff(k, names(z))] = NA; z[k] })) },
    charge = do.call(rbind, lapply(xs, charge_passed, area_cm2 = A, q_ml_uC_cm2 = input$ca_qml)),
    cottrell = do.call(rbind, lapply(xs, cottrell_fit, t_from = input$ca_tfrom, t_to = input$ca_tto, n = input$ca_n, C_mM = input$ca_C,
                                     D_cm2_s = if (input$ca_solve == "area") input$ca_D else NA, area_cm2 = if (input$ca_solve == "D") A else NA)),
    ume = do.call(rbind, lapply(xs, function(x) { d = x$data; i = mean(d$I_A[d$t_s >= max(d$t_s) - input$ca_last])
      cbind(file = x$file, ume_steady(i, input$ca_n, if (input$ca_umesolve == "D") NA else input$ca_D,
                                      if (input$ca_umesolve == "C") NA else input$ca_C, if (input$ca_umesolve == "r") NA else g$r_cm, input$ca_geom)) })))
  if (ca_eq()) return(r)
  t = req(ca_tbl()); r$conc = t$conc[match(r$file, t$file)]; r$include = t$include[match(r$file, t$file)]
  r[order(r$conc), ]
})
ca_y = reactive({ r = ca_res()
  if (input$ca_mode == "multistep") r$I_mean_A else if (input$ca_mode == "steps") r[[paste0("I_step", input$ca_step, "_A")]] else NULL })
ca_line = reactive({
  y = ca_y(); if (is.null(y)) return(NULL)
  r = ca_res(); k = r$include
  blank = if (input$ca_mode == "multistep") r$I_sd_A[k & !is.na(r$conc) & r$conc == 0] else numeric()
  calibration_line(r$conc[k], y[k], input$ca_cfrom, input$ca_cto, if (length(blank) == 1) blank else NA)
})
# single-entity events (spikes) and steps: each file in the chosen time window
# a blank time box means the start or the end of the recording
ca_window = function(x) { d = x$data; t0 = input$ca_sfrom; t1 = input$ca_sto
  if (!isTRUE(is.finite(t0))) t0 = -Inf; if (!isTRUE(is.finite(t1))) t1 = Inf
  k = d$t_s >= t0 & d$t_s <= min(t1, max(d$t_s))
  validate(need(sum(k) > 20, paste(x$file, ": fewer than 20 points in the time window")))
  list(t = d$t_s[k], I = d$I_A[k]) }
ca_spk = reactive({ lapply(ca_x(), function(x) { w = ca_window(x)
  list(file = x$file, t = w$t, I = w$I, res = detect_spikes(w$t, w$I, input$ca_base, input$ca_k, input$ca_pol, input$ca_sep)) }) })
ca_est = reactive({ lapply(ca_x(), function(x) { w = ca_window(x)
  list(file = x$file, t = w$t, I = w$I, res = detect_steps(w$t, w$I, input$ca_ebin / 1000, input$ca_ewin, input$ca_ek, input$ca_edir)) }) })
ca_steps = reactive(do.call(rbind, lapply(ca_est(), function(z) if (nrow(z$res$steps)) data.frame(file = z$file, t_step_s = z$res$steps$t_step_s,
  height_pA = z$res$steps$height_A * 1e12, level_before_pA = z$res$steps$level_before_A * 1e12, level_after_pA = z$res$steps$level_after_A * 1e12))))
ca_events = reactive(do.call(rbind, lapply(ca_spk(), function(z) if (nrow(z$res$events)) cbind(file = z$file, z$res$events))))
output$ca_events = renderTable(if (input$ca_mode == "spikes") { e = ca_events()
  if (is.null(e)) data.frame(Note = "No events above the threshold.") else round_df(utils::head(e, 200)) } else if (input$ca_mode == "esteps") { s = ca_steps()
  if (is.null(s)) data.frame(Note = "No steps above the threshold.") else round_df(utils::head(s, 200)) })
output$ca_table = renderTable(round_df(ca_res()))
output$ca_fit = renderTable(if (input$ca_mode %in% c("multistep", "steps")) fmt_line(ca_line(), input$ca_unit))
ca_plot_obj = reactive({
  if (input$ca_mode %in% c("charge", "ume")) return(ca_traces_obj())
  if (input$ca_mode == "spikes") {
    z = ca_spk()[[1]]; e = z$res$events; u = current_unit(z$I)
    d = thin_trace(z$t, z$I * u$scale); b = thin_trace(z$t, z$res$baseline * u$scale, 2000)
    g = ggplot2::ggplot(d, ggplot2::aes(t, y)) + ggplot2::geom_line(linewidth = 0.3, colour = "grey30") +
      ggplot2::geom_line(data = b, colour = "#1A6FDF", linewidth = 0.6)
    if (nrow(e)) g = g + ggplot2::geom_point(data = data.frame(t = e$t_peak_s, y = z$I[findInterval(e$t_peak_s, z$t)] * u$scale), colour = "#C0392B", size = 1.8)
    return(g + ggplot2::labs(x = "Time (s)", y = paste0("Current (", u$label, ")"), subtitle = sprintf("%s: %d events (red), baseline (blue), threshold %s x noise", z$file, nrow(e), input$ca_k)) + g_theme())
  }
  if (input$ca_mode == "esteps") {
    z = ca_est()[[1]]; s = z$res$steps; u = current_unit(z$I); d = thin_trace(z$t, z$I * u$scale)
    g = ggplot2::ggplot(d, ggplot2::aes(t, y)) + ggplot2::geom_line(linewidth = 0.3, colour = "grey45")
    if (nrow(s)) g = g + ggplot2::geom_point(data = data.frame(t = s$t_step_s, y = s$level_after_A * u$scale), colour = "#C0392B", shape = 17, size = 2)
    return(g + ggplot2::labs(x = "Time (s)", y = paste0("Current (", u$label, ")"), subtitle = sprintf("%s: %d steps (red, at the level after each step), threshold %s x spread", z$file, nrow(s), input$ca_ek)) + g_theme())
  }
  if (input$ca_mode == "cottrell") {
    d = do.call(rbind, lapply(ca_x(), function(x) { k = x$data$t_s > 0
      data.frame(file = x$file, x = 1 / sqrt(x$data$t_s[k]), y = abs(x$data$I_A[k]) * 1e6) }))
    lo = 1 / sqrt(input$ca_tto); hi = 1 / sqrt(input$ca_tfrom)
    return(ggplot2::ggplot(d, ggplot2::aes(x, y, colour = file)) + ggplot2::geom_path(linewidth = 0.5) +
      ggplot2::annotate("rect", xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "#1A6FDF") +
      ggplot2::coord_cartesian(xlim = c(0, 2 * hi)) + ggplot2::scale_colour_manual(values = pal(length(unique(d$file))), name = NULL) +
      ggplot2::labs(x = expression(italic(t)^{-1/2}*" (s"^{-1/2}*")"), y = expression("|"*italic(i)*"| (µA)"), subtitle = "Shaded: fitted time window") + g_theme())
  }
  r = ca_res(); r$y = ca_y(); r$used = r$include & !is.na(r$conc) & r$conc > 0 & r$conc >= input$ca_cfrom & r$conc <= input$ca_cto
  cal_plot_obj(r, ca_line(), input$ca_unit, "Current at the end of the step (µA)")
})
output$ca_plot = renderPlot(ca_plot_obj())
# current traces, thinned for long recordings and in a unit that suits the current
ca_traces_obj = reactive({ xs = ca_x(); u = current_unit(unlist(lapply(xs, function(x) range(x$data$I_A))))
  d = do.call(rbind, lapply(xs, function(x) cbind(file = x$file, thin_trace(x$data$t_s, x$data$I_A * u$scale))))
  g = ggplot2::ggplot(d, ggplot2::aes(t, y, colour = file)) + ggplot2::geom_path(linewidth = 0.4) +
    ggplot2::scale_colour_manual(values = pal(length(xs)), name = NULL) + ggplot2::labs(x = "Time (s)", y = paste0("Current (", u$label, ")")) +
    g_theme() + ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
  if (labels_on() && identical(input$ca_mode, "multistep")) { r = range(d$y, na.rm = TRUE)
    g = g + step_labels(xs[[1]], r[2] + 0.06 * diff(r), diff(r), input$ca_avg, input$ca_skip) }
  g })
output$ca_traces = renderPlot(ca_traces_obj())
# Exploration of chronoamperometry files: zoom, event explorer, file check and comparison (R/explore.R)
ca_names = reactive(vapply(ca_x(), `[[`, "", "file"))
observe({ f = ca_names()
  updateSelectInput(session, "ca_zfile", choices = f, selected = if (isTRUE(input$ca_zfile %in% f)) input$ca_zfile else f[1])
  updateSelectInput(session, "ca_ref", choices = f, selected = if (isTRUE(input$ca_ref %in% f)) input$ca_ref else f[length(f)]) })
ca_zr = reactiveVal(NULL)
observeEvent(input$ca_zbrush, { b = input$ca_zbrush; if (!is.null(b) && b$xmax > b$xmin) ca_zr(c(b$xmin, b$xmax)) })
observeEvent(input$ca_zdbl, ca_zr(NULL)); observeEvent(input$ca_zfile, ca_zr(NULL), ignoreInit = TRUE)
ca_zx = reactive({ xs = ca_x(); i = match(input$ca_zfile, ca_names()); xs[[if (is.na(i)) 1 else i]] })
# steps or spikes of one file, as marks for the zoom, only when that analysis is chosen
ca_marks = function(file) tryCatch({
  if (input$ca_mode == "esteps") { z = Filter(function(z) z$file == file, ca_est()); s = z[[1]]$res$steps; data.frame(t = s$t_step_s, I = s$level_after_A, kind = "step") }
  else if (input$ca_mode == "spikes") { z = Filter(function(z) z$file == file, ca_spk())[[1]]; e = z$res$events
    data.frame(t = e$t_peak_s, I = z$I[findInterval(e$t_peak_s, z$t)], kind = "spike") }
  else NULL }, error = function(e) NULL)
ca_zoom_obj = reactive({
  x = ca_zx(); t = x$data$t_s; r = ca_zr(); if (is.null(r)) r = range(t)
  k = trace_window(t, r[1], r[2]); u = current_unit(x$data$I_A[k]); d = thin_trace(t[k], x$data$I_A[k] * u$scale, 4000)
  g = ggplot2::ggplot(d, ggplot2::aes(t, y)) + ggplot2::geom_line(linewidth = 0.3, colour = "grey30")
  m = ca_marks(x$file)
  if (!is.null(m) && nrow(m)) { m = m[m$t >= r[1] & m$t <= r[2], ]
    if (nrow(m)) g = g + ggplot2::geom_point(data = data.frame(t = m$t, y = m$I * u$scale), colour = "#C0392B",
                                             shape = if (m$kind[1] == "step") 17 else 16, size = 2) }
  g + ggplot2::coord_cartesian(xlim = r) + ggplot2::labs(x = "Time (s)", y = paste0("Current (", u$label, ")"),
    subtitle = sprintf("%s: %.6g to %.6g s, %s points%s", x$file, r[1], r[2], format(length(k), big.mark = ","),
                       if (length(k) <= 8000) ", every point shown" else ", drawn from the highest and lowest point of each short interval")) + g_theme()
})
output$ca_zoom = renderPlot(ca_zoom_obj())

# events of the chosen analysis, with the time window of each file
ca_xev = reactive({
  validate(need(input$ca_mode %in% c("spikes", "esteps"), "Choose Single-entity spikes or Single-entity steps (top of this tab) to explore the events."))
  # a file with no events adds no rows
  one = function(z, t, amp, dur, q) if (length(t)) data.frame(file = z$file, t0 = min(z$t), t1 = max(z$t), t = t, amp = amp, dur = dur, q = q)
  e = if (input$ca_mode == "spikes") lapply(ca_spk(), function(z) with(z$res$events, one(z, t_peak_s, amplitude_A, duration_s, charge_C)))
      else lapply(ca_est(), function(z) with(z$res$steps, one(z, t_step_s, height_A, rep(NA_real_, length(t_step_s)), rep(NA_real_, length(t_step_s)))))
  e = do.call(rbind, e)
  if (is.null(e)) data.frame(file = character(), t0 = numeric(), t1 = numeric(), t = numeric(), amp = numeric(), dur = numeric(), q = numeric()) else e
})
ca_xword = reactive(if (input$ca_mode == "spikes") "Spikes" else "Steps")
# legends only where they tell something: the file key with two or more files, the up/down key when both occur
ca_xlegend = function(d, shape = FALSE) ggplot2::guides(colour = if (length(unique(d$file)) > 1) "legend" else "none",
                                                        fill = if (length(unique(d$file)) > 1) "legend" else "none", shape = if (shape) "legend" else "none")
output$ca_xnote = renderUI({ e = ca_xev(); helpText(sprintf("%d %s in %d file(s).", nrow(e), tolower(ca_xword()), length(unique(e$file)))) })
ca_x_rate_obj = reactive({ e = ca_xev(); validate(need(nrow(e) > 0, "No events to show."))
  r = do.call(rbind, lapply(split(e, e$file), function(s) cbind(file = s$file[1], event_rate(s$t, s$t0[1], s$t1[1], input$ca_xbin))))
  ggplot2::ggplot(r, ggplot2::aes(t_mid_s, count, colour = file)) + ggplot2::geom_step(direction = "mid", linewidth = 0.6) +
    ggplot2::scale_colour_manual(values = pal(length(unique(r$file))), name = NULL) +
    ggplot2::labs(x = "Time (s)", y = sprintf("%s per %g s", ca_xword(), input$ca_xbin), subtitle = "How often events occur over time") + g_theme() + ca_xlegend(r) })
output$ca_x_rate = renderPlot(ca_x_rate_obj())
ca_x_amp_t_obj = reactive({ e = ca_xev(); validate(need(nrow(e) > 0, "No events to show.")); u = current_unit(e$amp)
  ggplot2::ggplot(e, ggplot2::aes(t, abs(amp) * u$scale, colour = file, shape = amp > 0)) + ggplot2::geom_point(size = 1.4, alpha = 0.8) +
    ggplot2::scale_y_log10() + ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 6), labels = c(`TRUE` = "up", `FALSE` = "down"), name = NULL) +
    ggplot2::scale_colour_manual(values = pal(length(unique(e$file))), name = NULL) +
    ggplot2::labs(x = "Time (s)", y = sprintf("|%s| (%s, log scale)", if (input$ca_mode == "spikes") "Amplitude" else "Height", u$label),
                  subtitle = "Size of each event over time") + g_theme() + ca_xlegend(e, shape = length(unique(e$amp > 0)) > 1) })
output$ca_x_amp_t = renderPlot(ca_x_amp_t_obj())
ca_x_amp_dur_obj = reactive({ e = ca_xev()
  validate(need(input$ca_mode == "spikes", "Durations are measured for spikes. Steps are lasting changes and have no duration."), need(nrow(e) > 0, "No events to show."))
  u = current_unit(e$amp)
  ggplot2::ggplot(e, ggplot2::aes(dur * 1e3, amp * u$scale, colour = file)) + ggplot2::geom_point(size = 1.4, alpha = 0.8) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() + ggplot2::scale_colour_manual(values = pal(length(unique(e$file))), name = NULL) +
    ggplot2::labs(x = "Duration at half height (ms, log scale)", y = sprintf("Amplitude (%s, log scale)", u$label), subtitle = "Amplitude against duration") + g_theme() + ca_xlegend(e) })
output$ca_x_amp_dur = renderPlot(ca_x_amp_dur_obj())
ca_x_hist_obj = reactive({ e = ca_xev(); validate(need(nrow(e) > 0, "No events to show.")); u = current_unit(e$amp)
  ggplot2::ggplot(e, ggplot2::aes(abs(amp) * u$scale, fill = file)) + ggplot2::geom_histogram(bins = 30, colour = "white", position = "identity", alpha = 0.7) +
    ggplot2::scale_x_log10() + ggplot2::scale_fill_manual(values = pal(length(unique(e$file))), name = NULL) +
    ggplot2::labs(x = sprintf("|%s| (%s, log scale)", if (input$ca_mode == "spikes") "Amplitude" else "Height", u$label), y = ca_xword(),
                  subtitle = "Distribution of event sizes") + g_theme() + ca_xlegend(e) })
output$ca_x_hist = renderPlot(ca_x_hist_obj())

# file check: currents shown in one unit that suits them
ca_check_df = reactive(do.call(rbind, lapply(ca_x(), file_check)))
# the file check as shown in the tab and the report: currents in one unit that suits them, named in the headers
ca_check_table = function() { f = ca_check_df(); u = current_unit(c(f$min_A, f$max_A)); L = function(x) sprintf("%s (%s)", x, u$label)
  o = data.frame(File = f$file, Points = format(f$points, big.mark = ","), `Duration (s)` = signif(f$duration_s, 6), `Sampling (Hz)` = signif(f$sampling_Hz, 6),
                 Min = signif(f$min_A * u$scale, 5), Max = signif(f$max_A * u$scale, 5), Noise = signif(f$noise_A * u$scale, 3),
                 `Net change` = signif(f$net_change_A * u$scale, 4), Gaps = f$gaps, `At limit (%)` = signif(f$at_limit_pct, 2),
                 `Potential (V)` = f$E_V, `Filter (Hz)` = f$filter_Hz, check.names = FALSE)
  names(o)[names(o) %in% c("Min", "Max", "Noise", "Net change")] = L(c("Min", "Max", "Noise", "Net change"))
  o = o[, colSums(!is.na(o)) > 0, drop = FALSE]
  # numbers as text with four significant figures, so small currents are not rounded to 0.00 by the table
  o[] = lapply(o, function(v) if (is.numeric(v)) vapply(v, function(x) if (is.na(x)) "" else format(signif(x, 4), scientific = abs(x) > 0 && (abs(x) < 1e-3 || abs(x) >= 1e7), drop0trailing = TRUE), "") else v)
  o }
output$ca_check = renderTable(ca_check_table())
output$ca_check_notes = renderUI({ f = ca_check_df(); n = file_check_notes(f)
  tagList(if (length(n)) div(class = "alert alert-warning", HTML(paste(htmltools::htmlEscape(n), collapse = "<br>"))) else
            helpText("No gaps in time, no clipping and no sampling problem found.")) })

# comparison with a reference file
ca_cmp = reactive({ xs = ca_x(); validate(need(length(xs) >= 2, "Load two or more files to compare them."))
  compare_traces(xs, ref = max(1, match(input$ca_ref, ca_names())), n_bins = 1000, start_zero = isTRUE(input$ca_cmp0)) })
ca_cmp_plot_obj = reactive({ d = ca_cmp(); u = current_unit(d$I_A)
  ggplot2::ggplot(d, ggplot2::aes(t_s, I_A * u$scale, colour = file)) + ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = pal(length(unique(d$file))), name = NULL) +
    ggplot2::labs(x = "Time (s)", y = paste0("Current (", u$label, ")"), subtitle = if (isTRUE(input$ca_cmp0)) "Each trace shifted to start at zero" else "Bin means") + g_theme() })
output$ca_cmp_plot = renderPlot(ca_cmp_plot_obj())
ca_cmp_diff_obj = reactive({ d = ca_cmp(); d = d[d$file != input$ca_ref, ]; u = current_unit(d$diff_A)
  ggplot2::ggplot(d, ggplot2::aes(t_s, diff_A * u$scale, colour = file)) + ggplot2::geom_hline(yintercept = 0, colour = "grey60") + ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = pal(length(unique(d$file))), name = NULL) +
    ggplot2::labs(x = "Time (s)", y = paste0("Difference (", u$label, ")"), subtitle = sprintf("Each file minus %s", input$ca_ref)) + g_theme() })
output$ca_cmp_diff = renderPlot(ca_cmp_diff_obj())

output$ca_csv = downloadHandler("chronoamperometry.csv", function(f) write.csv(if (input$ca_mode == "spikes") ca_events() else if (input$ca_mode == "esteps") ca_steps() else if (ca_eq()) ca_res() else add_line(ca_res(), ca_line()), f, row.names = FALSE))
output$ca_png = downloadHandler("chronoamperometry.png", function(f) fig_save(f, ca_plot_obj(), width = 6.5, height = 4.5, dpi = 300))

# Chronopotentiometry
cp_tbl = conc_table("cp_files", "cp_unit", "cp_conc")
cp_x = reactive(req(safe_read(req(input$cp_files))))
cp_res = reactive({
  r = do.call(rbind, lapply(cp_x(), cp_transition))
  t = req(cp_tbl()); r$conc = t$conc[match(r$file, t$file)]; r$include = t$include[match(r$file, t$file)]
  r[order(r$conc), ]
})
cp_area = reactive(area_of("cp"))
cp_sand = reactive({ r = cp_res(); r = r[r$include & !is.na(r$tau_s), ]
  if (!nrow(r)) return(NULL)
  sand_table(r, r$conc, input$cp_n, cp_area()) })
output$cp_table = renderTable(round_df(cp_res()[, c("file", "conc", "current_A", "run_s", "tau_s", "E_quarter_V", "step_mV", "note")]))
output$cp_sand = renderTable({
  s = req(cp_sand()); f = sand_fit(s, input$cp_n, cp_area())
  if (is.null(f)) return(round_df(s))
  data.frame(Files = f$n_files, `Slope (A s^1/2 per mM)` = sprintf("%.4g", f$slope_A_s05_per_mM),
             `Intercept (A s^1/2)` = sprintf("%.3g", f$intercept_A_s05), `R²` = sprintf("%.4f", f$r2),
             `D from slope (cm²/s)` = if (is.na(f$D_from_slope_cm2_s)) "give the electrode size" else sprintf("%.3g", f$D_from_slope_cm2_s),
             check.names = FALSE) })
cp_plot_obj = reactive({
  xs = cp_x(); r = cp_res()
  g = traces_plot(xs, "t_s", "E_V", "Time (s)", "Potential (V)", 1)
  m = r[!is.na(r$tau_s), ]
  if (nrow(m)) g = g + ggplot2::geom_point(data = m, ggplot2::aes(tau_s, E_at_tau_V), inherit.aes = FALSE, size = 2.6, shape = 4) +
    ggplot2::labs(subtitle = "Crosses: transition times")
  g + ggplot2::coord_cartesian(xlim = c(0, min(max(r$run_s), 4 * max(c(m$tau_s, 5), na.rm = TRUE))))
})
output$cp_plot = renderPlot(cp_plot_obj())
cp_sand_plot_obj = reactive({
  s = req(cp_sand())
  ggplot2::ggplot(s, ggplot2::aes(conc_mM, i_sqrt_tau * 1e6)) + ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#1A6FDF") +
    ggplot2::labs(x = sprintf("Concentration (%s)", input$cp_unit), y = expression(italic(i)*tau^{1/2}*" (µA s"^{1/2}*")")) + g_theme()
})
output$cp_sand_plot = renderPlot(cp_sand_plot_obj())
output$cp_csv = downloadHandler("chronopotentiometry.csv", function(f) {
  s = cp_sand(); write.csv(if (is.null(s)) cp_res() else merge(cp_res(), s, all = TRUE), f, row.names = FALSE) })
output$cp_png = downloadHandler("chronopotentiometry.png", function(f) fig_save(f, cp_plot_obj(), width = 7, height = 4.5, dpi = 300))

# RRDE
rr_cal = reactive({
  xs = req(safe_read(req(input$rr_cal_files)))
  r = do.call(rbind, lapply(xs, function(x) rrde_collection(x, c(input$rr_cal_from, input$rr_cal_to), rpm_from_name(x$file))))
  r[order(r$rpm), ]
})
output$rr_cal_table = renderTable(round_df(rr_cal()))
output$rr_levich = renderTable({ l = levich(rr_cal()); if (is.null(l)) return(NULL)
  data.frame(`Levich: rotation rates` = l$n, `Slope (A per (rad/s)^1/2)` = sprintf("%.4g", l$slope_A_per_rad_s), `R²` = sprintf("%.4f", l$r2), check.names = FALSE) })
rr_cal_plot_obj = reactive({
  xs = safe_read(input$rr_cal_files)
  d = do.call(rbind, lapply(xs, function(x) rbind(data.frame(file = x$file, electrode = "Disk", E_V = x$data$E_V, I = x$data$I_A * 1e6),
                                                  data.frame(file = x$file, electrode = "Ring", E_V = x$data$E_V, I = x$data$I2_A * 1e6))))
  ggplot2::ggplot(d, ggplot2::aes(E_V, I, colour = file)) + ggplot2::geom_path(linewidth = 0.5) +
    ggplot2::facet_wrap(~ electrode, scales = "free_y", ncol = 1) +
    ggplot2::annotate("rect", xmin = input$rr_cal_from, xmax = input$rr_cal_to, ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "#1A6FDF") +
    ggplot2::scale_colour_manual(values = pal(length(xs)), name = NULL) +
    ggplot2::labs(x = "Disk potential (V)", y = "Current (µA)", subtitle = sprintf("Shaded: window for N. Mean N = %.3f", mean(rr_cal()$N))) + g_theme()
})
output$rr_cal_plot = renderPlot(rr_cal_plot_obj())
rr_N = reactive(if (!is.na(input$rr_N)) input$rr_N else mean(rr_cal()$N))
rr_fe_fits = reactive({
  xs = req(safe_read(req(input$rr_fe_files)))
  lapply(xs, rrde_fe, N = rr_N(), disk_window = c(input$rr_fe_from, input$rr_fe_to),
         ring_baseline_window = if (isTRUE(input$rr_use_base)) c(input$rr_base_from, input$rr_base_to))
})
output$rr_fe_table = renderTable(round_df(do.call(rbind, lapply(rr_fe_fits(), `[[`, "summary"))))
rr_fe_plot_obj = reactive({
  d = do.call(rbind, lapply(rr_fe_fits(), function(z) cbind(file = z$summary$file, z$curve[, c("E_V", "FE_pct")])))
  d = d[is.finite(d$FE_pct) & d$FE_pct < 200, ]
  ggplot2::ggplot(d, ggplot2::aes(E_V, FE_pct, colour = file)) + ggplot2::geom_path(linewidth = 0.5) +
    ggplot2::annotate("rect", xmin = input$rr_fe_from, xmax = input$rr_fe_to, ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "#1A6FDF") +
    ggplot2::scale_colour_manual(values = pal(length(unique(d$file))), name = NULL) +
    ggplot2::labs(x = "Disk potential (V)", y = "Faradaic efficiency (%)", subtitle = sprintf("N = %.3f. Shaded: averaging window", rr_N())) + g_theme()
})
output$rr_fe_plot = renderPlot(rr_fe_plot_obj())
output$rr_csv = downloadHandler("rrde.csv", function(f) {
  a = tryCatch(rr_cal(), error = function(e) NULL); b = tryCatch(do.call(rbind, lapply(rr_fe_fits(), `[[`, "summary")), error = function(e) NULL)
  out = if (is.null(a)) b else if (is.null(b)) a else merge(a, b, all = TRUE)
  write.csv(out, f, row.names = FALSE) })

# RRDE: electron number and peroxide yield over the same disk window as the faradaic efficiency
output$rr_n_table = renderTable(round_df(do.call(rbind, lapply(rr_fe_fits(), function(z) {
  d = z$curve; w = d[d$E_V >= min(input$rr_fe_from, input$rr_fe_to) & d$E_V <= max(input$rr_fe_from, input$rr_fe_to), ]
  v = rrde_n_peroxide(w$I_A, w$I2_A - z$summary$ring_baseline_A, z$summary$N)
  data.frame(file = z$summary$file, N = z$summary$N, window = z$summary$window, n_mean = mean(v$n), n_sd = stats::sd(v$n),
             H2O2_mean_pct = mean(v$H2O2_pct), H2O2_sd_pct = stats::sd(v$H2O2_pct))
}))))

# Double-layer capacitance
dl_res = reactive(cdl_from_cvs(req(safe_read(req(input$dl_files))), input$dl_E, input$dl_Cs))
output$dl_points = renderTable(round_df(dl_res()$points))
output$dl_fit = renderTable({ f = dl_res()$fit
  if (is.null(f)) return(data.frame(Note = "At least 2 scan rates are needed."))
  round_df(cbind(f, roughness_factor = f$ECSA_cm2 / area_of("dl"))) })
dl_plot_obj = reactive({ p = dl_res()$points
  ggplot2::ggplot(p, ggplot2::aes(v_V_s, I_charging_A * 1e6)) + ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#1A6FDF") +
    ggplot2::labs(x = "Scan rate (V/s)", y = "Charging current (µA)", subtitle = sprintf("Read at %g V", input$dl_E)) + g_theme() })
output$dl_plot = renderPlot(dl_plot_obj())
output$dl_csv = downloadHandler("double_layer.csv", function(f) write.csv(cbind(dl_res()$points, dl_res()$fit), f, row.names = FALSE))

# Chronocoulometry (Anson)
cc_qdl = reactive({ if (is.null(input$cc_blank)) return(NA)
  b = safe_read(input$cc_blank)[[1]]; f = anson_fit(b, input$cc_tfrom, input$cc_tto, 1, 1)
  if (is.null(f)) NA else f$intercept_C })
cc_res = reactive({
  xs = req(safe_read(req(input$cc_files))); A = area_of("cc")
  do.call(rbind, lapply(xs, anson_fit, t_from = input$cc_tfrom, t_to = input$cc_tto, n = input$cc_n, C_mM = input$cc_C,
                        D_cm2_s = if (input$cc_solve == "area") input$cc_D else NA, area_cm2 = if (input$cc_solve == "D") A else NA,
                        Q_dl_C = cc_qdl()))
})
output$cc_table = renderTable(round_df(cbind(cc_res(), Q_dl_blank_C = cc_qdl())))
cc_plot_obj = reactive({
  xs = req(safe_read(req(input$cc_files)))
  d = do.call(rbind, lapply(xs, function(x) { q = if (!is.null(x$data$Q_C)) x$data$Q_C else
      c(0, cumsum(diff(x$data$t_s) * (x$data$I_A[-1] + x$data$I_A[-nrow(x$data)]) / 2))
    data.frame(file = x$file, x = sqrt(x$data$t_s), y = abs(q) * 1e6) }))
  ggplot2::ggplot(d, ggplot2::aes(x, y, colour = file)) + ggplot2::geom_path(linewidth = 0.5) +
    ggplot2::annotate("rect", xmin = sqrt(input$cc_tfrom), xmax = sqrt(input$cc_tto), ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "#1A6FDF") +
    ggplot2::scale_colour_manual(values = pal(length(xs)), name = NULL) +
    ggplot2::labs(x = expression(italic(t)^{1/2}*" (s"^{1/2}*")"), y = expression("|"*italic(Q)*"| (µC)"), subtitle = "Anson plot. Shaded: fitted window") + g_theme() })
output$cc_plot = renderPlot(cc_plot_obj())
output$cc_csv = downloadHandler("chronocoulometry.csv", function(f) write.csv(cc_res(), f, row.names = FALSE))

# Rotating disk: Levich from the limiting current, Koutecky-Levich at chosen potentials
rd_x = reactive(req(safe_read(req(input$rd_files))))
rd_lim = reactive({ xs = rd_x()
  r = do.call(rbind, lapply(xs, function(x) { d = x$data; k = d$E_V >= min(input$rd_from, input$rd_to) & d$E_V <= max(input$rd_from, input$rd_to)
    data.frame(file = x$file, rpm = rpm_from_name(x$file), I_lim_A = mean(d$I_A[k])) }))
  r[order(r$rpm), ] })
rd_lev = reactive({ r = rd_lim(); r = r[!is.na(r$rpm), ]; if (nrow(r) < 2) return(NULL)
  levich_fit(r$rpm, r$I_lim_A, input$rd_n, area_of("rd"), input$rd_C, input$rd_nu, if (input$rd_solve == "n") input$rd_D else NA) })
output$rd_lev = renderTable({ l = rd_lev(); if (is.null(l)) return(data.frame(Note = "At least 2 files with a rotation rate in the name are needed."))
  round_df(cbind(l, area_cm2 = area_of("rd"))) })
rd_lev_plot_obj = reactive({ r = rd_lim()
  ggplot2::ggplot(r, ggplot2::aes(sqrt(rpm_to_omega(rpm)), abs(I_lim_A) * 1e6)) + ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#1A6FDF") +
    ggplot2::labs(x = expression(omega^{1/2}*" ((rad/s)"^{1/2}*")"), y = expression("|"*italic(i)[L]*"| (µA)")) + g_theme() })
output$rd_lev_plot = renderPlot(rd_lev_plot_obj())
rd_kl = reactive({ xs = rd_x(); E = suppressWarnings(as.numeric(strsplit(input$rd_Ekl, "[,; ]+")[[1]])); E = E[!is.na(E)]
  rpm = vapply(xs, function(x) rpm_from_name(x$file), 0); if (sum(!is.na(rpm)) < 2 || !length(E)) return(NULL)
  pts = do.call(rbind, lapply(E, function(e) data.frame(E_V = e, rpm = rpm, I_A = vapply(xs, function(x)
    stats::approx(x$data$E_V, x$data$I_A, e, ties = mean)$y, 0))))
  fits = do.call(rbind, lapply(E, function(e) { p = pts[pts$E_V == e & !is.na(pts$rpm), ]
    cbind(E_V = e, koutecky_levich(p$rpm, p$I_A, if (input$rd_solve == "D") input$rd_n else NA, area_of("rd"), input$rd_C, input$rd_nu,
                                   if (input$rd_solve == "n") input$rd_D else NA)) }))
  list(points = pts, fits = fits) })
output$rd_kl = renderTable({ k = rd_kl(); if (is.null(k)) return(data.frame(Note = "Give at least one potential and two rotation rates.")); round_df(k$fits) })
rd_kl_plot_obj = reactive({ k = req(rd_kl()); p = k$points[!is.na(k$points$rpm), ]
  ggplot2::ggplot(p, ggplot2::aes(1 / sqrt(rpm_to_omega(rpm)), 1 / abs(I_A) / 1e3, colour = factor(E_V))) + ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, fullrange = TRUE) + ggplot2::expand_limits(x = 0) +
    ggplot2::scale_colour_manual(values = pal(length(unique(p$E_V))), name = "E (V)") +
    ggplot2::labs(x = expression(omega^{-1/2}*" ((rad/s)"^{-1/2}*")"), y = expression("1/|"*italic(i)*"| (mA"^{-1}*")")) + g_theme() })
output$rd_kl_plot = renderPlot(rd_kl_plot_obj())
rd_traces_obj = reactive(traces_plot(rd_x(), "E_V", "I_A", "Potential (V)", "Current (µA)"))
output$rd_traces = renderPlot(rd_traces_obj())
output$rd_csv = downloadHandler("rotating_disk.csv", function(f) {
  l = tryCatch(rd_lev(), error = function(e) NULL); k = tryCatch(rd_kl()$fits, error = function(e) NULL)
  write.csv(rd_lim(), f, row.names = FALSE)
  if (!is.null(l)) suppressWarnings(write.table(cbind(analysis = "Levich", l), f, sep = ",", append = TRUE, row.names = FALSE))
  if (!is.null(k)) suppressWarnings(write.table(cbind(analysis = "Koutecky-Levich", k), f, sep = ",", append = TRUE, row.names = FALSE)) })

# Impedance
ei_e = reactive({ f = req(input$ei_file); r = tryCatch(read_eis(f$datapath, f$name, input$instrument %||% "auto"), error = function(e) {
  notify_error(e, duration = 15); NULL })
  if (!is.null(r) && nzchar(r$note)) showNotification(r$note, duration = 10)
  req(r) })
ei_d = reactive(ei_e()$data)
ei_zu = reactive(if (identical(input$ei_zunit, "auto")) auto_ohm(ei_d()$Z) else input$ei_zunit)
ei_zs = function(z) z / unit_factor(ei_zu()); ei_zl = function() unit_label(ei_zu())
ei_fs = function(f) f / unit_factor(input$ei_funit)
ei_ph = function(z) { a = Arg(z); if (identical(input$ei_phunit, "rad")) a else a * 180 / pi }
ei_phl = function() if (identical(input$ei_phunit, "rad")) "radians" else "degrees"

output$ei_valid = renderTable({ v = eis_validate(ei_e())
  data.frame(Check = c("Points", "Frequency range", "Repeated frequencies", "Rows left out", "Share of points with Z'' < 0",
                       "|Z| in the file vs from Z', Z'' (largest difference)", "Phase in the file vs from Z', Z'' (largest difference)", "Flags", "Note"),
             Result = c(v$points, v$freq_range, v$repeated_frequencies, v$missing_or_nonfinite, sprintf("%.0f%%", 100 * v$share_capacitive),
                        if (is.null(v$max_mod_diff_pct)) "no |Z| column" else sprintf("%.3f%%", v$max_mod_diff_pct),
                        if (is.null(v$max_phase_diff_deg)) "no phase column" else sprintf("%.2f degrees", v$max_phase_diff_deg), v$flags,
                        if (nzchar(v$note)) v$note else "none")) })
ei_ft = reactive(eis_features(ei_d()))
output$ei_feat = renderTable({ f = ei_ft()$features
  data.frame(Feature = c("Real part at the highest frequency", "Phase at the lowest frequency", "Slope of log|Z| over the lowest decade",
                         "Most negative phase", "at", "Maxima of -Z''"),
             Value = c(sprintf("%.4g ohm", f$high_f_real_ohm), sprintf("%.1f degrees", f$low_f_phase_deg), sprintf("%.2f", f$low_f_logZ_slope),
                       sprintf("%.1f degrees", f$min_phase_deg), sprintf("%.4g Hz", f$f_at_min_phase_Hz), f$local_maxima_of_minus_Zim)) })
output$ei_hints = renderUI(tags$ul(lapply(ei_ft()$hints, tags$li)))
nyq = function(meas, fit = NULL, sub = NULL) {
  g = ggplot2::ggplot() + ggplot2::geom_point(data = meas, ggplot2::aes(re, -im), size = 1.8)
  if (!is.null(fit)) g = g + ggplot2::geom_path(data = fit, ggplot2::aes(re, -im), colour = "#F14040", linewidth = 0.7)
  g + ggplot2::coord_equal() + ggplot2::labs(x = sprintf("Z′ (%s)", ei_zl()), y = sprintf("−Z″ (%s)", ei_zl()), subtitle = sub) + g_theme()
}
ei_meas = reactive({ d = ei_d(); data.frame(f = ei_fs(d$freq_Hz), re = Re(ei_zs(d$Z)), im = Im(ei_zs(d$Z)), mod = Mod(ei_zs(d$Z)), ph = ei_ph(d$Z)) })
ei_nyq_raw_obj = reactive(nyq(ei_meas(), sub = "Data as read"))
output$ei_nyq_raw = renderPlot(ei_nyq_raw_obj())
ei_kkr = reactive(lin_kk(ei_d()))
output$ei_kk = renderTable({ k = ei_kkr()
  data.frame(`RC elements (M)` = k$M, mu = sprintf("%.3f", k$mu), `Largest residual (% of |Z|)` = sprintf("%.2f", k$max_abs_res_pct),
             `RMS residual (% of |Z|)` = sprintf("%.2f", k$rms_res_pct), check.names = FALSE) })
res_plot = function(r, sub) {
  d = rbind(data.frame(f = ei_fs(r$freq_Hz), res = r$res_re_pct, part = "real"), data.frame(f = ei_fs(r$freq_Hz), res = r$res_im_pct, part = "imaginary"))
  ggplot2::ggplot(d, ggplot2::aes(f, res, colour = part)) + ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_point(size = 1.8) + ggplot2::geom_line(linewidth = 0.4) + ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = c(real = "#1A6FDF", imaginary = "#F14040"), name = NULL) +
    ggplot2::labs(x = sprintf("Frequency (%s)", input$ei_funit), y = "Residual (% of |Z|)", subtitle = sub) + g_theme() }
ei_kk_plot_obj = reactive(res_plot(ei_kkr()$residuals, "Linear Kramers-Kronig test (Schonleber et al. 2014)"))
output$ei_kk_plot = renderPlot(ei_kk_plot_obj())
ei_circs = reactive({
  own = trimws(strsplit(input$ei_custom %||% "", ";")[[1]]); own = own[nzchar(own)]
  ok = vapply(own, function(s) !is.null(tryCatch(eis_compile(s), error = function(e) NULL)), TRUE); bad = own[!ok]
  if (length(bad)) showNotification(paste("Not a valid circuit:", paste(bad, collapse = ", ")), type = "error", duration = 12)
  unique(c(input$ei_circuits, setdiff(own, bad))) })
ei_all = eventReactive(input$ei_run, {
  d = ei_d(); cs = req(ei_circs())
  withProgress(message = "Fitting circuits", value = 0, {
    fits = lapply(seq_along(cs), function(i) { incProgress(1 / length(cs), detail = cs[i])
      tryCatch(eis_fit(d, cs[i], input$ei_weight, input$ei_start, input$ei_bounds), error = function(e) {
        showNotification(sprintf("%s: %s", cs[i], conditionMessage(e)), type = "warning", duration = 10); NULL }) })
  })
  names(fits) = cs; fits = fits[!vapply(fits, is.null, TRUE)]
  rank = do.call(rbind, lapply(fits, function(f) data.frame(circuit = f$circuit, notation = f$notation, parameters = f$n_params,
    chi2 = f$chi2, chi2_red = f$chi2_red, nrmse_pct = f$nrmse_pct, AIC = f$aic, AICc = f$aicc, BIC = f$bic,
    worst_rel_se_pct = suppressWarnings(max(f$params$rel_se_pct, na.rm = TRUE)), at_bound = any(f$params$at_bound), converged = f$converged)))
  rank = rank[order(rank$AICc), ]; rank$dAICc = rank$AICc - min(rank$AICc)
  rank$flag = ifelse(!rank$converged, "optimizer did not converge",
               ifelse(!is.finite(rank$worst_rel_se_pct) | rank$worst_rel_se_pct > 100, "a parameter is not determined (error above 100%)",
               ifelse(rank$at_bound, "a parameter is at a bound", "")))
  rank$recommended = FALSE; ok = rank$flag == ""; if (any(ok)) rank$recommended[which(ok)[1]] = TRUE
  list(rank = rank, fits = fits, settings = list(file = input$ei_file$name, weighting = input$ei_weight, start = input$ei_start,
       bounds = input$ei_bounds, circuits = cs, when = format(Sys.time(), "%Y-%m-%d %H:%M"), points = nrow(d)))
})
observeEvent(ei_all(), { r = ei_all()$rank
  updateSelectInput(session, "ei_show", choices = r$circuit, selected = if (any(r$recommended)) r$circuit[r$recommended] else r$circuit[1]) })
output$ei_rank = renderTable({ r = ei_all()$rank
  r$recommended = ifelse(r$recommended, "recommended", "")
  round_df(r[, c("recommended", "circuit", "notation", "parameters", "AICc", "dAICc", "BIC", "chi2_red", "nrmse_pct", "worst_rel_se_pct", "flag")]) })
ei_sel = reactive(ei_all()$fits[[req(input$ei_show)]])
output$ei_fitstats = renderTable({ f = ei_sel()
  data.frame(Statistic = c("Weighting", "Weighted chi-squared", "Reduced chi-squared", "RMSE", "Normalised RMSE", "AIC", "AICc", "BIC", "Converged"),
             Value = c(f$weighting, sprintf("%.4g", f$chi2), sprintf("%.4g", f$chi2_red), sprintf("%.4g ohm", f$rmse_ohm),
                       sprintf("%.2f%% of |Z|", f$nrmse_pct), sprintf("%.1f", f$aic), sprintf("%.1f", f$aicc), sprintf("%.1f", f$bic),
                       if (f$converged) "yes" else "no")) })
output$ei_params = renderTable({ p = ei_sel()$params
  # capacitances in the chosen unit; everything else in base units
  isC = p$unit == "F"; k = unit_factor(input$ei_cunit)
  data.frame(Parameter = p$parameter, Value = ifelse(isC, sprintf("%.4g", p$value / k), sprintf("%.4g", p$value)),
             `Standard error` = ifelse(isC, sprintf("%.3g", p$se / k), sprintf("%.3g", p$se)),
             Unit = ifelse(isC, unit_label(input$ei_cunit), p$unit), `Error (%)` = sprintf("%.1f", p$rel_se_pct),
             `At bound` = ifelse(p$at_bound, "yes", ""), check.names = FALSE) })
output$ei_derived = renderTable({
  v = eis_derived(ei_sel(), area_of("ei"), input$ei_C, input$ei_n, T_K("ei"))
  if (!length(v)) return(data.frame(Note = "No derived quantities for this circuit."))
  k = unit_factor(input$ei_cunit)
  data.frame(Quantity = names(v), Value = vapply(names(v), function(n) if (grepl("_F$", n)) sprintf("%.4g %s", v[[n]] / k, unit_label(input$ei_cunit)) else sprintf("%.4g", v[[n]]), ""),
             Note = ifelse(grepl("^k0|^D_from", names(v)), "equal concentrations of both forms, at the formal potential",
                    ifelse(grepl("HsuMansfeld", names(v)), "C = (Q R^(1-n))^(1/n), R in parallel with the CPE",
                    ifelse(grepl("Brug", names(v)), "C = Q^(1/n) (Rs R/(Rs + R))^((1-n)/n)", ifelse(grepl("_area_", names(v)), "resistance x geometric area", ""))))) })
ei_curves = reactive({ d = ei_d(); f = ei_sel()
  list(meas = ei_meas(), fit = data.frame(f = ei_fs(d$freq_Hz), re = Re(ei_zs(f$fit_Z)), im = Im(ei_zs(f$fit_Z)), mod = Mod(ei_zs(f$fit_Z)), ph = ei_ph(f$fit_Z))) })
ei_nyq_obj = reactive({ k = ei_curves(); nyq(k$meas, k$fit, sprintf("Points: data. Line: %s", ei_sel()$notation)) })
bode = function(what, ylab, logy) { k = ei_curves()
  g = ggplot2::ggplot() + ggplot2::geom_point(data = k$meas, ggplot2::aes(f, .data[[what]]), size = 1.8) +
    ggplot2::geom_line(data = k$fit, ggplot2::aes(f, .data[[what]]), colour = "#F14040", linewidth = 0.7) + ggplot2::scale_x_log10() +
    ggplot2::labs(x = sprintf("Frequency (%s)", input$ei_funit), y = ylab) + g_theme()
  if (logy) g + ggplot2::scale_y_log10() else g }
ei_mag_obj = reactive(bode("mod", sprintf("|Z| (%s)", ei_zl()), TRUE))
ei_ph_obj = reactive(bode("ph", sprintf("Phase (%s)", ei_phl()), FALSE))
ei_res_obj = reactive(res_plot(ei_sel()$residuals, sprintf("Fit residuals, %s", ei_sel()$notation)))
ei_diag_obj = reactive({ g = eis_diagram(ei_sel()$circuit)
  ggplot2::ggplot() + ggplot2::geom_segment(data = g$wires, ggplot2::aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.6) +
    ggplot2::geom_rect(data = g$boxes, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white", colour = "black") +
    ggplot2::geom_text(data = g$boxes, ggplot2::aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label), size = 4) +
    ggplot2::coord_equal() + ggplot2::theme_void() + ggplot2::labs(subtitle = ei_sel()$notation) })
output$ei_nyq = renderPlot(ei_nyq_obj()); output$ei_bode_mag = renderPlot(ei_mag_obj()); output$ei_bode_ph = renderPlot(ei_ph_obj())
output$ei_res_plot = renderPlot(ei_res_obj()); output$ei_diagram = renderPlot(ei_diag_obj())
output$ei_interp = renderUI({ t = eis_interpret(ei_sel(), ei_kkr(), ei_ft())
  r = ei_all()$rank; close = r$circuit[r$dAICc < 2 & r$circuit != ei_sel()$circuit]
  if (length(close)) t = c(t, sprintf("Within 2 AICc units of the best circuit, so not separated by the data: %s.", paste(close, collapse = ", ")))
  tagList(h6(sprintf("Circuit: %s", ei_sel()$notation)), tags$ul(lapply(t, tags$li))) })
output$ei_csv = downloadHandler("eis_fit.csv", function(f) {
  a = ei_all(); write.csv(do.call(rbind, lapply(a$fits, function(z) cbind(circuit = z$circuit, notation = z$notation, weighting = z$weighting, z$params))), f, row.names = FALSE) })
output$ei_settings = downloadHandler("eis_fit_settings.txt", function(f) { s = ei_all()$settings
  writeLines(c("ElectroFast impedance fit", sprintf("File: %s (%d frequencies)", s$file, s$points), sprintf("Fitted: %s", s$when),
               sprintf("Weighting: %s", s$weighting), sprintf("Starting values given: %s", if (nzchar(s$start)) s$start else "none (from the spectrum)"),
               sprintf("Bounds given: %s", if (nzchar(s$bounds)) s$bounds else "none (CPE exponents 0 to 1, other parameters positive)"),
               "Optimizer: Nelder-Mead then BFGS on log parameters, 12 starts", sprintf("Circuits: %s", paste(s$circuits, collapse = "; ")),
               "", "Ranking:", utils::capture.output(print(ei_all()$rank[, c("circuit", "AICc", "dAICc", "flag", "recommended")], row.names = FALSE))), f) })
output$ei_png = downloadHandler("eis_figures.zip", function(f) {
  dd = file.path(tempdir(), "eis_png"); unlink(dd, recursive = TRUE); dir.create(dd)
  fig_save(file.path(dd, "nyquist.png"), ei_nyq_obj(), width = 6, height = 5, dpi = 300)
  fig_save(file.path(dd, "bode_magnitude.png"), ei_mag_obj(), width = 6.5, height = 4, dpi = 300)
  fig_save(file.path(dd, "bode_phase.png"), ei_ph_obj(), width = 6.5, height = 4, dpi = 300)
  fig_save(file.path(dd, "fit_residuals.png"), ei_res_obj(), width = 6.5, height = 4, dpi = 300)
  fig_save(file.path(dd, "kramers_kronig.png"), res_plot(ei_kkr()$residuals, "Linear Kramers-Kronig test"), width = 6.5, height = 4, dpi = 300)
  fig_save(file.path(dd, "circuit.png"), ei_diag_obj(), width = 6, height = 2.5, dpi = 300)
  # the zip package needs no zip program on the computer; utils::zip is the fallback
  if (requireNamespace("zip", quietly = TRUE)) zip::zipr(f, list.files(dd, full.names = TRUE)) else utils::zip(f, list.files(dd, full.names = TRUE), flags = "-j9X") })

# Equations and sources
# Equations grouped by technique: the quantity on the left, the typeset formula (equation_display, rendered by
# MathJax) on its own line, its conditions below in smaller text, and the source last in small muted text
output$eq_table = renderUI({ d = equation_sources(); x = equation_display
  groups = unique(d$Technique)
  withMathJax(div(class = "eq-list",
    tags$style(HTML(".eq-list { max-width: 1100px; }
      .eq-list h6 { margin: 1.6rem 0 0.4rem; padding-bottom: 0.3rem; border-bottom: 2px solid #1A6FDF; text-transform: uppercase; letter-spacing: 0.04em; }
      .eq-row { display: grid; grid-template-columns: 15rem 1fr; gap: 0 1.5rem; padding: 0.8rem 0; border-bottom: 1px solid #E3E7EB; }
      .eq-q { font-weight: 600; }
      .eq-f { font-size: 1.08rem; overflow-x: auto; }
      .eq-n { color: #3d4650; margin-top: 0.25rem; }
      .eq-s { color: #6b7682; font-size: 0.82rem; margin-top: 0.35rem; }")),
    lapply(groups, function(gname) tagList(h6(gname),
      lapply(which(d$Technique == gname), function(i) div(class = "eq-row",
        div(class = "eq-q", d$Quantity[i]),
        div(if (nzchar(x$formula[i])) div(class = "eq-f", HTML(paste0("\\(\\displaystyle ", x$formula[i], "\\)"))),
            if (nzchar(x$note[i])) div(class = "eq-n", HTML(x$note[i])),
            div(class = "eq-s", d$Source[i])))))))) })

# Technique tabs for the HTML report and the slide deck: one section per tab that has data,
# each with a few lines of text, its plots and its tables. A tab without data is left out.
tp = function(f) tryCatch(f(), error = function(e) NULL)
tp_section = function(title, text = NULL, plots = list(), tables = list()) {
  plots = Filter(Negate(is.null), lapply(plots, tp)); tables = Filter(function(x) !is.null(x) && NROW(x) > 0, lapply(tables, tp))
  if (!length(plots) && !length(tables)) return(NULL)
  list(title = title, text = tp(text), plots = plots, tables = tables)
}
technique_parts = function() Filter(Negate(is.null), list(
  tp_section("Surface area (CV peak)", function() sprintf("Charge density %g µC/cm², window %g to %g V, last %s sweep.", input$sa_q, input$sa_from, input$sa_to, input$sa_dir),
             list(Voltammograms = sa_plot_obj), list(Results = function() round_df(sa_res()))),
  tp_section("Double-layer capacitance", function() sprintf("Charging current read at %g V, with C_s = %g µF/cm².", input$dl_E, input$dl_Cs),
             list(`Charging current against scan rate` = dl_plot_obj), list(Fit = function() round_df(dl_res()$fit), Points = function() round_df(dl_res()$points))),
  tp_section("CV: peaks, scan rate, kinetics", function() scan_rate_reading(sr_rate()),
             list(Voltammograms = sr_cv_plot_obj, `Peak currents against v and √v` = sr_lin_plot_obj, `Scan rate (log-log)` = sr_plot_obj,
                  `Peak current against √v (Randles–Ševčík)` = sr_rs_plot_obj),
             list(Peaks = function() round_df(sr_pk()[, c("file", "v_V_s", "Epa_V", "ipa_A", "Epc_V", "ipc_A", "dEp_mV", "regime", "k0_nicholson_cm_s", "k0_klingler_kochi_cm_s", "k0_K_table_cm_s", "k0_note")]),
                  `Scan-rate fits` = function() round_df(sr_rate()$fits), `Log-log slopes` = function() round_df(sr_rate()$loglog),
                  `Area, D or coverage` = function() round_df(sr_area_tab()))),
  tp_section("Tafel", function() sprintf("Fitted window %g to %g V.", input$tf_from, input$tf_to), list(`log I against E` = tf_plot_obj), list(`Chosen window` = function() round_df(tf_one()))),
  tp_section("DPV / SWV calibration", function() paste(unlist(fmt_line(pv_line(), input$pv_unit)), collapse = "; "),
             list(Calibration = pv_plot_obj, Voltammograms = pv_traces_obj), list(Points = function() round_df(pv_pts()[, c("file", "conc", "E_V", "I_A", "used")]))),
  tp_section("Open-circuit potential", NULL, list(Traces = oc_plot_obj), list(Results = function() round_df(oc_res()))),
  tp_section("Chronoamperometry", function() c(if (input$ca_mode %in% c("multistep", "steps")) paste(unlist(fmt_line(ca_line(), input$ca_unit)), collapse = "; ") else
               sprintf("Analysis: %s.", c(charge = "charge passed", cottrell = "Cottrell", ume = "microelectrode steady state",
                                          spikes = "single-entity spikes", esteps = "single-entity steps")[[input$ca_mode]]),
               tryCatch(file_check_notes(ca_check_df()), error = function(e) NULL)),
             list(Results = ca_plot_obj, `Current traces` = ca_traces_obj, `Events over time` = ca_x_rate_obj, `Event size over time` = ca_x_amp_t_obj,
                  `Amplitude against duration` = ca_x_amp_dur_obj, `Event sizes` = ca_x_hist_obj),
             c(list(Results = function() round_df(ca_res()), `File check` = ca_check_table),
               if (isTRUE(input$ca_mode == "spikes")) list(Events = function() round_df(utils::head(ca_events(), 200))),
               if (isTRUE(input$ca_mode == "esteps")) list(Steps = function() round_df(utils::head(ca_steps(), 200))))),
  tp_section("Chronocoulometry", NULL, list(`Anson plot` = cc_plot_obj), list(`Anson fit` = function() round_df(cc_res()))),
  tp_section("Chronopotentiometry", NULL, list(`Potential traces` = cp_plot_obj, `Sand line` = cp_sand_plot_obj),
             list(`Transition times` = function() round_df(cp_res()[, c("file", "conc", "current_A", "tau_s", "E_quarter_V", "note")]))),
  tp_section("Rotating disk", NULL, list(Voltammograms = rd_traces_obj, Levich = rd_lev_plot_obj, `Koutecky-Levich` = rd_kl_plot_obj),
             list(Levich = function() round_df(rd_lev()), `Koutecky-Levich` = function() round_df(rd_kl()$fits))),
  tp_section("RRDE", NULL, list(`Collection efficiency` = rr_cal_plot_obj, `Faradaic efficiency` = rr_fe_plot_obj),
             list(`Collection efficiency` = function() round_df(rr_cal()), `Faradaic efficiency` = function() round_df(do.call(rbind, lapply(rr_fe_fits(), `[[`, "summary"))))),
  tp_section("Impedance (EIS)", function() eis_interpret(ei_sel(), ei_kkr(), ei_ft()),
             list(`Circuit` = ei_diag_obj, Nyquist = ei_nyq_obj, `Bode magnitude` = ei_mag_obj, `Bode phase` = ei_ph_obj, `Fit residuals` = ei_res_obj, `Kramers-Kronig residuals` = ei_kk_plot_obj),
             list(`Circuit ranking` = function() round_df(ei_all()$rank[, c("circuit", "parameters", "AICc", "dAICc", "nrmse_pct", "worst_rel_se_pct", "flag", "recommended")]),
                  Parameters = function() round_df(ei_sel()$params)))))
