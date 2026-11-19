# Server side of the Sensor performance menu (UI in R/app_performance_ui.R, calculations in R/performance.R).
# Sourced inside the server function with local = TRUE, after app_methods_server.R, whose g_theme, round_df and
# tp_section it uses.

perf_read = function(id) {
  f = req(input[[id]])
  tryCatch(utils::read.csv(f$datapath, check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE),
           error = function(e) { notify_error(e, duration = 12); req(FALSE) })
}
# a calculation error becomes a message, and the tab stays empty
perf_try = function(f) tryCatch(f(), error = function(e) { notify_error(e, duration = 12); req(FALSE) })
perf_col = "#1A6FDF"
for (k in names(perf_templates)) local({ k = k
  id = c(sensitivity = "ps", selectivity = "pl", specificity = "pc", stability = "pt", precision = "pp")[[k]]
  output[[paste0(id, "_template")]] = downloadHandler(sprintf("%s_template.csv", k), function(f) utils::write.csv(perf_templates[[k]], f, row.names = FALSE))
})
lim = function(v) if (is.null(v) || is.na(v)) NULL else v

# Sensitivity
ps_res = reactive(perf_try(function() perf_sensitivity(perf_read("ps_file"), lim(input$ps_from) %||% -Inf, lim(input$ps_to) %||% Inf, lim(input$ps_area) %||% NA)))
ps_tab = reactive({ f = ps_res()$fit; cu = input$ps_cunit; su = input$ps_sunit
  data.frame(Standards = f$n, `Linear range` = sprintf("%g to %g %s", f$from, f$to, cu),
             Sensitivity = sprintf("%.4g ± %.2g %s/%s", f$slope, f$slope_se, su, cu), Intercept = sprintf("%.4g %s", f$intercept, su),
             `R²` = sprintf("%.4f", f$r2), LOD = sprintf("%.3g %s", f$lod, cu), LOQ = sprintf("%.3g %s", f$loq, cu),
             `Noise s from` = if (f$blanks >= 2) sprintf("%d blanks", f$blanks) else "residual SD of the line",
             `Sensitivity per area` = if (is.na(f$slope_per_area)) "" else sprintf("%.4g %s %s⁻¹ cm⁻²", f$slope_per_area, su, cu), check.names = FALSE) })
ps_plot_obj = reactive({ r = ps_res(); f = r$fit
  ggplot2::ggplot(r$points, ggplot2::aes(concentration, signal)) + ggplot2::geom_point(ggplot2::aes(shape = used), size = 2.6) +
    ggplot2::geom_abline(intercept = f$intercept, slope = f$slope, colour = perf_col) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "In the fit", `FALSE` = "Not in the fit"), name = NULL) +
    ggplot2::labs(x = sprintf("Concentration (%s)", input$ps_cunit), y = sprintf("Signal (%s)", input$ps_sunit),
                  subtitle = sprintf("Sensitivity %.4g %s/%s, R² %.4f, LOD %.3g %s", f$slope, input$ps_sunit, input$ps_cunit, f$r2, f$lod, input$ps_cunit)) + g_theme() })
output$ps_plot = renderPlot(ps_plot_obj())
output$ps_table = renderTable(ps_tab())
output$ps_csv = downloadHandler("sensitivity.csv", function(f) utils::write.csv(ps_res()$fit, f, row.names = FALSE))

# Selectivity
pl_res = reactive(perf_try(function() perf_selectivity(perf_read("pl_file"), input$pl_tol %||% 5)))
pl_plot_obj = reactive({ d = pl_res(); d$lab = trimws(paste(d$interferent, d$interferent_concentration)); tol = input$pl_tol %||% 5
  ggplot2::ggplot(d, ggplot2::aes(stats::reorder(lab, abs(change_pct)), change_pct)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -tol, ymax = tol, fill = "#E8EEF6") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") + ggplot2::geom_col(width = 0.6, fill = perf_col) +
    ggplot2::coord_flip() + ggplot2::labs(x = NULL, y = "Change of the analyte signal (%)", subtitle = sprintf("Shaded: within ±%g%%", tol)) + g_theme() })
output$pl_plot = renderPlot(pl_plot_obj())
pl_tab = reactive({ d = pl_res(); d$within_tolerance = ifelse(d$within_tolerance, "yes", "no"); round_df(d) })
output$pl_table = renderTable(pl_tab())
output$pl_csv = downloadHandler("selectivity.csv", function(f) utils::write.csv(pl_res(), f, row.names = FALSE))

# Specificity
pc_res = reactive(perf_try(function() perf_specificity(perf_read("pc_file"))))
pc_plot_obj = reactive({ d = pc_res()
  ggplot2::ggplot(d, ggplot2::aes(stats::reorder(species, cross_reactivity_pct), cross_reactivity_pct, fill = target)) +
    ggplot2::geom_col(width = 0.6) + ggplot2::scale_fill_manual(values = c(`TRUE` = perf_col, `FALSE` = "#8C96A3"), labels = c(`TRUE` = "Target", `FALSE` = "Other species"), name = NULL) +
    ggplot2::coord_flip() + ggplot2::labs(x = NULL, y = "Response relative to the target (%)") + g_theme() })
output$pc_plot = renderPlot(pc_plot_obj())
pc_tab = reactive({ d = pc_res(); d$target = ifelse(d$target, "yes", "no"); names(d)[names(d) == "signal"] = sprintf("Signal (%s)", input$pc_sunit)
  names(d)[names(d) == "selectivity_ratio"] = "Selectivity ratio (target / species)"; round_df(d) })
output$pc_table = renderTable(pc_tab())
output$pc_csv = downloadHandler("specificity.csv", function(f) utils::write.csv(pc_res(), f, row.names = FALSE))

# Stability
pt_res = reactive(perf_try(function() perf_stability(perf_read("pt_file"), input$pt_limit %||% 90)))
pt_plot_obj = reactive({ r = pt_res(); lim = input$pt_limit %||% 90
  ggplot2::ggplot(r$points, ggplot2::aes(time, retention_pct)) + ggplot2::geom_hline(yintercept = lim, linetype = "dashed", colour = "grey40") +
    ggplot2::geom_line(colour = perf_col) + ggplot2::geom_point(colour = perf_col, size = 2.6) +
    ggplot2::labs(x = sprintf("Time (%s)", input$pt_tunit), y = "Signal retained (% of the first)", subtitle = sprintf("Dashed: %g%% limit", lim)) + g_theme() })
output$pt_plot = renderPlot(pt_plot_obj())
pt_tab = reactive({ s = pt_res()$summary; u = input$pt_tunit
  data.frame(Measurements = s$times, `Last time` = sprintf("%g %s", s$last_time, u), `Retained at the last time` = sprintf("%.1f%%", s$retained_at_last_pct),
             `Change (straight line)` = sprintf("%.3g%% per %s", s$change_pct_per_unit, sub("s$", "", u)),
             `Below the limit from` = if (is.na(s$first_time_below_limit)) sprintf("not below %g%%", s$limit_pct) else sprintf("%.3g %s", s$first_time_below_limit, u), check.names = FALSE) })
output$pt_table = renderTable(pt_tab())
output$pt_csv = downloadHandler("stability.csv", function(f) utils::write.csv(pt_res()$points, f, row.names = FALSE))

# Speed
pq_x = reactive({ f = req(input$pq_file); perf_try(function() perf_trace(f$datapath, f$name, input$instrument %||% "auto", input$pq_channel)) })
observeEvent(input$pq_file, {
  x = tryCatch(perf_trace(input$pq_file$datapath, input$pq_file$name, input$instrument %||% "auto"), error = function(e) NULL)
  updateSelectInput(session, "pq_channel", choices = if (is.null(x$channels)) character(0) else x$channels)
})
pq_res = reactive({ x = pq_x(); ta = req(lim(input$pq_tadd))
  perf_try(function() perf_speed(x$t, x$y, ta, input$pq_pre %||% 2, input$pq_post %||% 10, input$pq_smooth %||% 0.2)) })
pq_plot_obj = reactive({ x = pq_x(); th = thin_trace(x$t, x$y * x$scale); d = data.frame(t = th$t, y = th$y)
  g = ggplot2::ggplot(d, ggplot2::aes(t, y)) + ggplot2::geom_line(linewidth = 0.4) + ggplot2::labs(x = "Time (s)", y = x$label) + g_theme()
  r = tryCatch(pq_res(), error = function(e) NULL)
  if (is.null(r)) return(g)
  g + ggplot2::geom_vline(xintercept = r$t_add_s, linetype = "dashed", colour = "grey40") +
    ggplot2::geom_hline(yintercept = c(r$baseline, r$final) * x$scale, colour = perf_col, linewidth = 0.4) +
    (if (!is.na(r$t90_s)) ggplot2::geom_vline(xintercept = r$t_add_s + r$t90_s, colour = perf_col)) +
    ggplot2::labs(subtitle = sprintf("Dashed: addition. Lines: baseline, final level and t90 (%s s)", format(signif(r$t90_s, 3)))) })
output$pq_plot = renderPlot(pq_plot_obj())
pq_tab = reactive({ r = pq_res(); x = pq_x(); u = sub("^.*\\((.*)\\)$", "\\1", x$label)
  data.frame(File = x$file, Channel = x$channel %||% "", `Addition at` = sprintf("%g s", r$t_add_s),
             Baseline = sprintf("%.4g %s", r$baseline * x$scale, u), `Final level` = sprintf("%.4g %s", r$final * x$scale, u),
             Change = sprintf("%.4g %s", r$change * x$scale, u), t90 = sprintf("%.3g s", r$t90_s), t95 = sprintf("%.3g s", r$t95_s), check.names = FALSE) })
output$pq_table = renderTable(pq_tab())
output$pq_csv = downloadHandler("speed.csv", function(f) utils::write.csv(cbind(file = pq_x()$file, pq_res()), f, row.names = FALSE))

# Precision
pp_res = reactive(perf_try(function() perf_precision(perf_read("pp_file"))))
pp_plot_obj = reactive({ d = perf_cols(perf_read("pp_file")); d$value = perf_num(d$value); d = d[!is.na(d$value), , drop = FALSE]
  x = if ("day" %in% names(d)) "day" else NULL; col = if ("electrode" %in% names(d)) "electrode" else NULL
  d$xv = if (is.null(x)) "All" else as.character(d[[x]]); d$cv = if (is.null(col)) "All" else as.character(d[[col]])
  g = ggplot2::ggplot(d, ggplot2::aes(xv, value, colour = cv)) + ggplot2::geom_hline(yintercept = mean(d$value), colour = "grey60") +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 1), size = 2.2) +
    ggplot2::labs(x = if (is.null(x)) NULL else "Day", y = sprintf("Value (%s)", input$pp_sunit), colour = if (is.null(col)) NULL else "Electrode") + g_theme()
  if (is.null(col)) g + ggplot2::scale_colour_manual(values = perf_col, guide = "none") else g + ggplot2::scale_colour_manual(values = pal(length(unique(d$cv)))) })
output$pp_plot = renderPlot(pp_plot_obj())
pp_tab = reactive({ r = pp_res(); d = r$components; names(d)[2:3] = c(sprintf("SD (%s)", input$pp_sunit), "RSD (%)"); round_df(d) })
output$pp_table = renderTable(pp_tab())
output$pp_notes = renderUI({ r = pp_res(); helpText(sprintf("%d values, mean %.4g %s.", r$n, r$mean, input$pp_sunit), if (length(r$notes)) tags$ul(lapply(r$notes, tags$li))) })
output$pp_csv = downloadHandler("precision.csv", function(f) utils::write.csv(pp_res()$components, f, row.names = FALSE))

# Report sections: one per entry with data
performance_parts = function() Filter(Negate(is.null), list(
  tp_section("Sensitivity", NULL, list(Calibration = ps_plot_obj), list(Results = ps_tab)),
  tp_section("Selectivity", function() sprintf("Change of the target signal with each interferent added (Thévenot et al., Pure Appl. Chem. 1999, 71, 2333, second procedure). Tolerated change ±%g%%.", input$pl_tol),
             list(`Change of the target signal` = pl_plot_obj), list(Interferents = pl_tab)),
  tp_section("Specificity", function() c("Response to each species alone at the target concentration (Thévenot et al., Pure Appl. Chem. 1999, 71, 2333, first procedure). It is given as a percentage of the target response (cross-reactivity) and as the ratio of the target signal to that of the species.",
                                          "Note: IUPAC recommends the term selectivity, which can be graded, and discourages specificity for the same concept. A method is either specific or not, and few methods, if any, are specific (Vessman et al., Pure Appl. Chem. 2001, 73, 1381)."),
             list(`Cross-reactivity` = pc_plot_obj), list(Species = pc_tab)),
  tp_section("Stability", NULL, list(Retention = pt_plot_obj), list(Summary = pt_tab, Points = function() round_df(pt_res()$points))),
  tp_section("Speed", NULL, list(Response = pq_plot_obj), list(Results = pq_tab)),
  tp_section("Precision", function() c(sprintf("%d values, mean %.4g %s.", pp_res()$n, pp_res()$mean, input$pp_sunit), pp_res()$notes),
             list(Values = pp_plot_obj), list(Components = pp_tab))))
