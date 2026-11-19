# ElectroFast: Shiny app for electrochemical data from CHI, PalmSens, EMF16 and Pine instruments.
# Modules: ECSA of gold (CV in H2SO4), ISE calibration (EMF16, one or more days), ISE
# selectivity (fixed interference and separate solution methods), precision, and the
# Techniques menu: surface area from any CV peak, scan rate, Tafel, DPV/SWV calibration,
# open-circuit potential, chronoamperometry (i-t, steps, multi-potential steps),
# chronopotentiometry (transition time, Sand) and RRDE (collection and faradaic efficiency).
# Run from the project folder: shiny::runApp()

suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(bslib))
# In the browser (Shinylive), officer::read_pptx() hangs unless xml2 is attached; attaching it costs nothing on the desktop
suppressPackageStartupMessages(library(xml2))

options(shiny.maxRequestSize = 300 * 1024^2)
source("R/read_chi_txt.R")
source("R/ecsa.R")
source("R/app_core.R")
source("R/read_pine_txt.R")
source("R/read_pssession.R")
source("R/read_other.R")
source("R/read_mpr.R")
source("R/read_abf.R")
source("R/methods.R")
source("R/electrochem.R")
source("R/rs_quasi_table.R")
source("R/provenance.R")
source("R/qc.R")
source("R/spikes.R")
source("R/steps.R")
source("R/explore.R")
source("R/app_guide.R")
source("R/app_theme.R")
source("R/app_panels.R")
source("R/diagnose.R")
source("R/eis.R")
source("R/report_export.R")
source("R/app_electrode.R")
source("R/app_methods_ui.R")
source("R/cell_diagram.R")
source("R/performance.R")
source("R/app_performance_ui.R")
source("R/table_names.R")

letters_from = function(txt) {
  x = toupper(trimws(strsplit(txt, "[,; ]+")[[1]]))
  x[nzchar(x)]
}
day_colours = function(days) setNames(rep(cycle_palette("#CC9900"), length.out = length(days)), days)
# numbers as text with 4 significant figures, so tables never round small values such as 1e-05 M to 0.00, and readable headers
round_df = function(d, digits = 4) {
  if (is.null(d)) return(d)   # an empty table (for example QC charts before a history is loaded) stays empty
  d[] = lapply(d, function(v) if (!is.numeric(v)) v else
    vapply(v, function(x) if (is.na(x)) NA_character_ else format(signif(x, digits), scientific = x != 0 && abs(x) < 1e-3), ""))
  names(d) = nice_names(names(d))   # readable headers for display; the data keep their code names
  d
}
std_cols = c("file", "day", "letter", "conc_M", "include", "note")

render_std = function(tbl) {
  DT::renderDT(DT::datatable(tbl()[, std_cols], rownames = FALSE, selection = "none",
                             editable = list(target = "cell", disable = list(columns = c(0, 2, 5))),
                             options = list(pageLength = 50, dom = "t")))
}
apply_edit = function(tbl, info) {
  d = tbl()
  col = std_cols[info$col + 1]
  # typing a concentration includes the standard, so one edit is enough
  if (col == "conc_M") {
    v = suppressWarnings(as.numeric(info$value))
    d$conc_M[info$row] = v
    d$include[info$row] = !is.na(v) && v > 0
    if (!is.na(v)) d$note[info$row] = "concentration entered by hand"
  }
  if (col == "day") d$day[info$row] = trimws(info$value)
  if (col == "include") d$include[info$row] = toupper(trimws(info$value)) %in% c("TRUE", "T", "1", "YES")
  tbl(d)
}

ui = page_navbar(
  title = electrofast_brand(), window_title = "ElectroFast",
  id = "main_nav",
  theme = electrofast_theme(),
  header = tags$head(tags$script(src = "session_keeper.js"), tags$link(rel = "icon", type = "image/png", href = "favicon.png")),
  methods_ui(),
  performance_ui(),
  nav_panel("QC charts",
    layout_sidebar(
      sidebar = sidebar(width = 340,
        fileInput("qc_file", "Load a QC history (CSV)", accept = ".csv"),
        downloadButton("qc_template", "Empty history (CSV)", class = "btn-sm mb-2"),
        actionButton("qc_add_session", "Add this session's results", class = "btn-sm mb-2"),
        helpText("Adds one row per day for each result that is loaded. These are the ISE slope and intercept (ISE calibration tab),",
                 "the background b (background-corrected fit) and the ECSA (ECSA tab)."),
        h6("Or add one value"),
        dateInput("qc_date", "Date", Sys.Date()),
        textInput("qc_quantity", "Quantity", "Check standard X (mV)"),
        numericInput("qc_value", "Value", NA),
        textInput("qc_electrode", "Electrode (optional)", ""),
        actionButton("qc_add_one", "Add", class = "btn-sm mb-2"),
        numericInput("qc_nref", "Reference period: first n sessions of each quantity (blank = all)", NA, min = 2, step = 1),
        helpText("Set the reference period to sessions known to be good. Left blank, every session (any bad ones too) sets the limits."),
        helpText("Limits are the mean and SD of the reference period. The rules (Westgard et al., Clin. Chem. 1981, 27, 493) are as thus.",
                 "1-3s: one point beyond 3 SD. 2-2s: two in a row beyond 2 SD on one side. R-4s: two in a row more than 4 SD apart.",
                 "4-1s: four in a row beyond 1 SD on one side. 10x: ten in a row on one side. One point beyond 2 SD is a warning.",
                 "The checks suggested for a broken rule are prompts, not diagnoses."),
        helpText("Nothing is stored on the server: download the history and load it again next time."),
        downloadButton("qc_download", "Download the history (CSV)")
      ),
      uiOutput("qc_status"),
      plotOutput("qc_plot", height = "560px", fill = FALSE),
      h6("Limits"), tableOutput("qc_limits"),
      h6("Sessions"), tableOutput("qc_table")
    )
  ),
  panel_cell(),
  nav_panel("Report",
    card(
      h5("Analysis report"),
      p("One HTML, PDF or PowerPoint report with every tab that has data loaded, the Techniques tabs included. It holds the settings, plots and tables, and a short summary of each result.",
        "Tabs without data are left out."),
      textInput("rep_title", "Report title", "ElectroFast analysis report"),
      textInput("rep_author", "Prepared by (optional)", ""),
      downloadButton("report_html", "Download report (HTML)"),
      downloadButton("report_pdf", "Download report (PDF)"),
      downloadButton("report_pptx", "Download slides (PowerPoint)"),
      downloadButton("report_bundle", "Download reproducibility bundle (zip)"),
      helpText("The PDF and the slides hold the same results as the HTML report. The slides are 16:9. They start with a title slide, which lists the electrolyte and electrodes typed on the Cell diagram tab, and a slide of main points. One slide per plot follows, with its result underneath, then the tables. The settings and files come last. The HTML report has a contents list and keeps long tables on one page."),
      helpText("Every report records the ElectroFast version, the R version, when it was made and a checksum of each file used. The bundle adds every setting, the package versions, every results table as CSV and a copy of the code. The analysis can then be traced and repeated.")
    )
  ),
  guide_ui(),
  panel_equations(),
  nav_panel("About",
    card(
      h5("What ElectroFast does"),
      p("ElectroFast analyses electrochemical data straight from the files the instruments write. Load the files on a tab,",
        "set the options, and read the results and plots next to them. The Report tab collects every tab with data into one report."),
      p("The Techniques menu has one entry per technique, named as in potentiostat software. Each entry holds the analyses",
        "that technique allows, for example surface area, peaks and kinetics under cyclic voltammetry."),
      p("The Sensor performance menu reports the sensitivity, selectivity, specificity, stability, speed and precision of a",
        "sensor from the results of any technique."),
      p("The Guide tab explains what each tab needs and what it gives. The Equations tab lists every equation and its source."),
      h5("Files"),
      p("ElectroFast reads files from CH Instruments, PalmSens, BioLogic, Gamry, Pine, Axon pCLAMP and Lawson Labs EMF16,",
        "and any other text or CSV file with a header row. Units are read from each file, never assumed."),
      p("The BioLogic and Gamry readers were checked on public test files from the yadg and gamry-parser projects.",
        "BioLogic binary (.mpr) files are not read: export the run from EC-Lab as text (.mpt)."),
      h5("Your data"),
      p("Files are used only while the page is open and are not kept on the server. A copy stays in this browser so that",
        "a page closed by mistake can be restored."),
      h5("Version"),
      p(sprintf("ElectroFast %s (%s). Every report records this version, the R version, the time and a checksum of each file.", electrofast_version, electrofast_version_date))
    )
  ),
  nav_spacer(),
  nav_item(div(class = "d-flex align-items-center gap-2", style = "padding-top: 0.6rem;",
    tags$label("Instrument", `for` = "instrument", class = "small fw-semibold mb-0", style = "color: #5B6470;"),
    selectInput("instrument", NULL, width = "250px",
      c("Auto-detect from the file" = "auto", "CH Instruments (.txt)" = "chi", "PalmSens (.pssession)" = "pssession",
        "BioLogic EC-Lab (.mpt)" = "biologic", "Gamry (.DTA)" = "gamry", "Pine AfterMath (.txt)" = "pine", "Axon pCLAMP (.abf)" = "abf",
        "Metrohm Autolab, other text or CSV" = "text")),
    # each user decides; untick for figures with no added text (for a paper, for example)
    div(style = "padding-left: 0.6rem;", checkboxInput("fig_credit", "Name on figures", TRUE)),
    # off by default, so figures stay as they are unless the user asks for the labels
    div(style = "padding-left: 0.6rem;", checkboxInput("fig_labels", "Label plots", FALSE))))
)

server = function(input, output, session) {

  # A plot on a hidden tab can report a size of 0, and redrawing it at 0 px stops it with
  # "invalid 'height' argument". Every plot below falls back to a default size instead.
  size_or = function(what, default) function() {
    info = getCurrentOutputInfo()
    v = if (is.function(info[[what]])) info[[what]]() else NULL
    if (is.null(v) || !is.finite(v) || v <= 0) default else v
  }
  renderPlot = function(expr, env = parent.frame()) shiny::renderPlot(substitute(expr), env = env, quoted = TRUE,
                                                                    width = size_or("width", 700), height = size_or("height", 450))

  # Techniques menu (surface area, scan rate, Tafel, DPV/SWV, OCP, CA, CP, RRDE)
  source("app_methods_server.R", local = TRUE)
  source("app_performance_server.R", local = TRUE)

  # ECSA
  ecsa_x = reactive({ req(input$ecsa_file); read_echem(input$ecsa_file$datapath, input$ecsa_file$name, input$instrument) })
  ecsa_geo = electrode_server("ecsa_el")
  ecsa_fit = reactive(ecsa_gold(ecsa_x(), q_ref_uC_cm2 = input$ecsa_qref,
                                geometric_area_cm2 = ecsa_geo()$area_cm2))
  ecsa_plot = reactive(plot_ecsa_cv(ecsa_fit(), ecsa_x(), title = input$ecsa_title, reference = input$ecsa_ref) + watermark())
  output$ecsa_plot = renderPlot(ecsa_plot())
  ecsa_res = reactive({
    r = ecsa_fit()$result
    data.frame(File = input$ecsa_file$name, `Cathodic sweep` = r$segment, `Peak potential (V)` = r$peak_E_V,
               `Baseline window (V)` = sprintf("%.3f to %.3f", r$lower_limit_V, r$upper_limit_V),
               `Q (µC)` = r$Q_uC, `ECSA (cm²)` = r$ECSA_cm2, `Roughness factor` = r$roughness_factor,
               `Charge change from previous cycle (%)` = r$change_last_cycle_pct, Note = r$note %||% "", check.names = FALSE)
  })
  observeEvent(ecsa_fit(), { n = ecsa_fit()$result$note; if (isTRUE(nzchar(n))) showNotification(n, type = "warning", duration = 15, id = "ecsa_inverted") })
  output$ecsa_table = renderTable(ecsa_res(), digits = 4)
  output$ecsa_png = downloadHandler("ecsa_cv.png", function(f) fig_save(f, ecsa_plot(), width = 7, height = 5, dpi = 300))
  output$ecsa_pdf = downloadHandler("ecsa_cv.pdf", function(f) fig_save(f, ecsa_plot(), width = 7, height = 5, device = echem_pdf))
  output$ecsa_csv = downloadHandler("ecsa_results.csv", function(f) write.csv(ecsa_res(), f, row.names = FALSE))

  # ISE calibration
  cal_std = reactiveVal(NULL)
  # the coding is read from the file names of each new upload; the user can still change it
  # files come only from the user's upload; nothing is bundled with the app
  cal_input = reactiveVal(NULL)
  observeEvent(input$cal_files, cal_input(input$cal_files))
  observeEvent(cal_input(), updateRadioButtons(session, "scheme", selected = detect_scheme(cal_input()$name)))
  output$sheet_template = downloadHandler("sample_sheet_template.csv",
    function(f) utils::write.csv(sample_sheet_template(), f, row.names = FALSE))
  cal_sheet = reactive({
    if (is.null(input$cal_sheet)) return(NULL)
    tryCatch(read_sample_sheet(input$cal_sheet$datapath), error = function(e) {
      notify_error(e, duration = 10); NULL })
  })
  observeEvent(list(cal_input(), input$scheme, input$cal_sheet), {
    st = standards_table(req(cal_input()), scheme = input$scheme, sheet = cal_sheet())
    n_unknown = sum(is.na(st$conc_M))
    if (n_unknown) showNotification(sprintf("%d file(s) need a concentration: type it in the conc_M cell of the Standards table, or upload a sample sheet.", n_unknown),
                                    type = "warning", duration = 12)
    cal_std(st)
  })
  observeEvent(input$scheme, {
    if (input$scheme == "protocol") {
      updateTextInput(session, "fit_letters", value = "E, F, G, H, I, J, K, L, M")
      updateTextInput(session, "cmp_letters", value = "G, H, I, J, K, L, M")
      updateTextInput(session, "chk_letters", value = "I, X, Z")
      # 1e-9 to 1e-4 M on the plateau and 1e-2 to 1e-1 M on the line, as in the July defaults
      updateTextInput(session, "fim_flat", value = "A, B, C, D, E, F, G")
      updateTextInput(session, "fim_nernst", value = "K, L, M")
    } else {
      updateTextInput(session, "fit_letters", value = "E, F, G, H, I")
      updateTextInput(session, "cmp_letters", value = "F, G, H, I")
      updateTextInput(session, "chk_letters", value = "")
      updateTextInput(session, "fim_flat", value = "A, B, C, D, E, F")
      updateTextInput(session, "fim_nernst", value = "H, I")
    }
  }, ignoreInit = TRUE)
  output$cal_std = render_std(reactive(req(cal_std())))
  observeEvent(input$cal_std_cell_edit, apply_edit(cal_std, input$cal_std_cell_edit))
  cal_pts = reactive({
    st = req(cal_std())
    # letter codes choose the fitted standards by letter, the other routes by the concentration range
    st$fit = fit_flags(st, letters_from(input$fit_letters), conc_range = c(input$fit_cfrom %||% 1e-5, input$fit_cto %||% 0.1))
    withProgress(message = "Reading EMF16 files", add_emf(st, input$read_s, input$bgI, input$zA))
  })
  cal_res = reactive(calibrate(cal_pts(), input$xaxis, input$temp, input$zA))
  cal_cmp = reactive({
    L = letters_from(input$cmp_letters)
    if (!length(L)) return(NULL)
    p = cal_pts(); p$fit = p$letter %in% L
    calibrate(p, input$xaxis, input$temp, input$zA)
  })
  cal_bg = reactive(if (isTRUE(input$bg_fit)) fit_background(cal_pts(), input$xaxis, input$temp, input$zA))
  cal_plot = reactive({
    p = cal_pts(); r = req(cal_res())
    x = input$xaxis
    lines = do.call(rbind, lapply(seq_len(nrow(r$per_day)), function(i) {
      d = p[p$day == r$per_day$day[i] & p$fit, ]
      xr = range(d[[x]]) + c(-0.2, 0.2)
      data.frame(day = r$per_day$day[i], x = xr, E_mV = r$per_day$intercept_mV[i] + r$per_day$slope_mV_dec[i] * xr)
    }))
    s = r$summary
    lab = sprintf("Slope %.2f ± %.2f mV/decade (n = %d days)\n%.1f%% of %.2f, R² %.4f to %.4f",
                  s$slope_mean, ifelse(is.na(s$slope_sd), 0, s$slope_sd), s$days, s$pct_nernst, r$S, s$r2_min, s$r2_max)
    bg = cal_bg()
    curves = if (!is.null(bg)) background_curve(bg, range(p[[x]]))
    g0 = ggplot2::ggplot(p, ggplot2::aes(.data[[x]], E_mV, colour = day))
    if (!is.null(curves)) g0 = g0 + ggplot2::geom_line(data = curves, ggplot2::aes(xval, E_mV), linetype = "dashed", linewidth = 0.45)
    g0 +
      ggplot2::geom_line(data = lines, ggplot2::aes(x, E_mV), linewidth = 0.6) +
      ggplot2::geom_point(ggplot2::aes(shape = fit), size = 2.6, stroke = 0.8) +
      ggplot2::scale_shape_manual(values = c(`TRUE` = 15, `FALSE` = 0), labels = c(`TRUE` = "Fitted", `FALSE` = "Not fitted"), name = NULL) +
      ggplot2::scale_colour_manual(values = day_colours(sort(unique(p$day))), name = "Day") +
      ggplot2::annotate("text", x = max(p[[x]]), y = min(p$E_mV), hjust = 1, vjust = 0, size = 3.4, label = lab) +
      ggplot2::scale_x_continuous(breaks = seq(0, ceiling(max(p[[x]])), 1)) +
      ggplot2::labs(x = if (x == "neg_log_a") expression("−log "*italic(a)) else expression("−log "*italic(c)),
                    y = "EMF (mV)", subtitle = sprintf("Potential at %g s in each standard, %g °C.
Solid: linear fit over the fitted standards.%s",
                    input$read_s, input$temp, if (!is.null(cal_bg())) " Dashed: background-corrected fit to all standards." else "")) +
      ggplot2::theme_classic(base_size = 12) + watermark()
  })
  output$cal_plot = renderPlot(cal_plot())
  fmt_summary = function(r) {
    if (is.null(r)) return(NULL)
    s = r$summary
    data.frame(Days = s$days, `Linear range` = s$range, `Slope (mV/decade)` = sprintf("%.2f ± %.2f", s$slope_mean, s$slope_sd),
               `% of Nernstian` = sprintf("%.1f", s$pct_nernst), `R²` = sprintf("%.4f to %.4f", s$r2_min, s$r2_max), check.names = FALSE)
  }
  output$cal_summary = renderTable(fmt_summary(cal_res()))
  output$cal_perday = renderTable(round_df(req(cal_res())$per_day[, c("day", "n", "from_M", "to_M", "slope_mV_dec", "r2", "pct_nernst")]))
  output$cal_cmp = renderTable(fmt_summary(cal_cmp()))
  output$cal_bg = renderTable({
    bg = req(cal_bg()); s = bg$summary
    data.frame(Days = s$days, `Slope (mV/decade)` = sprintf("%.2f ± %.2f", s$slope_mean, s$slope_sd),
               `% of Nernstian` = sprintf("%.1f", s$pct_nernst),
               `Background b (M)` = sprintf("%.2g ± %.2g", s$b_bkg_mean, s$b_bkg_sd),
               `R²` = sprintf("%.4f to %.4f", s$r2_min, s$r2_max), check.names = FALSE)
  })
  output$cal_bg_perday = renderTable(round_df(req(cal_bg())$per_day))
  output$cal_bg_csv = downloadHandler("ise_background_fit.csv", function(f) write.csv(req(cal_bg())$per_day, f, row.names = FALSE))
  output$cal_steps = renderTable(round_df(steps_table(cal_pts(), input$xaxis, input$temp, input$zA)))
  observe(updateSelectInput(session, "trace_day", choices = sort(unique(req(cal_std())$day))))
  output$trace_plot = renderPlot({
    st = req(cal_std()); st = st[st$include & st$day == req(input$trace_day), ]
    tr = do.call(rbind, lapply(seq_len(nrow(st)), function(i) cbind(label = st$letter[i], emf_trace(st$datapath[i]))))
    ggplot2::ggplot(tr, ggplot2::aes(time_s, emf_mV, colour = label)) +
      ggplot2::geom_line(linewidth = 0.5) + ggplot2::geom_vline(xintercept = input$read_s, linetype = "dashed") +
      ggplot2::scale_colour_manual(values = setNames(rep(cycle_palette("#CC9900"), length.out = nrow(st)), st$letter), name = "Standard") +
      ggplot2::labs(x = "Time (s)", y = "EMF (mV)", subtitle = "Dashed line: reading time") +
      ggplot2::theme_classic(base_size = 12) + watermark()
  })
  output$cal_csv = downloadHandler("ise_calibration.csv", function(f) {
    r = cal_res(); p = cal_pts()
    out = merge(p[, c("day", "letter", "conc_M", "neg_log_c", "neg_log_a", "E_mV", "fit")], r$per_day, by = "day", all.x = TRUE)
    write.csv(out, f, row.names = FALSE)
  })
  output$cal_png = downloadHandler("ise_calibration.png", function(f) fig_save(f, cal_plot(), width = 7.5, height = 5, dpi = 300))
  output$cal_pdf = downloadHandler("ise_calibration.pdf", function(f) fig_save(f, cal_plot(), width = 7.5, height = 5, device = echem_pdf))

  # Reading time
  rt_settle = reactive(withProgress(message = "Checking settling", settling_table(cal_pts(), input$rt_drift)))
  # possible causes: the signs measured on this tab, day by day, with the user's cut-offs
  dx = reactive({
    th = list(slope_lo = input$dx_slope_lo %||% 90, slope_hi = 200 - (input$dx_slope_lo %||% 90), low_end = input$dx_low %||% 0.6,
              background = input$dx_bg %||% 3, drift_mV = input$dx_drift %||% 5, r2 = input$dx_r2 %||% 0.98, unsettled = input$dx_unsettled %||% 0.2)
    diagnose_ise(cal_pts(), req(cal_res()), safely(cal_bg), safely(rt_settle), input$xaxis, input$read_s, th) })
  output$dx_signs = renderTable(round_df(dx()$signs))
  output$dx_causes = renderTable(dx()$causes)
  rt_slopes = reactive(withProgress(message = "Fitting at each reading time", slope_by_time(cal_pts(), x = input$xaxis, T_C = input$temp, z = input$zA, flag_mV = input$rt_flag)))
  output$rt_suggest = renderTable({
    s = req(suggest_read_time(rt_settle(), input$rt_share))
    data.frame(`Suggested reading time (s)` = ifelse(is.na(s$read_s), "none: not enough standards settled", as.character(s$read_s)),
               `Share settled by then` = sprintf("%.0f%%", 100 * s$share_settled), `Fitted standards` = s$fitted_standards,
               `Never settled` = s$never_settled, check.names = FALSE)
  })
  rt_plot_obj = reactive({
    r = req(rt_slopes())
    S = nernst_S(input$temp, input$zA)
    ggplot2::ggplot(r$long, ggplot2::aes(read_s, slope_mV_dec, colour = day)) +
      ggplot2::geom_hline(yintercept = c(S, -S), linetype = "dashed", colour = "grey60") +
      ggplot2::geom_vline(xintercept = input$read_s, linetype = "dotted") +
      ggplot2::geom_line(linewidth = 0.6) + ggplot2::geom_point(size = 2.4) +
      ggplot2::scale_colour_manual(values = day_colours(sort(unique(r$long$day))), name = "Day") +
      ggplot2::labs(x = "Reading time (s)", y = "Slope (mV/decade)", subtitle = "Dashed: Nernstian slope. Dotted: the reading time in use.") +
      ggplot2::theme_classic(base_size = 12) + watermark()
  })
  output$rt_plot = renderPlot(rt_plot_obj())
  output$rt_spread = renderTable(round_df(req(rt_slopes())$spread))
  output$rt_settle = renderTable(round_df(req(rt_settle())))
  output$rt_csv = downloadHandler("reading_time.csv", function(f) {
    st = rt_settle(); sl = rt_slopes()$long
    write.csv(merge(st, sl, by = "day", all = TRUE), f, row.names = FALSE)
  })

  # Precision
  # readings at the check concentration whose letter is not listed are left out of the
  # repeatability estimate only; the calibration itself is unchanged
  prec_pts = reactive({
    p = cal_pts(); L = letters_from(if (is.null(input$chk_letters)) "" else input$chk_letters)
    if (length(L)) p = p[!(abs(p$conc_M - input$chk_conc) < 1e-12 & !is.na(p$letter) & !p$letter %in% L), ]
    p
  })
  prec_chk = reactive(precision_check(prec_pts(), input$chk_conc, input$temp, input$zA))
  prec_cal = reactive(precision_calibration(cal_pts(), input$xaxis))
  output$prec_notes = renderUI({
    n = c(prec_chk()$notes, prec_cal()$notes)
    if (length(n)) div(class = "alert alert-warning", HTML(paste(htmltools::htmlEscape(n), collapse = "<br>")))
  })
  output$prec_metrics = renderTable({
    m = req(prec_chk()$metrics)
    data.frame(Metric = m$metric, `SD (mV)` = sprintf("%.2f", m$SD_mV),
               `Equivalent concentration error (%)` = sprintf("%.1f", m$conc_error_pct), check.names = FALSE)
  })
  prec_blup_obj = reactive({
    b = req(prec_chk()$blup)
    ggplot2::ggplot(b, ggplot2::aes(day, effect_mV)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.15) +
      ggplot2::geom_point(size = 3, colour = "#1A6FDF") +
      ggplot2::labs(x = "Day", y = "Day effect (mV)", subtitle = "Estimated offset of each day from the overall mean (BLUP, 95% interval)") +
      ggplot2::theme_classic(base_size = 12) + watermark()
  })
  output$prec_blup = renderPlot(prec_blup_obj())
  output$prec_cal = renderTable({
    t = req(prec_cal()$table)
    data.frame(Days = t$days, `Pooled slope (mV/decade)` = sprintf("%.2f ± %.2f (SE)", t$pooled_slope, t$pooled_slope_se),
               `Between-day SD of baseline (mV)` = sprintf("%.2f", t$between_day_baseline_sd),
               `Between-day SD of slope (mV/decade)` = ifelse(is.na(t$between_day_slope_sd), "not estimable", sprintf("%.2f", t$between_day_slope_sd)),
               check.names = FALSE)
  })
  output$prec_csv = downloadHandler("precision.csv", function(f) {
    m = prec_chk()$metrics; m$notes = paste(c(prec_chk()$notes, prec_cal()$notes), collapse = " | ")
    write.csv(m, f, row.names = FALSE)
  })

  # Fixed interference method
  fim_std = reactiveVal(NULL)
  observeEvent(list(input$fim_files, input$scheme), fim_std(standards_table(req(input$fim_files), scheme = if (isTRUE(input$scheme %in% c("july", "protocol"))) input$scheme else detect_scheme(input$fim_files$name))))
  output$fim_std = render_std(reactive(req(fim_std())))
  observeEvent(input$fim_std_cell_edit, apply_edit(fim_std, input$fim_std_cell_edit))
  fim_pts = reactive({
    st = req(fim_std()); st$fit = st$include
    withProgress(message = "Reading EMF16 files", add_emf(st, input$read_s, background_I = input$fim_cB, z = input$zA))
  })
  fim_res = reactive(fim(fim_pts(), letters_from(input$fim_flat), letters_from(input$fim_nernst),
                         input$fim_cB, input$zA, input$fim_zB, input$temp))
  output$fim_summary = renderTable({
    s = req(fim_res())$summary
    data.frame(Days = s$days, K = sprintf("%.3g ± %.2g", s$K_mean, s$K_sd),
               `log K` = sprintf("%.2f ± %.2f", s$lgK_mean, s$lgK_sd), check.names = FALSE)
  })
  output$fim_perday = renderTable(round_df(req(fim_res())$per_day))
  fim_plot_obj = reactive({
    p = fim_pts(); r = req(fim_res())
    lines = do.call(rbind, lapply(seq_len(nrow(r$per_day)), function(i) {
      d = p[p$day == r$per_day$day[i] & p$letter %in% letters_from(input$fim_nernst), ]
      m = lm(E_mV ~ neg_log_a, d)
      xn = c(min(p$neg_log_a) - 0.1, max(d$neg_log_a) + 2)
      rbind(data.frame(day = r$per_day$day[i], part = "flat", x = c(min(d$neg_log_a) + 0.5, max(p$neg_log_a) + 0.2), E_mV = r$per_day$flat_mV[i]),
            data.frame(day = r$per_day$day[i], part = "line", x = xn, E_mV = coef(m)[1] + coef(m)[2] * xn))
    }))
    ggplot2::ggplot(p, ggplot2::aes(neg_log_a, E_mV, colour = day)) +
      ggplot2::geom_line(data = lines, ggplot2::aes(x, E_mV, group = interaction(day, part)), linetype = "dashed") +
      ggplot2::geom_point(size = 2.6) +
      ggplot2::scale_colour_manual(values = day_colours(sort(unique(p$day))), name = "Day") +
      ggplot2::scale_x_continuous(breaks = 1:10) +
      ggplot2::labs(x = expression("−log "*italic(a)*" (primary ion)"), y = "EMF (mV)",
                    subtitle = sprintf("Interferent fixed at %g M, potential read at %g s", input$fim_cB, input$read_s)) +
      ggplot2::theme_classic(base_size = 12) + watermark()
  })
  output$fim_plot = renderPlot(fim_plot_obj())
  output$fim_csv = downloadHandler("fim_results.csv", function(f) write.csv(req(fim_res())$per_day, f, row.names = FALSE))

  # Separate solution method
  ssm_std = reactiveVal(NULL)
  observeEvent(list(input$ssm_files, input$scheme), ssm_std(standards_table(req(input$ssm_files), scheme = if (isTRUE(input$scheme %in% c("july", "protocol"))) input$scheme else detect_scheme(input$ssm_files$name))))
  output$ssm_std = render_std(reactive(req(ssm_std())))
  observeEvent(input$ssm_std_cell_edit, apply_edit(ssm_std, input$ssm_std_cell_edit))
  ssm_res = reactive({
    b = req(ssm_std()); b$fit = b$include
    ib = withProgress(message = "Reading EMF16 files", add_emf(b, input$read_s, input$bgI, input$ssm_zB))
    concs = as.numeric(strsplit(input$ssm_concs, "[,; ]+")[[1]])
    ssm(cal_pts(), ib, concs[!is.na(concs)], input$zA, input$ssm_zB, input$temp)
  })
  output$ssm_table = renderTable(round_df(req(ssm_res())))
  output$ssm_csv = downloadHandler("ssm_results.csv", function(f) write.csv(req(ssm_res()), f, row.names = FALSE))

  # ECSA tab: gold is measured here; other metals and carbon are pointed to the tab with the right method
  ecsa_metal = list(pt = c("pt_h", "Platinum: hydrogen desorption (210 µC/cm²) is the usual method. Oxide reduction (420) or CO stripping (420) also work."),
                    pd = c("pd_ox424", "Palladium: PdO reduction, 424 µC/cm² (hydrogen is absorbed into Pd, so the hydrogen region is not used)."),
                    rh = c("rh_h", "Rhodium: hydrogen desorption (221 µC/cm²) or oxide reduction (442)."),
                    ir = c("ir_h", "Iridium: hydrogen desorption (218 µC/cm²)."))
  output$ecsa_metal_help = renderUI(helpText(ecsa_metal[[req(input$ecsa_material)]][2]))
  # choosing a metal sets the window method below to that metal's usual reference
  observeEvent(input$ecsa_material, if (input$ecsa_material %in% names(ecsa_metal))
    updateSelectInput(session, "sa_preset", selected = ecsa_metal[[input$ecsa_material]][1]))
  observeEvent(input$ecsa_go_cv, bslib::nav_select("cv_sub", "Peaks and kinetics"))
  observeEvent(input$ecsa_go_ca, bslib::nav_select("main_nav", "Chronoamperometry"))
  observeEvent(input$ecsa_go_dl, bslib::nav_select("cv_sub", "Double-layer capacitance"))

  # QC charts: the history is kept in this session and downloaded by the user
  qc_hist = reactiveVal(qc_template())
  observeEvent(input$qc_file, qc_hist(tryCatch(qc_read(input$qc_file$datapath), error = function(e) {
    notify_error(e); qc_hist() })))
  as_qc_date = function(day) { d = as.Date(day, tryFormats = c("%m-%d-%Y", "%Y-%m-%d", "%d/%m/%Y"), optional = TRUE)
    ifelse(is.na(d), format(Sys.Date()), format(d)) }
  observeEvent(input$qc_add_session, {
    rows = list()
    r = safely(cal_res)
    if (!is.null(r)) rows = c(rows, list(
      data.frame(date = as_qc_date(r$per_day$day), quantity = "ISE slope (mV/decade)", value = r$per_day$slope_mV_dec, unit = "mV/decade", electrode = "", note = "ISE calibration tab"),
      data.frame(date = as_qc_date(r$per_day$day), quantity = "ISE intercept (mV)", value = r$per_day$intercept_mV, unit = "mV", electrode = "", note = "ISE calibration tab")))
    b = safely(cal_bg)
    if (!is.null(b)) rows = c(rows, list(data.frame(date = as_qc_date(b$per_day$day), quantity = "ISE background b (M)", value = b$per_day$b_bkg_M, unit = "M", electrode = "", note = "background-corrected fit")))
    e = safely(ecsa_res)
    if (!is.null(e)) rows = c(rows, list(data.frame(date = format(Sys.Date()), quantity = "ECSA (cm²)", value = e[["ECSA (cm²)"]], unit = "cm²", electrode = "", note = "ECSA tab")))
    if (!length(rows)) { showNotification("No results loaded yet: load files on the ISE calibration or ECSA tab first.", type = "warning"); return() }
    qc_hist(rbind(qc_hist(), do.call(rbind, rows)))
  })
  observeEvent(input$qc_add_one, {
    req(is.finite(input$qc_value), nzchar(input$qc_quantity))
    qc_hist(rbind(qc_hist(), data.frame(date = format(input$qc_date), quantity = input$qc_quantity, value = input$qc_value, unit = "",
                                        electrode = input$qc_electrode, note = "entered by hand")))
  })
  qc_nref = reactive(if (is.null(input$qc_nref) || is.na(input$qc_nref)) NA else input$qc_nref)
  qc_tab = reactive({ h = qc_hist(); validate(need(nrow(h) > 0, "Load a history, add this session's results, or add a value.")); qc_table(h, qc_nref()) })
  output$qc_status = renderUI({ t = qc_tab(); k = sum(grepl("1-3s|2-2s|R-4s|4-1s|10x", t$rule))
    p(sprintf("%d sessions, %d quantities, %d with a broken rule, %d warnings.", nrow(t), length(unique(t$quantity)), k, sum(grepl("warning", t$rule)))) })
  qc_plot_obj = reactive({ qc_tab(); qc_plot(qc_hist(), qc_nref()) + g_theme() })
  output$qc_plot = renderPlot(qc_plot_obj())
  output$qc_limits = renderTable(round_df(qc_limits(qc_hist(), qc_nref())))
  output$qc_table = renderTable(round_df(qc_tab()[, c("date", "quantity", "value", "z", "rule", "check", "electrode", "note")]))
  output$qc_template = downloadHandler("electrofast_qc_history_template.csv", function(f) write.csv(qc_template(), f, row.names = FALSE))
  output$qc_download = downloadHandler(function() sprintf("electrofast_qc_history_%s.csv", format(Sys.Date())),
                                       function(f) write.csv(qc_hist(), f, row.names = FALSE))

  # Cell diagram: drawn only from the entries on the Cell diagram tab
  cd_spec = reactive(list(layout = input$cd_layout, mode = input$cd_mode, we_kind = input$cd_we_kind, we = input$cd_we, re_kind = input$cd_re_kind, re = input$cd_re,
    ce_kind = input$cd_ce_kind, ce = input$cd_ce, electrolyte = input$cd_electrolyte, analyte = input$cd_analyte, analyte_ch = input$cd_analyte_ch, a_kind = input$cd_a_kind, a = input$cd_a,
    b_kind = input$cd_b_kind, b = input$cd_b, sol_a = input$cd_sol_a, sol_b = input$cd_sol_b, separator = input$cd_sep,
    rotation = input$cd_rot, purge = input$cd_purge, gas = input$cd_gas, stir = input$cd_stir, instrument = input$cd_instr,
    title = input$cd_title, subscripts = input$cd_sub, labels = input$cd_labels, style = input$cd_style, particles = input$cd_particles))
  cd_out = reactive(cell_diagram(cd_spec()))
  output$cd_plot = renderPlot(cd_out()$plot)
  cd_save = function(dev) function(f) fig_save(f, cd_out()$plot, width = 10, height = 5.6, dpi = 300, device = dev)
  output$cd_png = downloadHandler("cell_diagram.png", cd_save("png"))
  output$cd_pdf = downloadHandler("cell_diagram.pdf", cd_save(echem_pdf))
  # svglite named here, not "svg": ggsave needs it for SVG, and the deploy and desktop builds find packages by name
  output$cd_svg = downloadHandler("cell_diagram.svg", cd_save(svglite::svglite))
  output$cd_has_photo = reactive(!is.null(input$cd_photo)); outputOptions(output, "cd_has_photo", suspendWhenHidden = FALSE)
  output$cd_photo_img = renderImage({ p = req(input$cd_photo)
    list(src = p$datapath, contentType = p$type, style = "max-width:100%; height:auto;", alt = "Photo of the cell") }, deleteFile = FALSE)

  # Report: each part is collected only if its tab has data. A tab without data stops at req(),
  # which is caught here so that part is simply left out.
  safely = function(f) tryCatch(f(), error = function(e) NULL)
  report_parts = function() list(
    # the ISE settings are listed only when ISE files were loaded, so a report on other data does not show them
    settings = c(if (!is.null(input$cal_files) || !is.null(input$cal_sheet)) list(`Reading time (s)` = input$read_s, `Temperature (°C)` = input$temp, `Primary ion charge` = input$zA,
                    `Background ionic strength (M)` = input$bgI, `x-axis` = input$xaxis, `Concentrations from` = if (!is.null(input$cal_sheet)) "Sample sheet" else unname(scheme_labels[input$scheme %||% "names"]),
                    `Standards in the fit` = if (isTRUE(input$scheme %in% c("july", "protocol"))) input$fit_letters else sprintf("%g to %g M", input$fit_cfrom %||% 1e-5, input$fit_cto %||% 0.1),
                    `Second range` = if (isTRUE(input$scheme %in% c("july", "protocol"))) input$cmp_letters else ""),
                 as.list(provenance_lines(input_files(reactiveValuesToList(input))))),
    ecsa = safely(function() list(table = ecsa_res(), plot = ecsa_plot())),
    cal = safely(function() list(std = cal_std(), res = req(cal_res()), cmp = cal_cmp(), bg = cal_bg(), plot = cal_plot(), dx = safely(function() dx()$causes),
                                 steps = steps_table(cal_pts(), input$xaxis, input$temp, input$zA))),
    rt = safely(function() list(slopes = req(rt_slopes()), suggest = suggest_read_time(rt_settle(), input$rt_share), plot = rt_plot_obj())),
    prec = safely(function() list(chk = prec_chk(), cal = prec_cal(), notes = c(prec_chk()$notes, prec_cal()$notes),
                                  plot = safely(prec_blup_obj))),
    fim = safely(function() list(res = req(fim_res()), plot = fim_plot_obj(), cB = input$fim_cB)),
    ssm = safely(function() req(ssm_res())),
    # the cell diagram leads the technique sections when the user ticked "Add the diagram to the report"
    techniques = c(if (isTRUE(input$cd_report)) list(list(title = "Cell", text = cd_out()$key$Part, plots = list(`Cell diagram` = cd_out()$plot), tables = list())),
                   technique_parts(), performance_parts(), Filter(Negate(is.null), list(safely(function() { t = qc_tab()
      list(title = "QC charts", text = sprintf("%d sessions, %d with a broken Westgard rule.", nrow(t), sum(grepl("1-3s|2-2s|R-4s|4-1s|10x", t$rule))),
           plots = list(`Control charts` = qc_plot_obj()),
           tables = list(Limits = round_df(qc_limits(qc_hist(), qc_nref())), Sessions = round_df(t[, c("date", "quantity", "value", "z", "rule", "check")]))) })))))
  # The same parts as slide sections (ISE, ECSA, reading time, precision, selectivity), followed by the technique tabs
  pptx_sections = function(P) Filter(Negate(is.null), c(list(
    if (!is.null(P$ecsa)) list(title = "ECSA (gold)", text = sprintf("ECSA %.4f cm², charge %.2f µC, roughness factor %.2f.", P$ecsa$table[["ECSA (cm²)"]], P$ecsa$table[["Q (µC)"]], P$ecsa$table[["Roughness factor"]]),
                               plots = list(Voltammogram = P$ecsa$plot), tables = list(Results = round_df(P$ecsa$table))),
    if (!is.null(P$cal)) { s = P$cal$res$summary
      list(title = "ISE calibration", text = c(sprintf("Slope %s mV/decade over %s, %.1f%% of Nernstian, %d day(s). R² %.4f to %.4f.", if (is.na(s$slope_sd)) sprintf("%.2f", s$slope_mean) else sprintf("%.2f ± %.2f", s$slope_mean, s$slope_sd), s$range, s$pct_nernst, s$days, s$r2_min, s$r2_max),
                                               if (!is.null(P$cal$bg)) sprintf("Background-corrected fit: %.2f mV/decade, background %.2g M.", P$cal$bg$summary$slope_mean, P$cal$bg$summary$b_bkg_mean)),
           plots = list(Calibration = P$cal$plot), tables = Filter(Negate(is.null), list(`Each day` = round_df(P$cal$res$per_day[, c("day", "n", "slope_mV_dec", "intercept_mV", "r2", "pct_nernst")]), Steps = round_df(P$cal$steps), `Possible causes` = P$cal$dx))) },
    if (!is.null(P$rt)) list(title = "Reading time", text = NULL, plots = list(`Slope by reading time` = P$rt$plot), tables = list(Spread = round_df(P$rt$slopes$spread))),
    if (!is.null(P$prec)) list(title = "Precision", text = P$prec$notes, plots = Filter(Negate(is.null), list(`Day effects` = P$prec$plot)), tables = Filter(Negate(is.null), list(`Check standard` = P$prec$chk$metrics, `Calibration across days` = P$prec$cal$table))),
    if (!is.null(P$fim)) list(title = "Selectivity: fixed interference", text = sprintf("log K = %.2f (interferent %g M).", P$fim$res$summary$lgK_mean, P$fim$cB), plots = list(FIM = P$fim$plot), tables = list(`Each day` = round_df(P$fim$res$per_day))),
    if (!is.null(P$ssm)) list(title = "Selectivity: separate solution", text = NULL, plots = list(), tables = list(Results = round_df(P$ssm)))),
    P$techniques))
  output$report_pptx = downloadHandler(
    filename = function() sprintf("ElectroFast_slides_%s.pptx", format(Sys.Date())),
    content = function(file) withProgress(message = "Building the slides", {
      # the title slide lists the cell only from what was typed on the Cell diagram tab
      det = c(Electrolyte = input$cd_electrolyte %||% "", `Working electrode` = input$cd_we %||% "", `Reference electrode` = input$cd_re %||% "",
              `Counter electrode` = input$cd_ce %||% "")
      P = report_parts(); build_pptx(file, input$rep_title, input$rep_author %||% "", P$settings, pptx_sections(P), det) }))
  output$report_pdf = downloadHandler(
    filename = function() sprintf("ElectroFast_report_%s.pdf", format(Sys.Date())),
    content = function(file) withProgress(message = "Building the PDF", {
      P = report_parts(); build_pdf(file, input$rep_title, input$rep_author %||% "", P$settings, pptx_sections(P)) }))
  output$report_bundle = downloadHandler(
    filename = function() sprintf("ElectroFast_bundle_%s.zip", format(Sys.time(), "%Y-%m-%d_%H%M")),
    content = function(file) withProgress(message = "Building the bundle", {
      P = report_parts(); build_bundle(file, reactiveValuesToList(input), pptx_sections(P)) }))
  output$report_html = downloadHandler(
    filename = function() sprintf("ElectroFast_report_%s.html", format(Sys.Date())),
    content = function(file) {
      # Pandoc is missing in the browser version (Shinylive): build the same report directly in R there
      if (!rmarkdown::pandoc_available()) return(withProgress(message = "Building the report", {
        P = report_parts(); build_html(file, input$rep_title, input$rep_author %||% "", P$settings, pptx_sections(P)) }))
      tmp = file.path(tempdir(), "electrofast_report.Rmd")
      file.copy("report/electrofast_report.Rmd", tmp, overwrite = TRUE)
      withProgress(message = "Building the report", {
        parts = report_parts()
        rmarkdown::render(tmp, output_file = file, quiet = TRUE, envir = new.env(parent = globalenv()),
                          params = list(parts = parts, title = input$rep_title, author = input$rep_author))
      })
    })
}

shinyApp(ui, server)
