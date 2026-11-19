# User interface of the Sensor performance menu: sensitivity, selectivity, specificity, stability, speed and
# precision from tables of results of any technique (calculations in R/performance.R, server side in
# app_performance_server.R). The ISE methods keep their own tabs under Techniques > Potentiometry.

perf_table_input = function(p, what) tagList(
  fileInput(paste0(p, "_file"), paste(what, "(CSV, one row per measurement)"), accept = ".csv"),
  downloadButton(paste0(p, "_template"), "Empty template (CSV)", class = "btn-sm mb-2"))
perf_units = function(p, sig = "µA") layout_columns(textInput(paste0(p, "_sunit"), "Signal unit", sig))
perf_result = function(p, h = "400px", note = NULL) tagList(
  card(plotOutput(paste0(p, "_plot"), height = h, fill = FALSE)),
  card(card_header("Results"), tableOutput(paste0(p, "_table")), uiOutput(paste0(p, "_notes")), downloadButton(paste0(p, "_csv"), "Results (CSV)", class = "btn-sm"),
       if (!is.null(note)) tags$p(class = "small text-muted mt-3 mb-0", note)))

perf_panel_sensitivity = function() nav_panel("Sensitivity",
  layout_sidebar(
    sidebar = sidebar(width = 330,
      perf_table_input("ps", "Calibration results"),
      helpText("Columns: concentration and signal (peak current, charge, potential or any other result). Rows with",
               "concentration 0 are blanks. Two or more blanks give the noise s for the LOD, otherwise the residual SD of the line is used."),
      layout_columns(textInput("ps_cunit", "Concentration unit", "µM"), textInput("ps_sunit", "Signal unit", "µA")),
      layout_columns(numericInput("ps_from", "Linear range from", NA, min = 0), numericInput("ps_to", "to", NA, min = 0)),
      numericInput("ps_area", "Electrode area for the slope per area (cm², optional)", NA, min = 0),
      helpText("Sensitivity is the slope of the calibration line. LOD = 3.3 s/slope and LOQ = 10 s/slope",
               "(ICH Q2(R2), 2023, section 3.2.3.3).")),
    perf_result("ps")))

perf_panel_selectivity = function() nav_panel("Selectivity",
  layout_sidebar(
    sidebar = sidebar(width = 330,
      perf_table_input("pl", "Signals of the target without and with each interferent"),
      helpText("Columns: interferent, signal_analyte (the target alone) and signal_mixture (the same target",
               "concentration with the interferent added), and optionally interferent_concentration. Repeated rows of one",
               "interferent are averaged."),
      numericInput("pl_tol", "Tolerated change (%)", 5, min = 0),
      helpText("The change is 100 × (with − without) / without. Within the tolerance, the interferent does not interfere",
               "at that concentration. Add each interferent at its expected concentration to the target at the middle of its",
               "expected range: the result holds for those concentrations (Thévenot et al., Pure Appl. Chem. 1999, 71, 2333, section 5.2).",
               "The response to each species alone is on the Specificity tab.",
               "For an ISE, use Techniques > Potentiometry > Selectivity (IUPAC methods for the selectivity coefficient).")),
    perf_result("pl")))

perf_panel_specificity = function() nav_panel("Specificity",
  layout_sidebar(
    sidebar = sidebar(width = 330,
      perf_table_input("pc", "Responses to the target and to each other species alone"),
      helpText("Columns: species, signal and target (yes for the target, no for the others), all at the same concentration."),
      perf_units("pc"),
      helpText("Cross-reactivity is the signal of each species as a percentage of the target signal. The selectivity ratio is",
               "the target signal divided by the signal of the species at the same concentration (Thévenot et al., Pure Appl.",
               "Chem. 1999, 71, 2333, section 5.2, which also recommends comparing a full calibration curve of each species",
               "with that of the target).")),
    perf_result("pc", note = "Note: this tab measures cross-reactivity, the response to each other species alone. Together with the Selectivity tab (the target with each interferent added), it describes how selective the sensor is. IUPAC recommends the term selectivity, which can be graded, and discourages specificity for the same concept. A method is either specific or not, and few methods, if any, are specific (Vessman et al., Pure Appl. Chem. 2001, 73, 1381).")))

perf_panel_stability = function() nav_panel("Stability",
  layout_sidebar(
    sidebar = sidebar(width = 330,
      perf_table_input("pt", "Signal over time"),
      helpText("Columns: time and signal, for the same standard measured at intervals (days of storage, hours of",
               "operation, or cycles). Repeated rows at one time are averaged."),
      layout_columns(textInput("pt_tunit", "Time unit", "days"), textInput("pt_sunit", "Signal unit", "µA")),
      numericInput("pt_limit", "Retention limit (%)", 90, min = 0, max = 100),
      helpText("Retention is the signal as a percentage of the first one. The time below the limit is interpolated between",
               "the two measurements on either side.")),
    perf_result("pt")))

perf_panel_speed = function() nav_panel("Speed",
  layout_sidebar(
    sidebar = sidebar(width = 330,
      file_in("pq_file", "One recording with an addition (EMF16, or an export from a potentiostat)", accept = file_types),
      selectInput("pq_channel", "EMF16 channel", choices = NULL),
      numericInput("pq_tadd", "Time of the addition (s)", NA, min = 0),
      layout_columns(numericInput("pq_pre", "Baseline before it (s)", 2, min = 0), numericInput("pq_post", "Final level over the last (s)", 10, min = 0)),
      numericInput("pq_smooth", "Running median over (s)", 0.2, min = 0),
      helpText("The baseline is the median before the addition, the final level the median at the end of the recording.",
               "t90 and t95 are the times after the addition at which the smoothed signal has covered 90% and 95% of the change.",
               "Record long enough after the addition for the signal to settle.")),
    perf_result("pq")))

perf_panel_precision = function() nav_panel("Precision",
  layout_sidebar(
    sidebar = sidebar(width = 330,
      perf_table_input("pp", "Repeated results"),
      helpText("Columns: value (one result per row, for example the signal of one standard), and any of day, electrode",
               "and laboratory. A column with only one level is ignored."),
      perf_units("pp"),
      helpText("The terms of ICH Q2(R2) are as thus. Repeatability is the precision under the same conditions over a short time.",
               "Intermediate precision covers days, analysts and equipment within one laboratory. Reproducibility is the",
               "precision between laboratories. Electrode-to-electrode variation, often called reproducibility in sensor",
               "papers, is part of the intermediate precision. Each component is a variance from a mixed model (REML, lme4).",
               "Five or more levels of each grouping are usually needed. For an ISE, the check-standard analysis is under",
               "Techniques > Potentiometry > Precision.")),
    perf_result("pp")))

performance_ui = function() nav_menu("Sensor performance",
  perf_panel_sensitivity(), perf_panel_selectivity(), perf_panel_specificity(), perf_panel_stability(), perf_panel_speed(), perf_panel_precision())
