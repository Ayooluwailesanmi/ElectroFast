# Look of the app: one accent colour (the blue of the plots), a near-white page, one sans-serif face, small labels,
# square-ish corners and no shadows, digits that line up in tables, and a hint in empty plot areas.

electrofast_theme = function() {
  # local = TRUE downloads the font files in R, which R cannot do inside a browser (Shinylive, os "emscripten"); there the page loads them
  font = tryCatch(bslib::font_google("Source Sans 3", local = !identical(R.version$os, "emscripten")), error = function(e) NULL)
  bs_theme(version = 5,
    bg = "#FBFBFA", fg = "#1F2328", primary = "#1A6FDF", secondary = "#5B6470",
    base_font = if (is.null(font)) bslib::font_collection("Segoe UI", "Helvetica Neue", "Arial", "sans-serif") else
      bslib::font_collection(font, "Segoe UI", "Arial", "sans-serif"),
    font_scale = 0.92, "border-radius" = "3px", "border-radius-lg" = "4px", "border-radius-sm" = "2px",
    "navbar-light-bg" = "#FFFFFF", "bslib-sidebar-bg" = "#F3F4F6") |>
    bs_add_rules("
      .navbar { border-bottom: 1px solid #E3E5E8; }
      .navbar-brand { font-weight: 700; letter-spacing: 0.01em; }
      .card, .bslib-card { box-shadow: none !important; border-color: #E3E5E8; }
      .bslib-sidebar-layout { --_sidebar-bg: #F3F4F6; --bslib-sidebar-bg: #F3F4F6; --_sidebar-border: 1px solid #E3E5E8; }
      .bslib-sidebar-layout > .sidebar, .bslib-sidebar-layout > aside.sidebar { background-color: #F3F4F6 !important; }
      .bslib-sidebar-layout > .sidebar > .sidebar-content { gap: 0.35rem !important; }
      .bslib-sidebar-layout > .sidebar .form-group { margin-bottom: 0.35rem !important; }
      .sidebar .control-label, .sidebar label { font-size: 0.82rem; font-weight: 600; color: #3A4048; margin-bottom: 0.2rem; }
      .sidebar .help-block, .sidebar .form-text { font-size: 0.78rem; color: #6A717B; line-height: 1.35; margin-top: 0.1rem; }
      .sidebar .shiny-input-container { margin-bottom: 0.55rem; }
      .sidebar .shiny-input-container .progress { margin-top: 0.3rem; margin-bottom: 0.2rem; height: 1.05rem; font-size: 0.72rem; }
      .sidebar .form-control, .sidebar .selectize-input { font-size: 0.88rem; padding-top: 0.3rem; padding-bottom: 0.3rem; min-height: 0; }
      .sidebar .btn { font-size: 0.85rem; }
      .sidebar .btn-default, .sidebar .shiny-download-link { padding: 0.25rem 0.6rem; border-color: #C9CED4; background: #FFFFFF; color: #2B3138; }
      .card-header { background: #FFFFFF; }
      .card-header:not(:has(.nav)) { font-size: 0.8rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.04em; color: #5B6470; }
      .nav-tabs .nav-link, .navbar .nav-link { font-size: 0.92rem; }
      .card-body > h6 { font-size: 0.8rem; text-transform: uppercase; letter-spacing: 0.04em; color: #5B6470; margin-top: 0.8rem; }
      table td, table th { font-variant-numeric: tabular-nums; }
      table.table { font-size: 0.86rem; }
      h5, h6 { font-weight: 600; }
      .shiny-plot-output:not(:has(img)) { background: #F7F8F9; border: 1px dashed #D5D9DE; border-radius: 3px; position: relative; }
      .shiny-plot-output:not(:has(img))::before { content: 'Load files in the sidebar to see the plot here.';
        position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); color: #8A919A; font-size: 0.9rem; }
    ") |>
    bs_add_rules(technique_icons())
}

# Small line drawings of the curve each technique gives, shown before the tab names (CSS only, matched on the tab
# names, so no tab code changes). Same blue as the plots.
technique_icons = function(col = "#1A6FDF") {
  path = c(
    "Cyclic voltammetry" = "M2 11 C5 11 6 4 8 6 C9 7 10 8 14 8 M14 8 C11 9 10 12 8 10 C6 8 5 12 2 11",
    "Peaks and kinetics" = "M2 11 C5 11 6 4 8 6 C9 7 10 8 14 8 M14 8 C11 9 10 12 8 10 C6 8 5 12 2 11",
    "Linear sweep voltammetry" = "M2 13 C8 13 9 9 14 3",
    "Pulse voltammetry (DPV, SWV)" = "M2 12 L5 12 C7 12 7 4 8 4 C9 4 9 12 11 12 L14 12",
    "Potentiometry" = "M2 13 L14 3",
    "Equations" = "M3 4 L13 4 M3 8 L13 8 M3 12 L9 12",
    "CV: peaks, scan rate, kinetics" ="M2 11 C5 11 6 4 8 6 C9 7 10 8 14 8 M14 8 C11 9 10 12 8 10 C6 8 5 12 2 11",
    "Surface area (CV peak)" = "M2 12 L5 12 C7 12 7 4 8 4 C9 4 9 12 11 12 L14 12",
    "ECSA" = "M2 12 L5 12 C7 12 7 4 8 4 C9 4 9 12 11 12 L14 12",
    "Double-layer capacitance" = "M3 5 L13 5 L13 11 L3 11 Z",
    "Tafel" = "M2 13 L6 9 L14 3",
    "DPV / SWV calibration" = "M2 12 L5 12 C7 12 7 4 8 4 C9 4 9 12 11 12 L14 12",
    "Open-circuit potential" = "M2 5 C5 9 8 8 14 8",
    "Chronoamperometry" = "M2 14 L4 14 L4 3 C5 9 8 11 14 12",
    "Chronocoulometry" = "M2 14 L4 14 L4 10 C6 6 9 5 14 4",
    "Chronopotentiometry" = "M2 11 C6 10 8 9 9 7 C10 5 10 3 14 3",
    "Rotating disk (Levich, Koutecky-Levich)" = "M2 12 C7 12 7 4 14 4",
    "RRDE" = "M2 12 C7 12 7 4 14 4 M2 14 C7 14 7 9 14 9",
    "Impedance (EIS)" = "M2 13 A4 4 0 0 1 10 13 L14 8",
    "ISE calibration" = "M2 13 L14 3",
    "ISE" = "M2 13 L14 3",
    "Selectivity" = "M2 6 L6 6 L14 13",
    "Precision" = "M2 8 L14 8 M5 6 L5 10 M8 5 L8 11 M11 6 L11 10",
    # Sensor performance (Selectivity and Precision share the drawings above with the Potentiometry tabs)
    "Sensitivity" = "M2 2 L2 14 L14 14 M4 12 L13 4",
    "Specificity" = "M4 14 L4 3 M8 14 L8 11 M12 14 L12 12",
    "Stability" = "M2 4 L14 6 M2 10 L4 10 M6 10 L8 10 M10 10 L12 10",
    "Speed" = "M2 13 L5 13 C6 5 8 4 14 4",
    # top bar
    "Cell diagram" = "M4 3 L4 13 L12 13 L12 3 M4 8 L12 8 M7 3 L7 10",
    "Report" = "M4 2 L10 2 L13 5 L13 14 L4 14 Z M10 2 L10 5 L13 5 M6 8 L11 8 M6 11 L11 11",
    "Guide" = "M2 4 L7 3 L8 4 L9 3 L14 4 L14 13 L9 12 L8 13 L7 12 L2 13 Z M8 4 L8 13",
    "About" = "M8 1.5 A6.5 6.5 0 1 1 7.99 1.5 M8 7 L8 11.5 M8 4.6 L8 4.9",
    "QC charts" = "M2 3 L14 3 M2 13 L14 13 M2 8 L5 6 L8 10 L11 7 L14 9")
  svg = function(d) paste0("url(\"data:image/svg+xml,", utils::URLencode(sprintf(
    "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 16 16' width='16' height='16'><path d='%s' fill='none' stroke='%s' stroke-width='1.5' stroke-linecap='round' stroke-linejoin='round'/></svg>", d, col), reserved = TRUE), "\")")
  paste(sprintf(".nav-link[data-value=\"%1$s\"]::before, .dropdown-item[data-value=\"%1$s\"]::before, .dropdown-toggle[data-value=\"%1$s\"]::before { content: %2$s; display: inline-block; width: 16px; height: 16px; margin-right: 6px; vertical-align: -3px; }",
                names(path), vapply(path, svg, "")), collapse = "\n")
}

# The name in the top bar, with a small cyclic voltammogram before it
electrofast_brand = function() tags$span(
  HTML("<svg viewBox='0 0 16 16' width='22' height='22' style='vertical-align:-5px;margin-right:6px'><path d='M2 11 C5 11 6 4 8 6 C9 7 10 8 14 8 M14 8 C11 9 10 12 8 10 C6 8 5 12 2 11' fill='none' stroke='#1A6FDF' stroke-width='1.6' stroke-linecap='round'/></svg>"),
  "ElectroFast")
