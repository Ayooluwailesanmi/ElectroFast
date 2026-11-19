# Electrode panel shared by the tabs: the user picks a standard disk, types a disk or
# microelectrode diameter, gives the dimensions of a plate, strip, wire, band or
# (hemi)sphere, or enters the area directly. electrode_server() returns the geometric
# area (cm2), the radius (cm, for microelectrode equations) and a short description.
# (base R has %||% from 4.4.0; this keeps older R versions working)
if (!exists("%||%")) `%||%` = function(a, b) if (is.null(a)) b else a

electrode_ui = function(id, label = "Electrode", allow_none = FALSE, default = "preset") {
  ns = NS(id)
  shapes = c("Standard disk (pick a size)" = "preset", "Disk, enter the diameter" = "disk",
             "Microdisk, enter the diameter in µm" = "microdisk", "Plate or strip" = "rectangle",
             "Wire or rod" = "wire", "Band" = "band", "Hemisphere" = "hemisphere", "Sphere" = "sphere",
             "Enter the area (cm²)" = "area")
  if (allow_none) shapes = c("Not given" = "none", shapes)
  sizes = stats::setNames(as.character(electrode_sizes$d_mm), electrode_sizes$label)
  tagList(
    selectInput(ns("shape"), label, shapes, selected = default),
    conditionalPanel("input.shape == 'preset'", ns = ns, selectInput(ns("preset"), "Size", sizes, selected = "2")),
    conditionalPanel("['disk', 'hemisphere', 'sphere', 'wire'].indexOf(input.shape) >= 0", ns = ns,
                     numericInput(ns("d_mm"), "Diameter (mm)", 2, min = 0)),
    conditionalPanel("input.shape == 'microdisk'", ns = ns, numericInput(ns("d_um"), "Diameter (µm)", 10, min = 0)),
    conditionalPanel("['rectangle', 'wire', 'band'].indexOf(input.shape) >= 0", ns = ns,
                     numericInput(ns("len_mm"), "Length in the solution (mm)", 10, min = 0)),
    conditionalPanel("['rectangle', 'band'].indexOf(input.shape) >= 0", ns = ns,
                     numericInput(ns("wid_mm"), "Width (mm)", 5, min = 0)),
    conditionalPanel("input.shape == 'rectangle'", ns = ns,
                     radioButtons(ns("faces"), NULL, c("One face exposed" = 1, "Both faces" = 2), inline = TRUE),
                     checkboxInput(ns("edges"), "Add the immersed edges", FALSE),
                     conditionalPanel("input.edges", ns = ns, numericInput(ns("thick_mm"), "Thickness (mm)", 0.1, min = 0))),
    conditionalPanel("input.shape == 'wire'", ns = ns, checkboxInput(ns("wire_end"), "Add the flat end", TRUE)),
    conditionalPanel("input.shape == 'area'", ns = ns, numericInput(ns("area"), "Area (cm²)", 0.0314, min = 0)),
    textOutput(ns("area_txt"), container = function(...) tags$small(class = "text-muted", ...)))
}

electrode_server = function(id) moduleServer(id, function(input, output, session) {
  geo = reactive({
    s = input$shape; if (is.null(s) || s == "none") return(list(area_cm2 = NA_real_, r_cm = NA_real_, shape = "none", text = "No electrode size given"))
    d = switch(s, preset = as.numeric(input$preset), microdisk = input$d_um / 1000, input$d_mm)
    shape = switch(s, preset = "disk", microdisk = "disk", s)
    a = tryCatch(electrode_area(shape, d_mm = d, length_mm = input$len_mm, width_mm = input$wid_mm,
                                faces = as.numeric(input$faces %||% 1), edges = isTRUE(input$edges),
                                thickness_mm = input$thick_mm %||% 0, wire_end = isTRUE(input$wire_end), area_cm2 = input$area),
                 error = function(e) NA_real_)
    r = if (shape %in% c("disk", "hemisphere", "sphere") && !is.null(d)) d / 20 else NA_real_
    txt = sprintf("Geometric area %s cm²%s", format(signif(a, 4)), if (!is.na(r) && r < 0.0025) sprintf(", radius %g µm", r * 1e4) else "")
    list(area_cm2 = a, r_cm = r, shape = shape, text = txt)
  })
  output$area_txt = renderText(geo()$text)
  geo
})
