# Cell diagram: every layout draws, the labels hold only what was entered (a blank shows the role alone), every text in
# the drawing comes from the entries, the role names or the fixed terminal names, formula subscripts, and the tab,
# downloads and report section through the app. Run from the project folder: Rscript tests/check_cell_diagram.R
`%||%` = function(a, b) if (is.null(a)) b else a
source("R/cell_diagram.R")
chk = function(ok, what) { cat(if (isTRUE(ok)) "ok  " else "FAIL", what, "\n"); if (!isTRUE(ok)) quit(status = 1) }
texts = function(p) unname(unlist(lapply(p$layers, function(l) if (inherits(l$geom, "GeomText")) as.character(l$data$label))))
roles = c("Working electrode", "Reference electrode", "Counter electrode", "Electrolyte", "Analyte", "Gas purge", "Stir bar",
          "Rotating disk electrode", "Electrode in chamber A", "Electrode in chamber B", "Chamber A", "Chamber B", "Separator",
          "Analyte, chamber A", "Analyte, chamber B", "Reference and counter electrode")
# a drawn text is allowed only if it is a role, a terminal or one of the given names, or part of an entered value
# (the "+" between the electrolyte and the analyte is the only added text)
only_entered = function(d, extra = character(0)) all(vapply(setdiff(texts(d$plot), "+"), function(t)
  t %in% c(roles, "RE", "WE", "CE", "Potentiostat", "Potentiometer", extra) || any(grepl(t, d$key$Part, fixed = TRUE)), TRUE))

s1 = list(layout = "single", we = "Au disk, 2 mm", re = "Ag/AgCl", ce = "Pt wire", electrolyte = "0.5 M H2SO4", analyte = "50 mM ethanol",
          purge = TRUE, gas = "N2", stir = TRUE, instrument = "CHI760E", title = "Gold surface area", subscripts = TRUE)
d1 = cell_diagram(s1)
chk(setequal(d1$key$Part, c("Working electrode: Au disk, 2 mm", "Reference electrode: Ag/AgCl", "Counter electrode: Pt wire",
                            "Electrolyte: 0.5 M H₂SO₄", "Analyte: 50 mM ethanol", "Gas purge: N₂", "Stir bar")), "one compartment: labels are exactly the entries")
chk(only_entered(d1, c("CHI760E", "Gold surface area")), "one compartment: no text in the drawing beyond the entries, roles and terminal names")
chk(!any(grepl("^[0-9]+$", texts(d1$plot))), "no numbered markers")
chk(any(texts(d1$plot) == "CHI760E"), "instrument name kept as typed (no subscripts in names)")

blank = cell_diagram(list(layout = "single"))
chk(setequal(blank$key$Part, c("Working electrode", "Reference electrode", "Counter electrode", "Electrolyte")),
    "blank entries: roles only, no guessed materials, no analyte, purge or stir bar unless entered or ticked")
chk(any(texts(blank$plot) == "Potentiostat") && only_entered(blank), "blank instrument: the role name Potentiostat")

s2 = list(layout = "two", mode = "potentiometric", a_kind = "reference", a = "Ag/AgCl", b_kind = "reference", b = "Ag/AgCl",
          sol_a = "1 mM LiClO4", sol_b = "LiClO4 standards", separator = "LASPT sheet", analyte = "Li+ standards", analyte_ch = "B", subscripts = TRUE)
d2 = cell_diagram(s2)
chk(setequal(d2$key$Part, c("Electrode in chamber A: Ag/AgCl", "Electrode in chamber B: Ag/AgCl", "Chamber A: 1 mM LiClO₄",
                            "Chamber B: LiClO₄ standards", "Separator: LASPT sheet", "Analyte, chamber B: Li+ standards")), "two compartments, potentiometric: labels")
chk(only_entered(d2) && any(texts(d2$plot) == "Potentiometer") && !any(texts(d2$plot) %in% c("WE", "RE", "CE")),
    "potentiometric: potentiometer, no potentiostat terminal names, nothing beyond the entries")
chk("Analyte, chamber A: Li+ standards" %in% cell_diagram(modifyList(s2, list(analyte_ch = "A")))$key$Part, "analyte can be placed in chamber A")
d2b = cell_diagram(modifyList(s2, list(mode = "potentiostat", we = "GC disk", re = "SCE", ce = "Pt mesh", ce_kind = "plate")))
chk(all(c("WE", "RE", "CE") %in% texts(d2b$plot)) && length(grep("^(Working|Reference|Counter) electrode", d2b$key$Part)) == 3 && only_entered(d2b),
    "two compartments, potentiostat: three electrodes and terminals")

d3 = cell_diagram(list(layout = "rde", we = "Pt disk, 5 mm", rotation = "400 to 2500 rpm"))
chk("Rotation: 400 to 2500 rpm" %in% d3$key$Part && !any(grepl("Stir", d3$key$Part)) && only_entered(d3), "rotating disk: rotation entry, no stir bar")
chk("Electrolyte: 0.5 M H2SO4" %in% cell_diagram(modifyList(s1, list(subscripts = FALSE)))$key$Part, "subscripts off: text as typed")
chk(identical(cd_subscripts(c("0.1 M Mg(ClO4)2", "5 mM K3[Fe(CN)6]", "1 M KCl", "PBS pH 7.4")),
              c("0.1 M Mg(ClO₄)₂", "5 mM K₃[Fe(CN)₆]", "1 M KCl", "PBS pH 7.4")), "subscripts: brackets, concentrations and pH")
for (k in c("disk", "plate", "wire", "rod", "coil", "reference", "double", "pseudo"))
  invisible(cell_diagram(list(layout = "two", mode = "potentiometric", a_kind = k, b_kind = k)))
cat("ok   every electrode symbol draws\n")
for (k in cell_re_kinds) { d = cell_diagram(list(layout = "single", re_kind = k, re = "Ag/AgCl"))
  chk("Reference electrode: Ag/AgCl" %in% d$key$Part, paste("reference electrode choice", k, "draws, label unchanged")) }
tx = function(m) texts(cell_diagram(modifyList(s1, list(labels = m)))$plot)
chk(!any(tx("short") %in% roles) && all(c("Au disk, 2 mm", "Pt wire") %in% tx("short")), "short labels (default): descriptions only, no part names")
chk(!"Stir bar" %in% tx("short"), "short labels: a part with no description gets no label")
sep_blank = texts(cell_diagram(modifyList(s2, list(separator = "", labels = "short")))$plot)
chk(!"Separator" %in% sep_blank, "short labels: a blank separator gets no label (type its name to label it)")
chk("LASPT sheet" %in% texts(cell_diagram(modifyList(s2, list(labels = "short")))$plot), "short labels: the separator description when entered")
chk(identical(tx("short")[tx("short") %in% c("0.5 M H₂SO₄", "+", "50 mM ethanol")], c("0.5 M H₂SO₄", "+", "50 mM ethanol")),
    "short labels: electrolyte, a centred + on its own line, then the analyte")
only_an = texts(cell_diagram(modifyList(s1, list(electrolyte = "", labels = "short")))$plot)
chk("50 mM ethanol" %in% only_an && !"+" %in% only_an, "short labels: no + when only the analyte is entered")
chk("Analyte: 50 mM ethanol" %in% d1$key$Part, "report entry still names the analyte")
chk(all(c("Working electrode", "Counter electrode", "Stir bar") %in% tx("full")), "full labels: part names and descriptions")
chk(all(tx("none") %in% c("RE", "WE", "CE", "CHI760E", "Gold surface area")), "no labels: only the terminals, instrument and title")
chk(identical(cell_diagram(modifyList(s1, list(labels = "none")))$key$Part, d1$key$Part), "report entries are the same in every label mode")
for (sp in list(s1, s2, modifyList(s2, list(mode = "potentiostat")), list(layout = "rde", we = "Pt disk"))) {
  flat = cell_diagram(modifyList(sp, list(style = "flat"))); d3d = cell_diagram(modifyList(sp, list(style = "3d")))
  chk(identical(sort(texts(d3d$plot)), sort(texts(flat$plot))) && identical(d3d$key$Part, flat$key$Part),
      paste("3D style,", sp$layout, sp$mode %||% "", ": same labels and entries as flat"))
}
for (k in c("disk", "plate", "wire", "rod", "coil", "reference", "double", "pseudo"))
  invisible(cell_diagram(list(layout = "two", mode = "potentiometric", a_kind = k, b_kind = k, style = "3d")))
cat("ok   every electrode symbol draws in 3D style\n")
chk(is.null(s1$style) && !identical(cell_diagram(s1)$plot$layers[[1]]$data, cell_diagram(modifyList(s1, list(style = "3d")))$plot$layers[[1]]$data),
    "flat is the default; 3D style changes the drawing")
pc = list(layout = "pc", we = "6 µm Pt disk", re = "Ag/AgCl", analyte = "Carbon particles", labels = "full")
dpc = cell_diagram(pc)
chk(setequal(dpc$key$Part, c("Working electrode: 6 µm Pt disk", "Reference and counter electrode: Ag/AgCl", "Electrolyte", "Analyte: Carbon particles")),
    "patch-clamp layout: two electrodes, no counter electrode")
chk(all(c("Headstage", "Ground", "Patch-clamp amplifier") %in% texts(dpc$plot)) && !any(texts(dpc$plot) %in% c("WE", "RE", "CE")) &&
    only_entered(dpc, c("Headstage", "Ground", "Patch-clamp amplifier")), "patch-clamp layout: headstage and ground terminals, nothing beyond the entries")
set.seed(11); a1 = runif(3); set.seed(11); invisible(cell_diagram(modifyList(pc, list(particles = TRUE)))); a2 = runif(3)
chk(identical(a1, a2), "particles use a fixed pattern and leave R's random numbers untouched")
npoly = function(d) nrow(ggplot2::ggplot_build(d$plot)$data[[1]])
chk(npoly(cell_diagram(modifyList(pc, list(particles = TRUE)))) > npoly(dpc), "particles drawn only when ticked")
invisible(cell_diagram(modifyList(pc, list(style = "3d", particles = TRUE)))); cat("ok   patch-clamp layout draws in 3D style
")
f = tempfile(fileext = ".svg"); ggplot2::ggsave(f, d1$plot, width = 10, height = 5.6); chk(file.size(f) > 1000, "SVG export")
# PDF: on Windows Arial has no subscript digits (H₂SO₄ printed boxes), so echem_pdf uses Segoe UI there
source("R/report_export.R")
f = tempfile(fileext = ".pdf"); ggplot2::ggsave(f, d1$plot, width = 10, height = 5.6, device = echem_pdf); chk(file.size(f) > 1000, "PDF export")
if (.Platform$OS.type == "windows" && requireNamespace("pdftools", quietly = TRUE))
  chk(all(grepl("Segoe", pdftools::pdf_fonts(f)$name)), "PDF on Windows uses Segoe UI, which has subscript digits")
cat("cell diagram checks passed\n")

# through the app: tab, downloads and the report section
library(shiny)
testServer(shinyAppDir("."), {
  session$setInputs(cd_layout = "single", cd_we_kind = "disk", cd_we = "Au disk", cd_re_kind = "double", cd_re = "Ag/AgCl", cd_ce_kind = "coil", cd_ce = "Pt wire",
                    cd_electrolyte = "0.5 M H2SO4", cd_analyte = "50 mM ethanol", cd_analyte_ch = "B", cd_purge = FALSE, cd_stir = FALSE,
                    cd_instr = "", cd_title = "", cd_sub = TRUE, cd_report = FALSE)
  chk(all(c("Electrolyte: 0.5 M H₂SO₄", "Analyte: 50 mM ethanol") %in% cd_out()$key$Part), "app: diagram from the tab's entries, with the analyte")
  invisible(output$cd_plot)
  for (id in c("cd_png", "cd_pdf", "cd_svg")) { p = output[[id]]; chk(file.exists(p) && file.size(p) > 1000, paste("app: download", id)) }
  chk(!any(vapply(report_parts()$techniques, function(s) identical(s$title, "Cell"), TRUE)), "app: not in the report unless ticked")
  session$setInputs(cd_report = TRUE)
  secs = report_parts()$techniques
  chk(identical(secs[[1]]$title, "Cell") && identical(secs[[1]]$text, cd_out()$key$Part), "app: report section when ticked, first, with the entries")
})
cat("cell diagram app checks passed\n")
