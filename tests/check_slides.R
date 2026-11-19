# Builds the PowerPoint deck through the app's download button with the simulated example data (www/examples) and
# checks that it follows the laboratory's report decks: 16:9, a title slide with the cell typed on the Cell diagram
# tab, a main-points slide, a slide per plot with its result underneath, 13 pt Calibri in the tables, and the
# settings at the end. Run from the project folder: Rscript tests/check_slides.R [copy of the deck to keep]
library(shiny); app = shinyAppDir(".")
keep = commandArgs(TRUE)[1]
ex = tempfile(); utils::unzip("www/examples/ElectroFast_example_data_simulated.zip", exdir = ex)
ex = file.path(ex, "ElectroFast_simulated_data")
up = function(p, type = "text/csv") data.frame(name = basename(p), size = file.size(p), type = type, datapath = normalizePath(p))
cvs = list.files(file.path(ex, "CV_eChem"), "^CV_reversible", full.names = TRUE)
chk = function(ok, what) { cat(if (ok) "ok  " else "FAIL", what, "\n"); if (!ok) quit(status = 1) }

testServer(app, {
  session$setInputs(instrument = "auto", sr_dir = "cathodic", sr_use_window = FALSE, sr_species = "diff", sr_n = 1, sr_C = 1, sr_T = 25, sr_alpha = 0.5,
                    `sr_el-shape` = "area", `sr_el-area` = 0.01, sr_mix = FALSE, sr_files = up(cvs))
  session$setInputs(sr_probe = "other", sr_D = 1e-5, sr_solve = "area", sr_rev = "rev")
  session$setInputs(pp_file = up(file.path(ex, "precision_8days_6electrodes.csv")))
  session$setInputs(qc_file = up(file.path(ex, "qc_history_2_quantities.csv")), qc_nref = 12)
  session$setInputs(cd_electrolyte = "0.1 M KCl (simulated)", cd_we = "Disk, 0.01 cm² (simulated)", cd_re = "", cd_ce = "",
                    rep_title = "ElectroFast slides test", rep_author = "")
  f = output$report_pptx
  if (!is.na(keep)) file.copy(f, keep, overwrite = TRUE)
  d = officer::read_pptx(f); sz = officer::slide_size(d); s = officer::pptx_summary(d)
  chk(abs(sz$width - 13.333) < 0.01 && abs(sz$height - 7.5) < 0.01, sprintf("16:9 slides (%.2f x %.2f in)", sz$width, sz$height))
  n = length(d); cat(n, "slides\n")
  t1 = s$text[s$slide_id == 1]
  chk(any(t1 == "ElectroFast slides test") && any(grepl("^Electrolyte: 0.1 M KCl", t1)) && any(grepl("^Working electrode:", t1)), "title slide with the typed cell")
  chk(!any(grepl("Reference electrode|Counter electrode", t1)), "empty cell entries left out of the title slide")
  chk(any(s$text[s$slide_id == 2] == "Main points"), "main points on slide 2")
  titles = vapply(seq_len(n), function(i) s$text[s$slide_id == i][1], "")
  chk(grepl("^Settings and files used", titles[n]), "settings on the last slide")
  pics = unique(s$slide_id[s$content_type == "image"]); cat("plot slides:", length(pics), "\n")
  under = vapply(pics, function(i) any(grepl("^•  ", s$text[s$slide_id == i])), TRUE)
  chk(length(pics) > 0 && any(under), "a result under at least one plot")
  x = tempfile(); utils::unzip(f, exdir = x)
  xml = vapply(list.files(file.path(x, "ppt/slides"), "[.]xml$", full.names = TRUE), function(p) paste(readLines(p, warn = FALSE), collapse = ""), "")
  tb = unlist(regmatches(xml, gregexpr("<a:tbl>.*?</a:tbl>", xml, perl = TRUE)))
  runs = unlist(regmatches(tb, gregexpr("<a:r>.*?</a:r>", tb, perl = TRUE)))
  chk(length(runs) > 0 && all(grepl("sz=\"1300\"", runs)) && all(grepl("typeface=\"Calibri\"", runs)), sprintf("%d table runs all 13 pt Calibri", length(runs)))
  chk(!any(grepl("<a:t>[^<]*•", unlist(regmatches(tb, gregexpr("<a:t>[^<]*</a:t>", tb))))), "no bullets inside tables")
  st = s$text[s$slide_id == n]
  chk(!any(grepl("^Reading time|^Primary ion charge", st)), "no ISE settings without ISE files")
  chk(all(grepl('algn="l"', unlist(regmatches(tb, gregexpr("<a:pPr[^>]*>", tb))))), "table cells left-aligned")
  print(data.frame(slide = seq_len(n), title = substr(titles, 1, 80)), row.names = FALSE)
})
cat("slide checks passed\n")
