# The "Name on figures" box: ticked (the default), plots carry the small "ElectroFast" caption; unticked, no caption
# at all, also in a plot built for a download. Outside the app (scripts) the caption stays.
# Run from the project folder: Rscript tests/check_fig_credit.R
library(shiny); app = shinyAppDir(".")
chk = function(ok, what) { cat(if (ok) "ok  " else "FAIL", what, "\n"); if (!ok) quit(status = 1) }
ex = tempfile(); utils::unzip("www/examples/ElectroFast_example_data_simulated.zip", exdir = ex)
f = file.path(ex, "ElectroFast_simulated_data", "qc_history_2_quantities.csv")
testServer(app, {
  session$setInputs(qc_file = data.frame(name = basename(f), size = file.size(f), type = "text/csv", datapath = f), qc_nref = 12)
  chk(identical(qc_plot_obj()$labels$caption, "ElectroFast"), "ticked by default: caption present")
  session$setInputs(fig_credit = FALSE)
  chk(is.null(qc_plot_obj()$labels$caption), "unticked: no caption")
  chk(is.null(isolate(qc_plot_obj())$labels$caption), "unticked: no caption outside a reactive context (download)")
  session$setInputs(fig_credit = TRUE)
  chk(identical(qc_plot_obj()$labels$caption, "ElectroFast"), "ticked again: caption back")
})
source("R/app_core.R")
g = ggplot2::ggplot() + watermark()
chk(identical(g$labels$caption, "ElectroFast"), "outside the app: caption kept")
cat("figure credit checks passed\n")
