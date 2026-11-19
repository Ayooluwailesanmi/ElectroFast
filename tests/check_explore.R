# Exploration tools in the Chronoamperometry tab: file check, zoom, event explorer and file comparison, on simulated
# current-time files with known content (written as text files). Run from the project folder: Rscript tests/check_explore.R
source("R/explore.R"); source("R/steps.R")
chk = function(ok, what) { cat(if (isTRUE(ok)) "ok  " else "FAIL", what, "\n"); if (!isTRUE(ok)) quit(status = 1) }

# functions on known traces
set.seed(5); t = seq(0, 60, by = 1e-3); I = 1e-9 + rnorm(length(t), 0, 2e-12) + 1e-10 * (t > 30)
fc = file_check(list(file = "sim", data = data.frame(t_s = t, I_A = I), params = list()))
chk(abs(fc$noise_A - 2e-12) / 2e-12 < 0.05, "file check: noise of 2 pA recovered within 5%, the 100 pA step not counted as noise")
chk(abs(fc$net_change_A - 1e-10) < 2e-12 && fc$gaps == 0, "file check: net change of 100 pA, no gaps")
cl = file_check(list(file = "clipped", data = data.frame(t_s = t, I_A = pmin(I, 1.05e-9)), params = list()))
chk(cl$at_limit_pct > 10 && length(file_check_notes(cl)) == 1, "file check: a clipped trace is flagged")
gp = file_check(list(file = "gap", data = data.frame(t_s = t[-(20000:20100)], I_A = I[-(20000:20100)]), params = list()))
chk(gp$gaps == 1 && any(grepl("gap", file_check_notes(gp))), "file check: a gap in time is flagged")
chk(length(file_check_notes(fc)) == 0, "file check: nothing flagged on a clean trace")
er = event_rate(c(1, 2, 15, 16, 17, 55), 0, 60, 10)
chk(identical(er$count, c(2L, 3L, 0L, 0L, 0L, 1L)), "event rate: counts per 10 s bin")
k = trace_window(t, 10, 10.5); chk(abs(t[min(k)] - 10) < 1e-9 && abs(t[max(k)] - 10.5) < 1e-9, "zoom window: the points between the two times")
blank = 1e-9 + rnorm(length(t), 0, 2e-12)
cmp = compare_traces(list(list(file = "sample", data = data.frame(t_s = t, I_A = I)), list(file = "blank", data = data.frame(t_s = t, I_A = blank))), ref = 2, n_bins = 60)
s = cmp[cmp$file == "sample", ]
chk(abs(mean(s$diff_A[s$t_s < 29])) < 2e-12 && abs(mean(s$diff_A[s$t_s > 31]) - 1e-10) < 2e-12, "comparison: sample minus blank is 0 before and 100 pA after the step")

# through the app, with the two traces as text files
d = tempfile(); dir.create(d)
write.csv(data.frame(`Time (s)` = t, `Current (A)` = I, check.names = FALSE), file.path(d, "sample.csv"), row.names = FALSE)
write.csv(data.frame(`Time (s)` = t, `Current (A)` = blank, check.names = FALSE), file.path(d, "blank.csv"), row.names = FALSE)
f = file.path(d, c("sample.csv", "blank.csv"))
library(shiny)
suppressMessages(testServer(shinyAppDir("."), {
  session$setInputs(ca_files = data.frame(name = basename(f), size = file.size(f), type = "text/csv", datapath = normalizePath(f)), instrument = "auto",
                    ca_mode = "esteps", ca_sfrom = 0, ca_sto = 1e9, ca_ebin = 10, ca_ewin = 0.5, ca_ek = 6, ca_edir = "both", ca_xbin = 10,
                    ca_zfile = "sample.csv", ca_ref = "blank.csv", ca_cmp0 = FALSE, `ca_el-shape` = "preset", `ca_el-preset` = "2")
  chk(nrow(ca_check_df()) == 2 && all(ca_check_df()$gaps == 0), "app: file check for both files"); invisible(output$ca_check); invisible(output$ca_check_notes)
  g = ca_zoom_obj(); chk(grepl("every point shown", g$labels$subtitle) == FALSE, "app: whole 60 s trace thinned for the zoom view")
  session$setInputs(ca_zbrush = list(xmin = 29.9, xmax = 30.1)); g = ca_zoom_obj()
  chk(grepl("every point shown", g$labels$subtitle) && grepl("29.9", g$labels$subtitle), "app: zooming to 0.2 s shows every point")
  pts = ggplot2::ggplot_build(g)$data; chk(length(pts) == 2 && nrow(pts[[2]]) == 1, "app: the detected step is marked in the zoom")
  session$setInputs(ca_zdbl = list(x = 1, y = 1)); chk(is.null(ca_zr()), "app: double-click shows the whole trace again")
  e = ca_xev(); chk(sum(e$file == "sample.csv") == 1 && sum(e$file == "blank.csv") == 0, "app: event explorer has the one step, in the sample only")
  for (o in c("ca_x_rate", "ca_x_amp_t", "ca_x_hist")) invisible(output[[o]])
  chk(inherits(tryCatch(ca_x_amp_dur_obj(), error = function(e) e), "error"), "app: amplitude against duration is explained, not drawn, for steps")
  dd = ca_cmp(); s = dd[dd$file == "sample.csv", ]
  chk(abs(mean(s$diff_A[s$t_s > 31]) - 1e-10) < 2e-12, "app: comparison shows the 100 pA difference from the blank"); invisible(output$ca_cmp_plot); invisible(output$ca_cmp_diff)
  session$setInputs(ca_mode = "multistep"); chk(inherits(tryCatch(ca_xev(), error = function(e) e), "error"), "app: event explorer asks for the spikes or steps analysis")
}))
cat("exploration checks passed\n")
