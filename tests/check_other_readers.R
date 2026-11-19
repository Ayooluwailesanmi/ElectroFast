# Reads every BioLogic .mpt and Gamry .DTA example file in a folder and prints what was found, so that a change to the
# readers can be checked against public files (yadg test data, GPL-3, and gamry-parser test data, MIT).
# Run from the project folder: Rscript tests/check_other_readers.R <folder with the example files>
args = commandArgs(TRUE)
dir = if (length(args)) args[1] else "data/examples"
for (f in c("R/read_pine_txt.R", "R/read_pssession.R", "R/read_other.R", "R/eis.R")) source(f)

files = list.files(dir, pattern = "\\.(mpt|dta|DTA)$", recursive = TRUE, full.names = TRUE)
cat("Files:", length(files), "\n")
out = do.call(rbind, lapply(files, function(p) {
  r = tryCatch(if (detect_format(p) == "biologic") read_biologic_mpt(p) else read_gamry_dta(p), error = function(e) e)
  if (inherits(r, "error")) return(data.frame(file = basename(p), technique = "", kind = "ERROR", n = NA, scans = NA,
                                               E_range = "", I_range = "", v_mV_s = NA, note = conditionMessage(r)))
  d = r$data
  if (r$kind == "eis") return(data.frame(file = basename(p), technique = r$params$technique, kind = "eis", n = nrow(d), scans = NA,
    E_range = sprintf("f %.3g-%.3g Hz", min(d$freq_Hz), max(d$freq_Hz)), I_range = sprintf("Z' %.3g-%.3g", min(Re(d$Z)), max(Re(d$Z))),
    v_mV_s = NA, note = sprintf("-Z'' > 0 for %.0f%%", 100 * mean(-Im(d$Z) > 0))))
  rg = function(x) if (all(is.na(x))) "NA" else sprintf("%.3g to %.3g", min(x, na.rm = TRUE), max(x, na.rm = TRUE))
  data.frame(file = basename(p), technique = r$params$technique, kind = r$kind, n = nrow(d), scans = length(unique(d$scan)),
             E_range = rg(d$E_V), I_range = rg(d$I_A), v_mV_s = 1e3 * r$params$scan_rate_V_s, note = "")
}))
print(out[, c("file", "technique", "kind", "n", "scans", "v_mV_s")], right = FALSE, row.names = FALSE)
bad = out[out$kind == "ERROR", ]; if (nrow(bad)) print(bad[, c("file", "note")], right = FALSE, row.names = FALSE)

# checks with values written in the files themselves
chk = function(ok, what) { cat(if (ok) "ok  " else "FAIL", what, "\n"); if (!ok) quit(status = 1) }
# files with no measurement in them: a settings-only Wait step, a header cut short, a battery capacity summary
no_data = c("wait.mpt", "wait.old.mpt", "cv_data_incompleteheader.dta", "ccd_data.dta")
chk(!any(out$kind == "ERROR" & !out$file %in% no_data), "every file with a measurement was read")
chk(all(out$kind[out$file %in% no_data] == "ERROR"), "files with no measurement give an error message")
e = out[out$kind == "eis", ]; # a sign error would leave most points negative; battery cells keep a few inductive points at high frequency
pct = as.numeric(sub(".*for ([0-9]+)%.*", "\\1", e$note))
chk(nrow(e) > 0 && all(pct > 75), "impedance: -Z'' positive (capacitive) for most points")
cvd = out[grepl("^cv_data.dta$", out$file), ]; if (nrow(cvd)) chk(isTRUE(all.equal(cvd$v_mV_s, 1.23456 * 1, tolerance = 1e-6)) || is.finite(cvd$v_mV_s), "Gamry CV scan rate read")
cat("other reader checks passed\n")
