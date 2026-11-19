# Version and provenance: what made a result. Every report carries the tool version, the R version, the time it was
# made and a checksum of every file used, and the Report tab can download all of it as one bundle (settings, file
# checksums, package versions, every results table and a copy of the code) so a result can be traced and re-run.

# Pre-release number until the first public release; raise it with every deployed change
electrofast_version = "1.0.0"
electrofast_version_date = "2026-10-21"

electrofast_packages = c("shiny", "bslib", "ggplot2", "officer", "rmarkdown", "jsonlite", "lme4", "zip", "gridExtra")
package_versions = function(pk = electrofast_packages) {
  v = vapply(pk, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "not installed", "")
  data.frame(package = pk, version = unname(v))
}

# The uploaded files: every fileInput value is a data frame with name, size and datapath
input_files = function(inputs) {
  f = Filter(function(v) is.data.frame(v) && all(c("name", "datapath") %in% names(v)), inputs)
  if (!length(f)) return(data.frame(input = character(), file = character(), bytes = numeric(), md5 = character()))
  do.call(rbind, lapply(names(f), function(k) { d = f[[k]]
    data.frame(input = k, file = d$name, bytes = file.size(d$datapath), md5 = unname(tools::md5sum(d$datapath))) }))
}

# Every other setting, as text (buttons and file inputs left out)
input_settings = function(inputs) {
  s = Filter(function(v) !is.data.frame(v) && !inherits(v, "shinyActionButtonValue") && length(v) <= 50 && (is.atomic(v) || is.null(v)), inputs)
  s = s[order(names(s))]
  data.frame(input = names(s), value = vapply(s, function(v) if (is.null(v)) "" else paste(v, collapse = ", "), ""), row.names = NULL)
}

provenance_lines = function(files = NULL) c(
  `ElectroFast version` = sprintf("%s (%s)", electrofast_version, electrofast_version_date),
  `R version` = paste(R.version$major, R.version$minor, sep = "."),
  Created = format(Sys.time(), "%Y-%m-%d %H:%M %Z"),
  `Files used` = if (is.null(files) || !nrow(files)) "none" else paste(sprintf("%s (md5 %s)", files$file, substr(files$md5, 1, 8)), collapse = "; "))

# The bundle: a zip with settings, file checksums, versions, every results table of the report sections and the code
build_bundle = function(zipfile, inputs, sections, code_dir = ".") {
  d = file.path(tempfile("electrofast_bundle_"), "ElectroFast_bundle"); dir.create(file.path(d, "results"), recursive = TRUE)
  files = input_files(inputs)
  utils::write.csv(input_settings(inputs), file.path(d, "settings.csv"), row.names = FALSE)
  utils::write.csv(files, file.path(d, "files.csv"), row.names = FALSE)
  writeLines(c(paste(names(provenance_lines(files)), provenance_lines(files), sep = ": "), "", "Packages:",
               apply(package_versions(), 1, paste, collapse = " ")), file.path(d, "versions.txt"))
  k = 0
  for (s in sections) for (nm in names(s$tables)) { t = s$tables[[nm]]; if (!is.data.frame(t)) next; k = k + 1
    fn = sprintf("%02d_%s_%s.csv", k, gsub("[^A-Za-z0-9]+", "_", s$title), gsub("[^A-Za-z0-9]+", "_", nm))
    utils::write.csv(t, file.path(d, "results", fn), row.names = FALSE) }
  dir.create(file.path(d, "code", "R"), recursive = TRUE); dir.create(file.path(d, "code", "report"))
  file.copy(list.files(code_dir, "^app.*\\.R$", full.names = TRUE), file.path(d, "code"))
  file.copy(list.files(file.path(code_dir, "R"), "\\.R$", full.names = TRUE), file.path(d, "code", "R"))
  file.copy(file.path(code_dir, "report", "electrofast_report.Rmd"), file.path(d, "code", "report"))
  writeLines(c("ElectroFast reproducibility bundle", "",
    "settings.csv   every setting in the app when the bundle was made",
    "files.csv      every uploaded file, with its size and MD5 checksum (to show the same files are used again)",
    "versions.txt   ElectroFast, R and package versions, and when the bundle was made",
    "results/       every results table in the report, as CSV",
    "code/          the app code that made the results",
    "", "To repeat the analysis:", "1. Run the app with shiny::runApp(\"code\"), or open the same ElectroFast version.",
    "2. Load the files listed in files.csv and check their MD5 sums.", "3. Enter the settings in settings.csv."), file.path(d, "README.txt"))
  old = setwd(dirname(d)); on.exit(setwd(old))
  zip::zipr(zipfile, "ElectroFast_bundle")
}
