# Every tab and sub-tab of the running app (port 8421), opened in a real (headless) Chrome with no data loaded: no red
# error text may show, such as "Error: attempt to set an attribute on NULL". Validation messages (the grey "load
# files" kind) are allowed. Start the app first, then run from the project folder: Rscript tests/check_empty_tabs.R
if (!isTRUE(tryCatch(nchar(paste(readLines(url("http://127.0.0.1:8421/"), warn = FALSE), collapse = "")) > 0, error = function(e) FALSE))) {
  cat("Local app not running on port 8421: test skipped\n"); quit(status = 0) }
library(chromote); options(chromote.timeout = 60)
b = ChromoteSession$new(width = 1400, height = 900)
ev = function(js) tryCatch(b$Runtime$evaluate(js)$result$value, error = function(e) paste("ERR", conditionMessage(e)))
b$Page$navigate("http://127.0.0.1:8421/", wait_ = FALSE); Sys.sleep(8)
ev("(async () => { localStorage.clear(); return 1 })()")
# every link that switches a tab: top bar, menu entries and sub-tabs
tabs = ev("JSON.stringify([...document.querySelectorAll('a[data-bs-toggle=tab], a[data-toggle=tab]')].map(a => a.getAttribute('data-value')).filter(v => v))")
tabs = unique(jsonlite::fromJSON(tabs))
bad = character(0)
for (v in tabs) {
  ev(sprintf("(() => { const a = [...document.querySelectorAll('a[data-bs-toggle=tab], a[data-toggle=tab]')].find(x => x.getAttribute('data-value') === %s); if (a) a.click(); return 1 })()", jsonlite::toJSON(v, auto_unbox = TRUE)))
  Sys.sleep(1.2)
  e = ev("JSON.stringify([...document.querySelectorAll('.shiny-output-error:not(.shiny-output-error-validation)')].filter(x => x.offsetParent !== null && x.innerText.trim()).map(x => x.id + ': ' + x.innerText.trim().slice(0, 120)))")
  e = jsonlite::fromJSON(e); if (length(e)) bad = c(bad, sprintf("[%s] %s", v, e))
}
b$close()
cat(length(tabs), "tabs opened\n")
if (length(bad)) { cat("FAIL red error text with no data loaded:\n", paste(" ", unique(bad), collapse = "\n"), "\n"); quit(status = 1) }
cat("ok   no red error text on any tab with no data loaded\nempty-tab checks passed\n")
