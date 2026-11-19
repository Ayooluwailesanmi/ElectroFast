# Session restore in a real (headless) Chrome against the local app on port 8421: no bar on a fresh page, file and
# setting saved in the browser, the warning before leaving, Restore brings both back (and the server reads the file
# again), a second reopen offers the restored session, Discard clears everything. Each stage is a check: the script
# stops with status 1 at the first one that fails. Start the app first, then run from the project folder:
# Rscript tests/check_restore.R
if (!isTRUE(tryCatch(nchar(paste(readLines(url("http://127.0.0.1:8421/"), warn = FALSE), collapse = "")) > 0, error = function(e) FALSE))) { cat("Local app not running on port 8421: test skipped\n"); quit(status = 0) }
out = tempdir()
library(chromote)
set.seed(5); t = seq(0, 20, by = 1e-3); I = 1e-9 + rnorm(length(t), 0, 2e-12)
f = file.path(out, "restore_sim.csv"); write.csv(data.frame(`Time (s)` = t, `Current (A)` = I, check.names = FALSE), f, row.names = FALSE)
b = ChromoteSession$new(width = 1300, height = 900)
ev = function(js, p = FALSE) tryCatch(b$Runtime$evaluate(js, awaitPromise = p)$result$value, error = function(e) paste("ERR", conditionMessage(e)))
dialogs = character(0); b$Page$enable()
b$Page$javascriptDialogOpening(callback = function(m) { dialogs <<- c(dialogs, m$type); b$Page$handleJavaScriptDialog(accept = TRUE, wait_ = FALSE) }, wait_ = FALSE)
wait = function(js, secs = 30) { t0 = Sys.time(); repeat { v = ev(js); if (isTRUE(v) || (is.character(v) && nzchar(v) && !startsWith(v, "ERR"))) return(v)
  if (difftime(Sys.time(), t0, units = "secs") > secs) return("TIMEOUT"); Sys.sleep(0.5) } }
click = function(txt) ev(sprintf("(() => { const a = [...document.querySelectorAll('a')].find(x => x.textContent.trim() === '%s'); if (!a) return 'missing'; a.click(); return 'ok' })()", txt))
connected = "(() => !!(window.Shiny && Shiny.shinyapp && Shiny.shinyapp.isConnected()))()"
chk = function(ok, what, got = "") {
  cat(if (isTRUE(ok)) "ok  " else "FAIL", what, sprintf("(%s)", paste(got, collapse = " | ")), "\n")
  if (!isTRUE(ok)) { try(b$close(), silent = TRUE); quit(status = 1) }
}

b$Page$navigate("http://127.0.0.1:8421/", wait_ = FALSE); wait(connected)
ev("(async () => { localStorage.clear(); await new Promise(r => { const q = indexedDB.deleteDatabase('electrofast_files'); q.onsuccess = q.onerror = q.onblocked = () => r(1); }); return 'cleared' })()", TRUE)
b$Page$reload(); Sys.sleep(1); wait(connected); Sys.sleep(2)
v = ev("!!document.getElementById('electrofast-restore-bar')"); chk(identical(v, FALSE), "1. fresh page: no restore bar", v)

click("Techniques"); click("Chronoamperometry"); Sys.sleep(1)
root = b$DOM$getDocument(depth = -1)$root
b$DOM$setFileInputFiles(files = list(normalizePath(f)), nodeId = b$DOM$querySelector(root$nodeId, "#ca_files")$nodeId)
v = wait("(() => { const e = document.querySelector('#ca_files_progress .progress-bar'); return e && e.textContent.includes('complete') ? 'yes' : '' })()")
chk(identical(v, "yes"), "2. file uploaded", v)
ev("(() => { $('#ca_qml').val(512).trigger('change'); return 1 })()"); Sys.sleep(2)
v = ev("(() => { const s = JSON.parse(localStorage.getItem('electrofast_session_v1')); return JSON.stringify({ files: s.files, ca_qml: s.inputs.ca_qml }) })()")
chk(grepl("\"ca_files\"", v) && grepl("512", v), "   file (by its upload box) and setting saved in this browser", v)

# a real mouse click gives the page the user interaction that browsers require before they show the leaving prompt
b$Input$dispatchMouseEvent(type = "mousePressed", x = 650, y = 600, button = "left", clickCount = 1)
b$Input$dispatchMouseEvent(type = "mouseReleased", x = 650, y = 600, button = "left", clickCount = 1)
ev("location.reload()"); Sys.sleep(3); wait(connected); Sys.sleep(2)
chk(length(dialogs) > 0, "3. warning before leaving shown", dialogs)
v = wait("(() => document.getElementById('electrofast-restore-bar') ? document.getElementById('electrofast-restore-bar').innerText : '')()")
chk(is.character(v) && grepl("Restore", v), "   restore bar offered after reopening", substr(v, 1, 60))
v = ev("document.getElementById('ca_qml').value"); chk(!identical(v, "512"), "   setting not restored before Restore is clicked", v)
ev("(() => { document.getElementById('electrofast-restore').click(); return 1 })()")
Sys.sleep(6)
v1 = ev("(() => { const e = document.querySelector('#ca_files').closest('.input-group').querySelector('input[type=text]'); return e ? e.value : 'none' })()")
v2 = ev("document.getElementById('ca_qml').value")
chk(grepl("restore_sim.csv", v1) && identical(v2, "512"), "4. after Restore: file and setting back", c(v1, v2))
click("Techniques"); click("Chronoamperometry"); Sys.sleep(1); click("File check")
v = wait("(() => { const t = document.querySelector('#ca_check table'); return t ? t.innerText.split('\\n')[1] : '' })()", 30)
chk(is.character(v) && grepl("restore_sim", v), "   server reads the restored file (file check table)", substr(v, 1, 60))
ev("location.reload()"); Sys.sleep(3); wait(connected); Sys.sleep(2)
v = ev("!!document.getElementById('electrofast-discard')"); chk(isTRUE(v), "5. restored session offered again on the next visit", v)
ev("(() => { document.getElementById('electrofast-discard').click(); return 1 })()"); Sys.sleep(2)
v1 = ev("localStorage.getItem('electrofast_session_v1_pending') === null")
v2 = ev("(async () => { const d = await new Promise(r => { const q = indexedDB.open('electrofast_files', 1); q.onsuccess = () => r(q.result); }); return await new Promise(r => { const c = d.transaction('files').objectStore('files').count(); c.onsuccess = () => r(c.result); }); })()", TRUE)
chk(isTRUE(v1) && identical(as.numeric(v2), 0), "   after Discard: kept copy and saved files removed", c(v1, v2))
ev("location.reload()"); Sys.sleep(3); wait(connected); Sys.sleep(2)
v = ev("!!document.getElementById('electrofast-restore-bar')"); chk(identical(v, FALSE), "   reopened after Discard: no restore bar", v)
b$close()
cat("restore checks passed\n")
