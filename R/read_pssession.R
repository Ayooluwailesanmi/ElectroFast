# Reader for PalmSens PSTrace session files (.pssession). A session is a JSON document
# stored as UTF-16 text (checked on PSTrace 5.11 files, October 2026). Each measurement
# holds a DataSet of arrays; each array names its quantity and its unit, and the unit
# type carries the scale (for example PalmSens.Units.MicroAmpere), so currents are
# converted to A from the unit given in the file, never assumed.
# Voltammetry and amperometry come back in the same form as read_chi_txt() (data with
# E_V, I_A, t_s, Q_C and a scan column; params with the scan rate and the method
# settings), so every technique tab can use them. Impedance comes back as freq_Hz and Z
# (complex, imaginary part as measured: the file stores -Z'').

pssession_scale = function(type) {
  s = c(PicoAmpere = 1e-12, NanoAmpere = 1e-9, MicroAmpere = 1e-6, MilliAmpere = 1e-3, Ampere = 1,
        NanoCoulomb = 1e-9, MicroCoulomb = 1e-6, MilliCoulomb = 1e-3, Coulomb = 1,
        MicroVolt = 1e-6, MilliVolt = 1e-3, Volt = 1, Time = 1, Hertz = 1, ZRe = 1, ZIm = 1, Z = 1, Phase = 1)
  k = sub("^PalmSens[.]Units[.]", "", type)
  if (k %in% names(s)) s[[k]] else NA_real_   # unknown scale: the array is kept out, see get() below
}

read_pssession_json = function(path) {
  b = readBin(path, "raw", file.info(path)$size)
  txt = if (length(b) >= 2 && b[1] == as.raw(0xFF) && b[2] == as.raw(0xFE)) iconv(list(b), "UTF-16LE", "UTF-8") else rawToChar(b)
  txt = gsub("﻿", "", txt, fixed = TRUE)
  txt = substr(txt, 1, max(gregexpr("}", txt, fixed = TRUE)[[1]]))   # PSTrace adds a mark after the JSON object
  jsonlite::fromJSON(txt, simplifyVector = FALSE)
}

method_settings = function(m) {
  l = strsplit(m, "\r?\n")[[1]]; l = l[grepl("^[A-Z_0-9]+=", l)]
  stats::setNames(as.list(sub("^[^=]+=", "", l)), sub("=.*$", "", l))
}

read_pssession = function(path, name = basename(path), measurement = 1) {
  j = read_pssession_json(path)
  if (!identical(j$Type, "PalmSens.DataFiles.SessionFile")) stop(name, " is not a PalmSens session file")
  me = j$Measurements[[measurement]]
  st = method_settings(me$Method %||% j$MethodForMeasurement)
  arr = lapply(me$DataSet$Values, function(v) list(desc = v$Description %||% "", q = v$Unit$Q %||% "", type = v$Unit$Type %||% "",
    values = vapply(v$DataValues, function(d) if (is.null(d$V)) NA_real_ else as.numeric(d$V), 0) * pssession_scale(v$Unit$Type %||% "")))
  qs = vapply(arr, `[[`, "", "q")
  # an array is used only when its unit scale is known; an unknown scale stops the read
  # rather than returning values that could be off by a factor of a thousand
  get = function(q, d = NULL) {
    k = which(qs == q & (is.null(d) | vapply(arr, `[[`, "", "desc") %in% d)); if (!length(k)) return(NULL)
    v = arr[[k[1]]]; if (anyNA(v$values) && all(is.na(v$values))) stop(name, ": unknown unit '", v$type, "' for ", q)
    v$values }
  params = list(file = name, technique = st$METHOD_ID %||% me$Title, title = me$Title, settings = st,
                scan_rate_V_s = suppressWarnings(as.numeric(st$SCAN_RATE %||% NA)), date = me$TimeStamp)
  if (identical(st$METHOD_ID, "eis") || any(qs == "Z'")) {
    zi = get("-Z''"); zi = if (!is.null(zi)) -zi else get("Z''")
    d = data.frame(freq_Hz = get("Frequency"), Z = complex(real = get("Z'"), imaginary = zi))
    d = d[is.finite(d$freq_Hz) & d$freq_Hz > 0, ]
    return(list(file = name, kind = "eis", params = params, data = d[order(-d$freq_Hz), ], note = ""))
  }
  # voltammetry and amperometry: one block of arrays per scan, in the order stored
  desc = vapply(arr, `[[`, "", "desc")
  # a scan is a description shared by a potential and a current array; chronopotentiometry and
  # other single-curve runs label each array by its own quantity ("potential", "current"), so they form one block
  scans = Filter(function(s) all(c("Potential", "Current") %in% qs[desc == s]), unique(desc[qs %in% c("Potential", "Current")]))
  if (!length(scans)) scans = list(NULL)
  t_all = get("Time")
  blocks = lapply(seq_along(scans), function(i) {
    s = scans[[i]]; E = get("Potential", s); I = get("Current", s); Q = get("Charge", s)
    n = max(length(E), length(I))
    t = if (length(t_all) == n) t_all else rep(NA_real_, n)
    data.frame(scan = i, E_V = if (length(E)) E else NA, I_A = I, t_s = t, Q_C = if (length(Q) == n) Q else NA)
  })
  d = do.call(rbind, blocks)
  # sweep segments (a new one at every change of direction), as in the CHI reader, for the CV analyses
  dir = sign(diff(d$E_V)); dir[dir == 0] = NA; dir = c(dir[1], dir)
  for (i in seq_along(dir)[-1]) if (is.na(dir[i])) dir[i] = dir[i - 1]
  d$segment = cumsum(c(1, diff(dir) != 0))
  # time runs from zero in every scan; make it continuous so the whole run can be plotted against time
  if (length(blocks) > 1 && !anyNA(d$t_s)) { per = vapply(blocks, function(b) max(b$t_s, na.rm = TRUE), 0); d$t_s = d$t_s + c(0, cumsum(per))[d$scan] }
  list(file = name, kind = "voltammetry", params = params, data = d)
}

# A session file is recognised by its UTF-16 mark and the PalmSens type
is_pssession = function(path) {
  b = readBin(path, "raw", 200)
  length(b) >= 2 && b[1] == as.raw(0xFF) && b[2] == as.raw(0xFE) &&
    grepl("PalmSens", iconv(list(b[3:length(b)]), "UTF-16LE", "UTF-8", sub = ""), fixed = TRUE)
}
