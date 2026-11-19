if (!exists("%||%")) `%||%` = function(a, b) if (is.null(a)) b else a
# Readers for BioLogic EC-Lab text exports (.mpt), Gamry Framework files (.DTA) and any other delimited text
# export (for example Metrohm Autolab NOVA ASCII), returned in the same form as read_chi_txt(): data with
# E_V, I_A, t_s, Q_C, scan and segment; params with the technique, the scan rate (V/s) and the header settings.
# Units are read from each column name (A, mA, uA, nA; V, mV; s, min, h; C, mC, mA.h), never assumed.
# Layouts checked on the public example files of yadg (BioLogic, GPL-3) and gamry-parser (Gamry, MIT), October 2026.

unit_scale = function(u) {
  u = gsub("\u00b5|\u03bc", "u", trimws(u)); u = sub("^\\((.*)\\)$", "\\1", u)
  s = c(A = 1, mA = 1e-3, uA = 1e-6, nA = 1e-9, pA = 1e-12, V = 1, mV = 1e-3, uV = 1e-6, s = 1, ms = 1e-3, min = 60, h = 3600,
        C = 1, mC = 1e-3, uC = 1e-6, nC = 1e-9, "mA.h" = 3.6, "A.h" = 3600, "uA.h" = 3.6e-3, Hz = 1, kHz = 1e3, MHz = 1e6,
        Ohm = 1, ohm = 1, kOhm = 1e3, kohm = 1e3, MOhm = 1e6, Mohm = 1e6, "V vs. Ref." = 1, "V vs. Ref" = 1)
  if (u %in% names(s)) s[[u]] else NA_real_
}
# name and unit of a column: "Ewe/V", "<I>/mA", "Time (s)", "WE(1).Current (A)"
split_col = function(nm) {
  if (grepl("/", nm) && !grepl("\\(.*/.*\\)", nm)) return(c(name = sub("/[^/]*$", "", nm), unit = sub("^.*/", "", nm)))
  m = regmatches(nm, regexec("^(.*?)\\s*\\(([^()]*)\\)\\s*$", nm))[[1]]
  if (length(m) == 3) c(name = m[2], unit = m[3]) else c(name = nm, unit = "")
}
num_col = function(v, comma = FALSE) { v = trimws(as.character(v)); if (comma) v = sub(",", ".", v, fixed = TRUE); suppressWarnings(as.numeric(v)) }
add_segments = function(d) {
  if (all(is.na(d$E_V))) { d$segment = 1L; return(d) }
  dir = sign(diff(d$E_V)); dir[dir == 0] = NA; dir = c(dir[1], dir)
  for (i in seq_along(dir)[-1]) if (is.na(dir[i])) dir[i] = dir[i - 1]
  if (all(is.na(dir))) dir[] = 1
  d$segment = cumsum(c(1, diff(dir) != 0)); d
}
# pick the first column whose name matches a pattern, and convert it to SI with its own unit
pick = function(tab, nms, patterns, comma = FALSE, unit_default = NA) {
  for (p in patterns) { k = which(grepl(p, nms, ignore.case = TRUE, perl = TRUE)); if (length(k)) {
    nu = split_col(nms[k[1]]); f = unit_scale(nu[["unit"]]); if (is.na(f)) f = unit_default
    if (is.na(f)) stop("column '", nms[k[1]], "': unknown unit '", nu[["unit"]], "'")
    return(list(v = num_col(tab[[k[1]]], comma) * f, name = nms[k[1]])) } }
  NULL
}

# BioLogic EC-Lab text export (.mpt): "EC-Lab ASCII FILE", "Nb header lines : n", header text, tab-separated table
read_biologic_mpt = function(path, name = basename(path)) {
  L = iconv(readLines(path, warn = FALSE, encoding = "latin1"), "latin1", "UTF-8")
  n = as.integer(gsub("[^0-9]", "", L[2]))
  hdr = L[seq_len(n - 1)]; cols = strsplit(L[n], "\t")[[1]]; cols = cols[nzchar(trimws(cols))]
  if (length(L) <= n) stop(name, ": the file holds settings only, no data rows (technique: ", trimws(L[4]), ")")
  # trailing empty fields are dropped by strsplit, so short rows are padded with NA
  rows = strsplit(L[(n + 1):length(L)], "\t"); rows = rows[lengths(rows) >= 2]
  tab = as.data.frame(do.call(rbind, lapply(rows, `[`, seq_along(cols))), stringsAsFactors = FALSE); names(tab) = cols
  comma = any(grepl("^-?[0-9]+,[0-9]", tab[[1]][seq_len(min(5, nrow(tab)))])) || any(grepl("[0-9],[0-9]", unlist(tab[1, ])))
  kv = regmatches(hdr, regexec("^(\\S.*?)\\s{2,}(\\S.*)$", hdr)); kv = kv[lengths(kv) == 3]
  st = stats::setNames(lapply(kv, function(m) trimws(m[3])), vapply(kv, function(m) trimws(m[2]), ""))
  # several values when the technique has several sequences: the first is kept, as the scan rate of the first sweep
  first_tok = function(x) strsplit(trimws(x), "\\s+")[[1]][1]
  v = if (!is.null(st[["dE/dt"]])) num_col(first_tok(st[["dE/dt"]]), TRUE) * (if (grepl("^mV", first_tok(st[["dE/dt unit"]] %||% "mV/s"))) 1e-3 else 1) else NA
  if (isTRUE(v == 0)) v = NA
  # exports saved without the settings block have no technique line
  tech = trimws(hdr[4]); if (is.na(tech)) tech = "" else if (!nzchar(tech) || grepl(":", tech)) tech = trimws(hdr[grepl("^[A-Z][a-z]", hdr) & !grepl(":", hdr)][1] %||% "")
  params = list(file = name, technique = tech, settings = st, scan_rate_V_s = v, date = sub("^.*: ", "", grep("Acquisition started on", hdr, value = TRUE)[1] %||% ""))
  # impedance when Re(Z) is there and most rows have a frequency; Modulo Bat runs mix DC steps with a few
  # impedance points and stay as voltammetry
  fq = pick(tab, cols, "^freq/", comma)
  if (!is.null(fq) && any(grepl("^Re\\(Z\\)/", cols)) && mean(fq$v > 0, na.rm = TRUE) > 0.5) {
    f = fq$v; re = pick(tab, cols, "^Re\\(Z\\)/", comma)$v
    im = pick(tab, cols, "^-Im\\(Z\\)/", comma); im = if (!is.null(im)) -im$v else pick(tab, cols, "^Im\\(Z\\)/", comma)$v
    d = data.frame(freq_Hz = f, Z = complex(real = re, imaginary = im)); d = d[is.finite(d$freq_Hz) & d$freq_Hz > 0, ]
    return(list(file = name, kind = "eis", params = params, data = d, note = ""))
  }
  # battery runs (BCD) call the same potential Ecell
  E = pick(tab, cols, c("^Ewe/", "^<Ewe>/", "^Ecell/", "^Ewe-Ece/", "^control/V$"), comma); I = pick(tab, cols, c("^<I>/", "^I/"), comma)
  t = pick(tab, cols, "^time/", comma); Q = pick(tab, cols, c("^\\(Q-Qo\\)/", "^Q-Qo/"), comma); cyc = pick(tab, cols, "^cycle number$", comma, 1)
  d = data.frame(E_V = if (is.null(E)) NA else E$v, I_A = if (is.null(I)) NA else I$v, t_s = if (is.null(t)) NA else t$v - t$v[1],
                 Q_C = if (is.null(Q)) NA else Q$v, scan = if (is.null(cyc)) 1L else as.integer(cyc$v))
  d = add_segments(d[is.finite(d$t_s) | is.finite(d$E_V), ])
  list(file = name, kind = "voltammetry", params = params, data = d)
}

# Gamry Framework file (.DTA): "EXPLAIN", "TAG <technique>", parameter lines, then tables named CURVE, CURVEn or ZCURVE,
# each with a column-name row and a unit row that start with a tab
read_gamry_dta = function(path, name = basename(path)) {
  L = readLines(path, warn = FALSE, encoding = "latin1"); L = iconv(L, "latin1", "UTF-8")
  tag = sub("^TAG\t", "", grep("^TAG\t", L, value = TRUE)[1])
  par = strsplit(grep("^[A-Z][A-Z0-9_]*\t(QUANT|POTEN|IQUANT|SELECTOR|TOGGLE|LABEL|STRING|TWOPARAM)\t", L, value = TRUE), "\t")
  st = stats::setNames(lapply(par, function(x) x[3]), vapply(par, `[`, "", 1))
  v = if (!is.null(st[["SCANRATE"]])) as.numeric(st[["SCANRATE"]]) * (if (any(grepl("^SCANRATE\t.*mV/s", L))) 1e-3 else 1) else NA
  starts = grep("^[A-Z]*CURVE[0-9]*\tTABLE", L)
  tables = lapply(starts, function(s) {
    cols = strsplit(sub("^\t", "", L[s + 1]), "\t")[[1]]; units = strsplit(sub("^\t", "", L[s + 2]), "\t")[[1]]
    e = s + 3; while (e <= length(L) && startsWith(L[e], "\t")) e = e + 1
    rows = strsplit(sub("^\t", "", L[(s + 3):(e - 1)]), "\t")
    tab = as.data.frame(do.call(rbind, lapply(rows, `[`, seq_along(cols))), stringsAsFactors = FALSE); names(tab) = cols
    list(name = sub("\t.*$", "", L[s]), tab = tab, units = stats::setNames(units[seq_along(cols)], cols))
  })
  params = list(file = name, technique = tag, settings = st, scan_rate_V_s = v, date = st[["DATE"]] %||% "")
  # files saved with a European locale write decimal commas
  col = function(tb, k) { if (!k %in% names(tb$tab)) return(NULL); f = unit_scale(tb$units[[k]]); if (is.na(f)) f = 1
    x = tb$tab[[k]]; num_col(x, any(grepl("^-?[0-9]+,[0-9]", utils::head(x, 5)))) * f }
  z = Filter(function(tb) tb$name == "ZCURVE", tables)
  if (length(z)) { tb = z[[1]]
    d = data.frame(freq_Hz = col(tb, "Freq"), Z = complex(real = col(tb, "Zreal"), imaginary = col(tb, "Zimag")))
    return(list(file = name, kind = "eis", params = params, data = d[is.finite(d$freq_Hz) & d$freq_Hz > 0, ], note = "")) }
  cur = Filter(function(tb) grepl("^CURVE[0-9]*$", tb$name), tables)
  if (!length(cur)) cur = Filter(function(tb) tb$name != "ZCURVE" && any(c("Vf", "Im", "Voltage", "Current") %in% names(tb$tab)), tables)
  if (!length(cur)) stop(name, ": no potential or current table found (tables: ", paste(vapply(tables, `[[`, "", "name"), collapse = ", "), ")")
  d = do.call(rbind, lapply(seq_along(cur), function(i) { tb = cur[[i]]
    # pulse techniques: base potential and difference current; otherwise the measured potential and current
    E = col(tb, "Vstep") %||% col(tb, "Vf") %||% col(tb, "Vu") %||% col(tb, "Voltage"); I = col(tb, "Idif") %||% col(tb, "Im") %||% col(tb, "Current")
    n = nrow(tb$tab)
    data.frame(E_V = if (is.null(E)) rep(NA_real_, n) else E, I_A = if (is.null(I)) rep(NA_real_, n) else I,
               t_s = col(tb, "T") %||% rep(NA_real_, n), Q_C = col(tb, "Q") %||% rep(NA_real_, n), scan = i) }))
  # time restarts in every CURVE table of a CV: make it continuous
  if (length(cur) > 1 && !anyNA(d$t_s)) { last = tapply(d$t_s, d$scan, max); d$t_s = d$t_s + c(0, cumsum(last))[d$scan] }
  list(file = name, kind = "voltammetry", params = params, data = add_segments(d))
}

# Any other delimited text export with a header row: columns are found by their names and units
read_generic_text = function(path, name = basename(path)) {
  L = readLines(path, warn = FALSE, encoding = "UTF-8"); L = L[nzchar(trimws(L))]
  sep = names(which.max(c("\t" = sum(grepl("\t", L)), ";" = sum(grepl(";", L)), "," = sum(grepl(",", L)))))
  isnum = function(l) { x = trimws(strsplit(l, sep, fixed = TRUE)[[1]]); x = x[nzchar(x)]; length(x) >= 2 && mean(!is.na(suppressWarnings(as.numeric(sub(",", ".", x))))) > 0.8 }
  if (length(L) <= 5e4) {
    # the same test as isnum(), on all lines at once (line by line took 12 s for a 30000-row file): a line is data when
    # it has 2 or more non-empty fields and more than 80% of them are numbers
    toks = strsplit(L, sep, fixed = TRUE); id = rep(seq_along(toks), lengths(toks)); x = trimws(unlist(toks, use.names = FALSE))
    keep = nzchar(x); id = id[keep]; x = x[keep]
    okn = !is.na(suppressWarnings(as.numeric(sub(",", ".", x))))
    cnt = tabulate(id, nbins = length(L)); good = tabulate(id[okn], nbins = length(L))
    num = cnt >= 2 & good > 0.8 * cnt; first = which(num)[1]
    if (is.na(first) || first < 2) stop(name, ": no header row followed by numbers was found")
    cols = trimws(gsub('"', "", strsplit(L[first - 1], sep, fixed = TRUE)[[1]]))
    rows = toks[first:length(L)][num[first:length(L)]]
    tab = as.data.frame(do.call(rbind, lapply(rows, `[`, seq_along(cols))), stringsAsFactors = FALSE); names(tab) = cols
  } else {
    # long recordings (a 10 kHz trace has millions of rows): find the header in the first lines, then read the table
    # in one pass; checking every line on its own took minutes
    num = vapply(L[1:2000], isnum, TRUE); first = which(num)[1]
    if (is.na(first) || first < 2) stop(name, ": no header row followed by numbers was found")
    cols = trimws(gsub('"', "", strsplit(L[first - 1], sep, fixed = TRUE)[[1]]))
    tab = utils::read.table(text = L[first:length(L)], sep = sep, header = FALSE, colClasses = "character", fill = TRUE,
                            quote = "\"", comment.char = "", col.names = paste0("V", seq_len(max(length(cols), 2))), flush = TRUE)
    tab = tab[, seq_along(cols), drop = FALSE]; names(tab) = cols
    tab = tab[!is.na(suppressWarnings(as.numeric(sub(",", ".", tab[[1]])))), , drop = FALSE]   # rows of numbers only, as above
  }
  comma = sep != "," && any(grepl("[0-9],[0-9]", tab[[1]]))
  f = pick(tab, cols, c("^freq"), comma)
  if (!is.null(f)) return(NULL)                     # impedance: read_eis() handles it
  E = pick(tab, cols, c("WE\\(1\\)\\.Potential", "^Ewe", "Potential(?! applied)", "^E\\b", "^Voltage", "Potential applied"), comma, 1)
  I = pick(tab, cols, c("WE\\(1\\)\\.Current", "^<?I>?\\b", "Current"), comma, 1)
  t = pick(tab, cols, c("^Time", "^t\\b"), comma, 1); Q = pick(tab, cols, c("Charge", "^Q\\b"), comma, 1)
  sc = pick(tab, cols, c("^Scan$", "^Cycle"), comma, 1)
  if (is.null(E) && is.null(I)) stop(name, ": no potential or current column was recognised (columns: ", paste(cols, collapse = ", "), ")")
  n = nrow(tab)
  d = data.frame(E_V = E$v %||% rep(NA_real_, n), I_A = I$v %||% rep(NA_real_, n), t_s = t$v %||% rep(NA_real_, n), Q_C = Q$v %||% rep(NA_real_, n),
                 scan = if (is.null(sc)) 1L else as.integer(sc$v))
  v = NA
  if (!all(is.na(d$t_s)) && !all(is.na(d$E_V))) { k = which(diff(d$t_s) > 0 & is.finite(diff(d$E_V))); if (length(k) > 5) v = stats::median(abs(diff(d$E_V)[k] / diff(d$t_s)[k])) }
  list(file = name, kind = "voltammetry", params = list(file = name, technique = "text export", settings = list(columns = paste(cols, collapse = "; ")),
       scan_rate_V_s = if (is.finite(v) && v > 0) signif(v, 3) else NA), data = add_segments(d))
}

# which reader a file needs, from its first lines
detect_format = function(path) {
  if (is_pssession(path)) return("pssession")
  if (exists("is_mpr") && is_mpr(path)) return("mpr")
  if (exists("is_abf") && is_abf(path)) return("abf")
  L = readLines(path, n = 3, warn = FALSE)
  if (grepl("^EC-Lab ASCII FILE", L[1])) return("biologic")
  if (grepl("^EXPLAIN", L[1])) return("gamry")
  if (grepl("Plot Name:", L[1])) return("pine")
  # CHI starts with the date, the month written in full or shortened ("July 2, 2026", "Sept. 9, 2026")
  if (grepl("^[A-Z][a-z]+\\.? +[0-9]{1,2}, +[0-9]{4}", L[1])) return("chi")
  # numbers only, no header (CHI "data only" export): the CHI reader knows its column layouts
  if (grepl("^\\s*[-+]?[0-9.]+([eE][-+]?[0-9]+)?([\t,]\\s*[-+]?[0-9.]+([eE][-+]?[0-9]+)?)+\\s*$", L[1])) return("chi")
  "text"
}

# The instrument the user chose in the top bar ("auto" = recognise it from the file). A file that does not match the
# chosen instrument stops with a message saying what it looks like; "text" reads any delimited file with a header.
instrument_labels = c(auto = "Auto-detect", chi = "CH Instruments", pssession = "PalmSens", biologic = "BioLogic EC-Lab",
                      gamry = "Gamry", pine = "Pine AfterMath", abf = "Axon pCLAMP (.abf)", text = "Metrohm Autolab or other text/CSV")
instrument_format = function(path, name, instrument = "auto") {
  det = detect_format(path)
  fam = c(chi = "chi", pssession = "pssession", mpr = "biologic", biologic = "biologic", gamry = "gamry", pine = "pine", abf = "abf", text = "text")[[det]]
  if (is.null(instrument) || instrument == "auto" || instrument == fam) return(det)
  if (instrument == "text" && !det %in% c("pssession", "mpr", "abf")) return("text")
  stop(sprintf("%s does not look like a %s file (it looks like %s). Choose Auto-detect, or the instrument that recorded it.",
               name, instrument_labels[[instrument]], instrument_labels[[fam]]))
}
