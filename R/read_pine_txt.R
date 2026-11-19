# Reader for Pine Research AfterMath text exports (WaveDriver bipotentiostats).
# The export puts every trace side by side, four columns each (Index, Segment, x, y),
# under a "Trace Name:" row and a units row. Each trace appears twice, once as a trace
# and once as a result, and the two copies can use different units (mA in one, µA in
# the other), and so can different files of the same experiment. Units are therefore
# always read from the header, never assumed.
# K1 is the first working electrode (the disk of an RRDE), K2 the second (the ring).
# The measured traces ("vs Time") are used; "Applied Potential" is the programmed
# waveform at a much finer time step and is not a measurement.

pine_units = c("V" = 1, "mV" = 1e-3, "uV" = 1e-6,
               "A" = 1, "mA" = 1e-3, "uA" = 1e-6, "nA" = 1e-9, "pA" = 1e-12,
               "s" = 1, "ms" = 1e-3)

read_pine_txt = function(path) {
  lines = readLines(path, warn = FALSE, encoding = "UTF-8")
  lines = sub("\r$", "", lines)
  lines[1] = sub("^[^P]*Plot", "Plot", lines[1])          # drops a byte-order mark
  if (!grepl("^Plot Name:", lines[1])) stop(basename(path), " is not an AfterMath text export")
  split = function(x) strsplit(x, ",", fixed = TRUE)[[1]]
  name_row = split(lines[grep("^Trace Name:", lines)[1]])
  unit_row = split(lines[grep("^Index,", lines)[1]])
  starts = which(name_row %in% c("Trace Name:", "Result Name:"))
  first = starts[!duplicated(name_row[starts + 1])]       # first copy of each trace
  # all rows into one character matrix, padded to the widest row, so each column is read once
  rows = strsplit(lines[(grep("^Index,", lines)[1] + 1):length(lines)], ",", fixed = TRUE)
  width = max(lengths(rows))
  body = do.call(rbind, lapply(rows, function(r) c(r, rep("", width - length(r)))))

  # the unit is the text in brackets; the micro sign can arrive as U+00B5, U+03BC or a
  # mis-decoded pair, so anything before the letters is read as "u"
  unit = function(h) {
    u = sub("[)].*$", "", sub("^.*[(]", "", h))
    u = sub("^[^A-Za-z]+", "u", u)
    if (u %in% names(pine_units)) pine_units[[u]] else stop("Unknown unit '", h, "' in ", basename(path))
  }
  block = function(i) {
    num = function(k) suppressWarnings(as.numeric(body[, k]))
    x = num(i + 2); y = num(i + 3)
    ok = !is.na(x) & !is.na(y)
    data.frame(segment = num(i + 1)[ok], x = x[ok] * unit(unit_row[i + 2]), y = y[ok] * unit(unit_row[i + 3]))
  }
  traces = stats::setNames(lapply(first, block), name_row[first + 1])
  pick = function(name) traces[[name]]
  k1i = pick("K1 Current vs Time"); k1e = pick("K1 Potential vs Time")
  if (is.null(k1i) || is.null(k1e)) stop(basename(path), ": no K1 current and potential traces found")
  data = data.frame(t_s = k1i$x, E_V = k1e$y[seq_len(nrow(k1i))], I_A = k1i$y)
  k2i = pick("K2 Current vs Time"); k2e = pick("K2 Potential vs Time")
  if (!is.null(k2i) && nrow(k2i) == nrow(data)) { data$I2_A = k2i$y; data$E2_V = k2e$y }
  data$segment = assign_segments(data$E_V)
  dt = stats::median(diff(data$t_s)); dE = stats::median(abs(diff(data$E_V)))
  params = list(technique = "Pine AfterMath export", instrument = "Pine WaveDriver",
                scan_rate_V_s = round(dE / dt, 4), high_E_V = max(data$E_V), low_E_V = min(data$E_V),
                ring_E_V = if (!is.null(data$E2_V)) stats::median(data$E2_V) else NA)
  list(file = basename(path), layout = "pine_aftermath", params = params, data = data,
       traces = names(traces))
}

# One reader for any supported file: CHI text, Pine AfterMath text
read_echem = function(path, name = basename(path), instrument = "auto") {
  fmt = instrument_format(path, name, instrument)
  if (fmt == "pssession") {
    x = read_pssession(path, name)
    if (x$kind == "eis") stop(name, " is an impedance measurement: open it in the Impedance (EIS) tab")
    return(x)
  }
  x = switch(fmt, biologic = read_biologic_mpt(path, name), mpr = read_biologic_mpr(path, name), gamry = read_gamry_dta(path, name),
             pine = read_pine_txt(path), chi = read_chi_txt(path, name), abf = read_abf(path, name), text = read_generic_text(path, name))
  if (is.null(x) || identical(x$kind, "eis")) stop(name, " is an impedance measurement: open it in the Impedance (EIS) tab")
  x$file = name
  x
}
