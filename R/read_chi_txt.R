# Reader for CH Instruments text files.
# Two layouts are seen in CHI exports:
#   1. "Convert to text" export: header with date, technique, instrument and
#      parameters, CHI's own results, then a column line such as "Potential/V, Current/A"
#      or "Time/sec, Current/A", and comma-separated rows. Every technique is written this
#      way: CV and LSV, DPV, OCP, i-t, chronoamperometry, multi-potential steps,
#      chronopotentiometry, and bipotentiostat runs with a second current (i2/A).
#   2. Headerless, tab-separated, four columns: potential (V), current (A),
#      cumulative charge (C) and time (s). No parameters are stored, so the scan
#      rate is recovered from potential and time.

# CHI column labels and the names used everywhere else in the app
# (square wave: the difference current is the response, forward and reverse are kept)
chi_columns = c("Potential/V" = "E_V", "Current/A" = "I_A", "i1/A" = "I_A", "i2/A" = "I2_A",
                "Time/sec" = "t_s", "Charge/C" = "Q_C", "Diff(i/A)" = "I_A", "For(i/A)" = "I_fwd_A", "Rev(i/A)" = "I_rev_A")

# "Data only" text exports carry numbers alone. The columns CHI writes for each technique
# (checked against exports with headers, October 2026):
#   CV, LSV: potential, current, charge, time      SWV: potential, difference, forward, reverse current
#   DPV, NPV: potential, current                   i-t, CA: time, current, charge
#   OCP, CP: time, potential                       A.C. impedance: frequency, Z', Z'', |Z|, phase
# The technique is taken from the file name when it says (SWV, DPV, OCP, EIS, ...), otherwise
# from the number of columns and how they relate to each other.
chi_dataonly = function(m, name) {
  nm = toupper(name); has = function(p) grepl(paste0("(^|[^A-Z])(", p, ")([^A-Z]|$)"), nm)
  hint = if (has("EIS|IMPEDANCE|IMP")) "eis" else if (has("SWV|SWSV")) "swv" else if (has("DPV|DPSV|NPV")) "dpv" else if (has("OCP|OCPT")) "ocp" else
    if (has("CP|CHRONOPOT")) "cp" else if (has("CA|CC|I-T|IT|CHRONOAMP")) "ca" else if (has("CV|LSV")) "cv" else ""
  v = function(j) as.numeric(m[, j]); n = ncol(m)
  tech = if (n == 5) "eis" else if (n == 4) {
    if (hint == "swv" || stats::median(abs(v(2) - (v(3) - v(4)))) < 1e-3 * stats::median(abs(v(2))) + 1e-15) "swv" else "cv"
  } else if (n == 3) "ca" else if (n == 2) {
    if (hint %in% c("ocp", "cp", "dpv", "cv")) hint else if (max(abs(v(2))) > 0.01) "ocp" else "dpv"
  } else stop("Unrecognised data-only layout in ", name, ": ", n, " columns")
  how = if (nzchar(hint) && (hint == tech || (hint == "cp" && tech == "ocp"))) "file name" else "column layout"
  d = switch(tech,
    eis = data.frame(freq_Hz = v(1), Z = complex(real = v(2), imaginary = v(3))),
    swv = data.frame(E_V = v(1), I_A = v(2), I_fwd_A = v(3), I_rev_A = v(4)),
    cv  = data.frame(E_V = v(1), I_A = v(2), Q_C = v(3), t_s = v(4)),
    ca  = data.frame(t_s = v(1), I_A = v(2), Q_C = v(3)),
    dpv = data.frame(E_V = v(1), I_A = v(2)),
    ocp = , cp = data.frame(t_s = v(1), E_V = v(2)))
  label = c(eis = "A.C. Impedance", swv = "Square Wave Voltammetry", cv = "Cyclic or Linear Sweep Voltammetry", ca = "Chronoamperometry (i-t)",
            dpv = "Pulse Voltammetry (DPV or NPV)", ocp = "Open Circuit Potential", cp = "Chronopotentiometry")[[tech]]
  list(tech = tech, data = d, technique = paste0(label, " (data-only export, technique from the ", how, ")"))
}

read_chi_txt = function(path, name = basename(path)) {
  lines = sub("\r$", "", readLines(path, warn = FALSE))
  num = "^\\s*[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?\\s*$"

  hdr_i = grep("^(Potential/V|Time/sec)", lines)[1]
  if (!is.na(hdr_i)) {
    params = parse_chi_header(lines[1:(hdr_i - 1)])
    cols = trimws(strsplit(lines[hdr_i], ",")[[1]])
    body = lines[(hdr_i + 1):length(lines)]
    body = body[nzchar(trimws(body))]
    m = do.call(rbind, strsplit(body, ","))
    data = as.data.frame(apply(m, 2, as.numeric))
    names(data) = ifelse(cols %in% names(chi_columns), chi_columns[cols], make.names(cols))
    layout = "chi_export"
  } else {
    body = lines[nzchar(trimws(lines))]
    sep = if (grepl("\t", body[1])) "\t" else ","
    first = trimws(strsplit(body[1], sep)[[1]])
    if (length(first) < 2 || !all(grepl(num, first)))
      stop("Unrecognised CHI text layout in ", name)
    m = do.call(rbind, lapply(strsplit(body, sep), trimws))
    x = chi_dataonly(m, name)
    if (x$tech == "eis")
      return(list(file = name, kind = "eis", layout = "headerless", params = list(technique = x$technique), data = x$data, note = ""))
    data = x$data
    params = list(technique = x$technique, instrument = NA_character_, scan_rate_V_s = NA_real_)
    if (x$tech == "cv") {
      dt = stats::median(diff(data$t_s)); dE = stats::median(abs(diff(data$E_V)))
      params = c(params[1:2], list(scan_rate_V_s = round(dE / dt, 6), sample_interval_V = dE, high_E_V = max(data$E_V), low_E_V = min(data$E_V)))
    }
    layout = "headerless_tab"
  }

  # Sweeps only exist when potential is the swept variable (first column). In a
  # time-based file (OCP, chronopotentiometry) the potential drifts and has no sweeps.
  data$segment = if (names(data)[1] == "E_V") assign_segments(data$E_V) else 1L
  if (names(data)[1] == "E_V" && !is.null(params$segments) && !is.na(params$segments) &&
      max(data$segment) != params$segments)
    warning(basename(path), ": found ", max(data$segment), " segments, header says ", params$segments)

  list(file = name, layout = layout, params = params, data = data,
       chi_peaks = attr(params, "peaks"))
}

# Segments are the sweeps between turning points, where the direction of the
# potential step changes sign.
assign_segments = function(E) {
  d = sign(diff(E))
  d[d == 0] = NA
  d = zoo_fill(d)
  seg = c(1, 1 + cumsum(diff(d) != 0))
  c(seg, seg[length(seg)])
}

zoo_fill = function(x) {
  for (i in seq_along(x)[-1]) if (is.na(x[i])) x[i] = x[i - 1]
  if (is.na(x[1])) x[1] = x[which(!is.na(x))[1]]
  x
}

parse_chi_header = function(h) {
  get = function(key) {
    ln = grep(paste0("^", key), h, value = TRUE)[1]
    if (is.na(ln)) return(NA_character_)
    trimws(sub("^[^=:]*[=:]", "", ln))
  }
  getn = function(key) suppressWarnings(as.numeric(get(key)))
  p = list(
    date = h[1], technique = h[2],
    instrument = get("Instrument Model"),
    init_E_V = getn("Init E \\(V\\)"), high_E_V = getn("High E \\(V\\)"),
    low_E_V = getn("Low E \\(V\\)"), scan_rate_V_s = getn("Scan Rate \\(V/s\\)"),
    segments = getn("Segment ="), sample_interval_V = getn("Sample Interval \\(V\\)"),
    quiet_time_s = getn("Quiet Time"), sensitivity_A_V = getn("Sensitivity"),
    channel = suppressWarnings(as.integer(sub("^Channel ([0-9]+):.*", "\\1",
                                              grep("^Channel [0-9]+:", h, value = TRUE)[1])))
  )
  # Every "Name = value" setting, so each technique can find its own (E1 (V), T1 (s),
  # Cycle, Pulse Width (sec), Cathodic Current (A), 2nd E (V) and so on). Numbers are
  # stored as numbers; "30." and similar CHI spellings are read too.
  kv = grep("^[^=]+=[^=]+$", h, value = TRUE)
  key = trimws(sub("=.*$", "", kv))
  val = trimws(sub("^[^=]*=", "", kv))
  nv = suppressWarnings(as.numeric(val))
  p$settings = stats::setNames(lapply(seq_along(val), function(i) if (is.na(nv[i])) val[i] else nv[i]), key)

  # CHI writes its own peak results after the channel line; keep them so the
  # pipeline's peaks can be checked against the instrument software.
  start = grep("^Channel [0-9]+:", h)[1]
  peaks = NULL
  if (!is.na(start) && start < length(h)) {
    seg = NA
    rec = list()
    cur = list()
    for (ln in h[(start + 1):length(h)]) {
      if (grepl("^Segment [0-9]+:", ln)) seg = as.integer(sub("^Segment ([0-9]+):.*", "\\1", ln))
      else if (grepl("^Ep = ", ln)) cur = list(segment = seg, Ep_V = as.numeric(sub("^Ep = ([^V]+)V.*", "\\1", ln)))
      else if (grepl("^ip = ", ln)) cur$ip_A = as.numeric(sub("^ip = ([^A]+)A.*", "\\1", ln))
      else if (grepl("^Ah = ", ln)) { cur$Ah_C = as.numeric(sub("^Ah = ([^C]+)C.*", "\\1", ln)); rec[[length(rec) + 1]] = as.data.frame(cur) }
    }
    if (length(rec)) peaks = do.call(rbind, rec)
  }
  attr(p, "peaks") = peaks
  p
}

# One setting from the header, or a default when the file does not have it
chi_setting = function(x, key, default = NA) {
  v = x$params$settings[[key]]
  if (is.null(v)) default else v
}
