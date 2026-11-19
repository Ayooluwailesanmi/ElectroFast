# Reader for data files written by Lawson Labs EMF16_grf.
# Layout follows section F of the EMF16_grf manual: a header row starting with
# "Time / s", one column per channel, values in volts, an optional comment typed
# during the run, then a footer of run settings. Real EMF16 files are
# tab-separated and label channels with spaces ("1 - 1"); the manual's printed
# sample is space-separated ("1-7"), so both are handled.

read_emf16_grf = function(path) {
  lines = sub("\r$", "", readLines(path, warn = FALSE))

  hdr_i = grep("^\\s*Time\\s*/\\s*s", lines)[1]
  if (is.na(hdr_i)) stop("No 'Time / s' header found. Is this an EMF16_grf file?")
  tabbed = grepl("\t", lines[hdr_i])
  split_line = function(x) if (tabbed) strsplit(x, "\t") else strsplit(trimws(x), "[\t ]+")

  # Channel labels are "n-m" (spaces allowed around the dash); anything else on
  # the header line is the comment column title, not a channel.
  hdr = split_line(sub("^\\s*Time\\s*/\\s*s\\s*", "", lines[hdr_i]))[[1]]
  hdr = gsub("\\s", "", hdr)
  channels = hdr[grepl("^[0-9]+-[0-9]+$", hdr)]
  n_ch = length(channels)
  if (n_ch == 0) stop("No channel labels found in the header")

  body = lines[(hdr_i + 1):length(lines)]
  body = body[nzchar(trimws(body))]
  is_data = grepl("^\\s*[0-9]+\\.?[0-9]*\\s", body)
  footer = trimws(body[!is_data])

  # Fill a character matrix in one step (rows = data lines, columns = fields);
  # per-row loops made a 17,000-row file take several seconds.
  tok = split_line(body[is_data])
  lens = lengths(tok)
  m = matrix("", nrow = length(tok), ncol = max(lens, n_ch + 1))
  m[cbind(rep(seq_along(tok), lens), sequence(lens))] = trimws(unlist(tok, use.names = FALSE))
  wide = data.frame(time_s = as.numeric(m[, 1]))
  for (j in seq_len(n_ch)) wide[[channels[j]]] = as.numeric(m[, j + 1])
  # A comment typed during the run sits after the channel values on that row
  if (ncol(m) > n_ch + 1) {
    extra = trimws(gsub("\\s+", " ", do.call(paste, as.data.frame(m[, -(1:(n_ch + 1)), drop = FALSE]))))
    wide$comment = ifelse(nzchar(extra), extra, NA_character_)
  } else {
    wide$comment = NA_character_
  }

  # Long format is what the plateau and calibration steps use; mV because ISE
  # slopes are reported in mV per decade.
  long = data.frame(
    time_s = rep(wide$time_s, n_ch),
    channel = rep(channels, each = nrow(wide)),
    emf_mV = unlist(wide[channels], use.names = FALSE) * 1000
  )

  kv = regmatches(footer, regexec("^([^:]+):\\s*(.*)$", footer))
  kv = kv[lengths(kv) == 3]
  settings = setNames(trimws(vapply(kv, `[`, "", 3)), trimws(vapply(kv, `[`, "", 2)))

  # The logging interval comes from the footer; time is printed to 0.1 s, so at
  # fast intervals consecutive rows share a time stamp. Logging switched off and
  # on leaves a jump in time, flagged so later steps do not average across it.
  interval = suppressWarnings(as.numeric(settings["Applied Interval"]))
  if (is.na(interval)) interval = stats::median(diff(wide$time_s)[diff(wide$time_s) > 0])
  dt = diff(wide$time_s)
  gaps = which(dt > max(3 * interval, 0.3))

  if (any(abs(long$emf_mV) > 5000, na.rm = TRUE))
    warning("Values above 5 V after conversion. The file may already be in mV.")

  list(
    file = basename(path), channels = channels, wide = wide, long = long,
    comments = wide[!is.na(wide$comment), c("time_s", "comment")],
    settings = settings, interval_s = interval,
    gaps = data.frame(after_s = wide$time_s[gaps], resumes_s = wide$time_s[gaps + 1])
  )
}
