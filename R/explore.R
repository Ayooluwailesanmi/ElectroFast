# Exploration of current-time traces: a quality check of each file, the window of a trace shown in the zoom, events
# counted over time, and the comparison of files with a reference (for example a blank). These describe the data and
# help to judge the automatic results. They do not change any result.

# Indexes of the points of a sorted time vector between a and b
trace_window = function(t, a, b) {
  i = findInterval(c(a, b), t); lo = max(1, i[1]); hi = min(length(t), max(lo, i[2]))
  lo:hi
}

# Bin medians of a current trace: n_bins equal bins of points (the median keeps spikes out); returns time and current
bin_medians = function(t, I, n_bins = 300) {
  n = length(t); b = max(1, floor(n / n_bins)); nb = floor(n / b)
  M = matrix(I[seq_len(nb * b)], nrow = b)
  list(t = t[(seq_len(nb) - 1) * b + ceiling(b / 2)], I = apply(M, 2, stats::median), b = b, nb = nb, M = M)
}

# One row per file: what was recorded and whether it looks sound before any analysis
#   noise       robust SD of the current about local medians (bins of about 1000 points), so steps, drift and spikes do not count
#   net change  median of the last 1% minus the median of the first 1% of the trace
#   gaps        steps in time longer than 1.5 x the usual sampling interval
#   at limit    share of points equal to the largest or the smallest value: a flat top or bottom suggests a clipped
#               (saturated) amplifier or a range set too small
file_check = function(x) {
  d = x$data; t = d$t_s; I = d$I_A; n = length(t)
  dts = diff(t); dt = stats::median(dts)
  bm = bin_medians(t, I, max(1, min(n %/% 1000, 3000)))
  cols = unique(round(seq(1, bm$nb, length.out = min(bm$nb, 300))))
  res = bm$M[, cols, drop = FALSE] - rep(bm$I[cols], each = nrow(bm$M))
  noise = 1.4826 * stats::median(abs(res - stats::median(res)))
  k = max(1, n %/% 100); first = stats::median(I[seq_len(k)]); last = stats::median(I[(n - k + 1):n])
  at_limit = max(mean(I == max(I)), mean(I == min(I)))
  p = x$params
  data.frame(file = x$file, points = n, duration_s = t[n] - t[1], sampling_Hz = round(1 / dt, 6),
             min_A = min(I), max_A = max(I), noise_A = noise, net_change_A = last - first, gaps = sum(dts > 1.5 * dt),
             at_limit_pct = 100 * at_limit, E_V = if (!is.null(p$E_V)) p$E_V else NA_real_,
             filter_Hz = if (!is.null(p$filter_Hz)) p$filter_Hz else NA_real_, check.names = FALSE)
}

# Plain-language notes for the file check: only what was found
file_check_notes = function(fc) {
  out = character(0)
  for (i in seq_len(nrow(fc))) {
    r = fc[i, ]
    if (r$gaps > 0) out = c(out, sprintf("%s: %d gap(s) in time, longer than 1.5 times the sampling interval.", r$file, r$gaps))
    if (r$at_limit_pct > 0.1) out = c(out, sprintf("%s: %.2f%% of points sit at the largest or smallest value. The current may have been clipped (amplifier range too small).", r$file, r$at_limit_pct))
    if (!is.na(r$filter_Hz) && r$filter_Hz > 0 && r$sampling_Hz < 2.5 * r$filter_Hz)
      out = c(out, sprintf("%s: sampled at %.3g Hz, less than 2.5 times the %.3g Hz filter.", r$file, r$sampling_Hz, r$filter_Hz))
  }
  out
}

# Events (spike or step times) counted in bins of bin_s over the trace
event_rate = function(times, t0, t1, bin_s = 10) {
  br = seq(t0, t0 + max(1, ceiling((t1 - t0) / bin_s - 1e-9)) * bin_s, by = bin_s)   # whole bins from t0, the last one ending at or just after t1
  h = graphics::hist(times[times >= t0 & times <= t1], breaks = br, plot = FALSE)
  data.frame(t_mid_s = h$mids, count = h$counts, rate_per_s = h$counts / bin_s)
}

# Several traces on a common time grid (bin means over the time all files share), each optionally shifted to start at
# zero, and each minus a reference file. Bin means keep the comparison of levels; spikes are averaged within a bin.
compare_traces = function(xs, ref = 1, n_bins = 1000, start_zero = FALSE) {
  t0 = max(vapply(xs, function(x) min(x$data$t_s), 0)); t1 = min(vapply(xs, function(x) max(x$data$t_s), 0))
  if (!(t1 > t0)) stop("the files share no time range")
  br = seq(t0, t1, length.out = n_bins + 1); mid = (br[-1] + br[-length(br)]) / 2
  m = lapply(xs, function(x) { k = trace_window(x$data$t_s, t0, t1); g = findInterval(x$data$t_s[k], br, rightmost.closed = TRUE)
    v = rep(NA_real_, n_bins); s = tapply(x$data$I_A[k], g, mean); v[as.integer(names(s))] = s; v[seq_len(n_bins)] })
  if (start_zero) m = lapply(m, function(v) v - v[which(!is.na(v))[1]])
  d = do.call(rbind, lapply(seq_along(xs), function(i) data.frame(file = xs[[i]]$file, t_s = mid, I_A = m[[i]], diff_A = m[[i]] - m[[ref]])))
  d
}
