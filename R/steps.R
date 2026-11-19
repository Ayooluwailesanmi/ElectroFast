# Single-entity steps in a current-time trace: lasting changes of the current, as when a particle lands on an
# ultramicroelectrode and stays active (a "staircase"). Spikes, which come and go, are found by detect_spikes().
# The current is averaged in short bins. At each bin edge the step statistic is the median of the bins after the edge
# minus the median of the bins before it (win_s on each side); medians keep spikes from looking like steps. The spread
# of this statistic over the whole trace (1.4826 x median absolute deviation) sets the noise, so slow noise and drift
# raise the threshold instead of giving false steps. A step is a local extreme of the statistic beyond k x that spread,
# at least one window from the next; its height is the statistic there.

detect_steps = function(t, I, bin_s = 0.01, win_s = 0.5, k = 6, direction = c("both", "up", "down")) {
  direction = match.arg(direction)
  o = order(t); t = t[o]; I = I[o]; n = length(t)
  dt = stats::median(diff(t[seq_len(min(n, 10000))]))
  b = max(1, round(bin_s / dt)); nb = floor(n / b)
  empty = data.frame(t_step_s = numeric(), height_A = numeric(), level_before_A = numeric(), level_after_A = numeric())
  W = max(3, round(win_s / (b * dt)))
  if (nb < 2 * W + 3) return(list(steps = empty, spread_A = NA_real_, bin_s = b * dt, window_bins = W, t_bin = numeric(), stat = numeric()))
  m = .colMeans(I[seq_len(nb * b)], b, nb)
  tb = t[seq(1, nb * b, by = b)] + (b - 1) * dt / 2
  # difference of medians either side of the edge after bin j
  j = W:(nb - W)
  before = vapply(j, function(i) stats::median(m[(i - W + 1):i]), 0)
  after = vapply(j, function(i) stats::median(m[(i + 1):(i + W)]), 0)
  D = after - before
  spread = 1.4826 * stats::median(abs(D - stats::median(D)))
  s = switch(direction, both = abs(D), up = D, down = -D)
  cand = which(s > k * spread)
  # local extremes, at least one window apart: the largest first, then anything within W bins of it is dropped
  keep = integer(0)
  for (i in cand[order(-s[cand])]) if (!length(keep) || all(abs(i - keep) >= W)) keep = c(keep, i)
  keep = sort(keep)
  # the median statistic is flat near an edge, so each edge is placed where the bins within one window either side
  # split best into two levels (least squares on the bin means)
  edge_bin = vapply(j[keep], function(e) {
    lo = max(1, e - W); hi = min(nb, e + W); seg = m[lo:hi]; L = length(seg)
    sse = vapply(1:(L - 1), function(c) sum((seg[1:c] - mean(seg[1:c]))^2) + sum((seg[(c + 1):L] - mean(seg[(c + 1):L]))^2), 0)
    lo - 1 + which.min(sse) }, 0)
  edge_t = (tb[edge_bin] + tb[pmin(edge_bin + 1, nb)]) / 2
  # a step within one window of a larger step is the settling of that step (an overshoot after a landing, for
  # example), not a separate event, so it is dropped
  if (length(keep) > 1) {
    h = abs(D[keep]); ok = vapply(seq_along(keep), function(i) !any(abs(edge_t - edge_t[i]) < W * b * dt & h > h[i]), TRUE)
    keep = keep[ok]; edge_t = edge_t[ok]
  }
  steps = if (length(keep)) data.frame(t_step_s = edge_t, height_A = D[keep], level_before_A = before[keep], level_after_A = after[keep]) else empty
  list(steps = steps, spread_A = spread, bin_s = b * dt, window_bins = W, t_bin = (tb[j] + tb[j + 1]) / 2, stat = D)
}

# one row per trace: how many steps, how often, their sizes, and the threshold used (k x spread)
step_summary = function(file, t, res, k) {
  s = res$steps; dur = diff(range(t))
  data.frame(file = file, steps = nrow(s), frequency_Hz = nrow(s) / dur,
             median_height_pA = if (nrow(s)) stats::median(s$height_A) * 1e12 else NA_real_,
             largest_pA = if (nrow(s)) s$height_A[which.max(abs(s$height_A))] * 1e12 else NA_real_,
             threshold_pA = k * res$spread_A * 1e12, spread_pA = res$spread_A * 1e12, bin_ms = res$bin_s * 1e3, window_s = res$window_bins * res$bin_s)
}

# For drawing long traces: in each of n_bins bins the lowest and the highest point, in time order, so spikes and steps
# stay visible while a 30-million-point recording is drawn from about 2 x n_bins points
thin_trace = function(t, y, n_bins = 4000) {
  n = length(t); if (n <= 2 * n_bins) return(data.frame(t = t, y = y))
  b = ceiling(n / n_bins); nb = floor(n / b); M = matrix(y[seq_len(nb * b)], nrow = b)
  rest = if (nb * b < n) { r = (nb * b + 1):n; c(r[which.min(y[r])], r[which.max(y[r])]) }   # the part after the last full bin
  idx = sort(unique(c((seq_len(nb) - 1) * b + apply(M, 2, which.min), (seq_len(nb) - 1) * b + apply(M, 2, which.max), rest, n)))
  data.frame(t = t[idx], y = y[idx])
}

# a current unit that suits the size of the current: pA, nA, µA or mA, as a scale factor and a label
current_unit = function(I) {
  m = max(abs(I), na.rm = TRUE)
  if (!is.finite(m) || m >= 1e-3) list(scale = 1e3, label = "mA") else if (m >= 1e-6) list(scale = 1e6, label = "µA") else
    if (m >= 1e-9) list(scale = 1e9, label = "nA") else list(scale = 1e12, label = "pA")
}
