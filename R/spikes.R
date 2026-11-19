# Single-entity events in a current-time trace (for example particle impacts or collisions, or any other transient
# spikes): the baseline is a running median over a window much longer than an event, the noise is the robust SD of
# the baseline-corrected current (1.4826 x median absolute deviation), and an event is a run of points beyond
# k x noise. For each event: time of the peak, amplitude above the baseline, duration at half height, and charge
# (the corrected current integrated over the event, out to where it falls back inside 1 x noise).
# Current steps (a lasting change, as when particles stick and stay active) are not events here: the running median
# follows a step, so only its edge would show.

detect_spikes = function(t, I, base_s = 1, k = 5, polarity = c("positive", "negative", "both"), min_sep_s = 0) {
  polarity = match.arg(polarity)
  o = order(t); t = t[o]; I = I[o]; n = length(t)
  dt = stats::median(diff(t)); w = max(3, 2 * floor(base_s / dt / 2) + 1); w = min(w, if (n %% 2) n else n - 1)
  # very long traces (patch-clamp recordings hold millions of points): the running median is taken on every j-th point
  # and interpolated back, about 400 points per window; shorter traces use every point as before
  base = if (n > 2e5 && w > 801) { j = floor(w / 401); ks = seq(1, n, by = j); wk = 2 * floor(w / j / 2) + 1
    stats::approx(t[ks], as.numeric(stats::runmed(I[ks], min(wk, length(ks) - (1 - length(ks) %% 2)), endrule = "median")), xout = t, rule = 2)$y
  } else as.numeric(stats::runmed(I, w, endrule = "median"))
  r = I - base; noise = 1.4826 * stats::median(abs(r - stats::median(r)))
  s = switch(polarity, positive = r, negative = -r, both = abs(r))
  hit = s > k * noise
  empty = data.frame(t_peak_s = numeric(), amplitude_A = numeric(), duration_s = numeric(), charge_C = numeric(), t_start_s = numeric(), t_end_s = numeric())
  if (!any(hit)) return(list(events = empty, noise_A = noise, baseline = base, window_points = w, dt_s = dt))
  rl = rle(hit); ends = cumsum(rl$lengths); starts = ends - rl$lengths + 1
  runs = cbind(starts[rl$values], ends[rl$values])
  # runs closer than min_sep_s are one event
  if (nrow(runs) > 1 && min_sep_s > 0) {
    keep = list(runs[1, ]); for (i in 2:nrow(runs)) { last = keep[[length(keep)]]
      if (t[runs[i, 1]] - t[last[2]] < min_sep_s) keep[[length(keep)]] = c(last[1], runs[i, 2]) else keep[[length(keep) + 1]] = runs[i, ] }
    runs = do.call(rbind, keep)
  }
  # each run is widened to where the current falls back inside 1 x noise; runs that widen into the same stretch are
  # one event (without this merge, one event was listed once per run inside it)
  wide = t(apply(runs, 1, function(r) { a = r[1]; b = r[2]
    while (a > 1 && s[a - 1] > noise) a = a - 1
    while (b < n && s[b + 1] > noise) b = b + 1
    c(a, b) }))
  if (!is.matrix(wide)) wide = matrix(wide, ncol = 2)
  m = list(wide[1, ]); for (i in seq_len(nrow(wide))[-1]) { last = m[[length(m)]]
    if (wide[i, 1] <= last[2] + 1) m[[length(m)]] = c(last[1], max(last[2], wide[i, 2])) else m[[length(m) + 1]] = wide[i, ] }
  wide = do.call(rbind, m)
  ev = do.call(rbind, lapply(seq_len(nrow(wide)), function(i) {
    a = wide[i, 1]; b = wide[i, 2]
    j = a - 1 + which.max(s[a:b]); h = s[j] / 2
    l = j; while (l > a && s[l - 1] >= h) l = l - 1
    u = j; while (u < b && s[u + 1] >= h) u = u + 1
    q = if (b > a) sum(diff(t[a:b]) * (r[a:b][-1] + r[a:b][-(b - a + 1)]) / 2) else r[a] * dt
    data.frame(t_peak_s = t[j], amplitude_A = r[j], duration_s = t[u] - t[l] + dt, charge_C = q, t_start_s = t[a], t_end_s = t[b])
  }))
  list(events = ev, noise_A = noise, baseline = base, window_points = w, dt_s = dt)
}

spike_summary = function(file, t, res, k) {
  e = res$events; span = diff(range(t))
  data.frame(file = file, duration_s = span, noise_pA = 1e12 * res$noise_A, threshold_pA = 1e12 * k * res$noise_A, events = nrow(e),
             frequency_Hz = nrow(e) / span,
             median_amplitude_pA = if (nrow(e)) 1e12 * stats::median(e$amplitude_A) else NA,
             median_duration_ms = if (nrow(e)) 1e3 * stats::median(e$duration_s) else NA,
             median_charge_fC = if (nrow(e)) 1e15 * stats::median(e$charge_C) else NA)
}
