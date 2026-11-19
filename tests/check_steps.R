# Single-entity steps on a simulated staircase with known steps, drift, white noise and spikes; the trace thinning used
# for plots; the current unit. Run from the project folder: Rscript tests/check_steps.R
source("R/steps.R")
chk = function(ok, what) { cat(if (isTRUE(ok)) "ok  " else "FAIL", what, "\n"); if (!isTRUE(ok)) quit(status = 1) }
set.seed(3); fs = 10000; t = seq(0, 120, by = 1 / fs)
true_t = c(10.3, 25.7, 41.2, 58.9, 77.4, 96.1); true_h = c(250, 40, 120, 25, 80, -60) * 1e-12
I = 860e-12 + 0.5e-12 * t                                  # slow drift, 0.5 pA/s
for (i in seq_along(true_t)) I = I + true_h[i] * (t >= true_t[i])
I = I + rnorm(length(t), 0, 4e-12)                         # white noise, 4 pA
for (ts in runif(150, 0, 120)) { k = which(t >= ts & t < ts + 0.002); I[k] = I[k] + runif(1, 20, 300) * 1e-12 }   # 150 spikes of 2 ms
r = detect_steps(t, I, bin_s = 0.01, win_s = 0.5, k = 6)
mt = sapply(true_t, function(x) { d = abs(r$steps$t_step_s - x); if (length(d) && min(d) < 0.02) which.min(d) else NA })
chk(all(!is.na(mt)), "all 6 known steps found, each within 0.02 s")
chk(nrow(r$steps) == 6, "no false steps from 150 spikes, drift and noise")
chk(all(abs(r$steps$height_A[mt] - true_h) / abs(true_h) < 0.03), "step heights within 3%")
chk(nrow(detect_steps(t, I, k = 6, direction = "down")$steps) == 1, "direction down: only the one negative step")
flat = detect_steps(t, 860e-12 + rnorm(length(t), 0, 4e-12), k = 6)
chk(nrow(flat$steps) == 0, "no steps in noise alone")
s = step_summary("sim", t, r, 6); chk(s$steps == 6 && abs(s$threshold_pA - 6 * r$spread_A * 1e12) < 1e-9, "summary: count and threshold")
th = thin_trace(t, I, 2000)
chk(nrow(th) <= 2 * 2001 + 1 && max(th$y) == max(I) && min(th$y) == min(I), "thinned trace keeps the highest and lowest points")
chk(identical(current_unit(I)$label, "nA") && identical(current_unit(5e-11)$label, "pA") && identical(current_unit(2e-6)$label, "µA"), "current unit")
cat("step checks passed\n")
