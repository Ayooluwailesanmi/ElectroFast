# Spike detection on simulated traces with a known answer: 25 triangular spikes of known height and charge on a
# drifting baseline with Gaussian noise. Run from the project folder: Rscript tests/check_spikes.R
source("R/spikes.R")
chk = function(ok, what) { cat(if (isTRUE(ok)) "ok  " else "FAIL", what, "\n"); if (!isTRUE(ok)) quit(status = 1) }
set.seed(1)
dt = 1e-3; t = seq(0, 60, by = dt); n = length(t)
base = 2e-9 + 1e-10 * t / 60 + 2e-11 * sin(2 * pi * t / 20)        # 2 nA, slow drift
noise = 1e-12                                                        # 1 pA
I = base + rnorm(n, 0, noise)
times = sort(sample(seq(2, 58, by = 0.5), 25)); amp = runif(25, 15e-12, 40e-12); width = 0.010   # 10 ms triangles
for (i in seq_along(times)) { k = abs(t - times[i]) <= width / 2; I[k] = I[k] + amp[i] * (1 - abs(t[k] - times[i]) / (width / 2)) }
true_q = amp * width / 2                                             # triangle area

r = detect_spikes(t, I, base_s = 1, k = 5, polarity = "positive")
e = r$events; cat("found", nrow(e), "of 25 events; noise", signif(r$noise_A * 1e12, 3), "pA\n")
chk(nrow(e) == 25, "every spike found, none extra")
m = sapply(times, function(x) which.min(abs(e$t_peak_s - x)))
chk(all(abs(e$t_peak_s[m] - times) <= 2 * dt), "peak times within 2 ms")
chk(abs(r$noise_A / noise - 1) < 0.1, "noise estimated within 10%")
ea = e$amplitude_A[m] / amp; eq = e$charge_C[m] / true_q
cat(sprintf("amplitude ratio %.3f to %.3f, charge ratio %.3f to %.3f, duration at half height %.1f ms (true 5)\n",
            min(ea), max(ea), min(eq), max(eq), 1e3 * stats::median(e$duration_s)))
chk(all(abs(ea - 1) < 0.15), "amplitudes within 15% (noise of 1 pA on 15 to 40 pA)")
chk(abs(stats::median(eq) - 1) < 0.1, "median charge within 10%")

# negative spikes and the empty case
r2 = detect_spikes(t, 2 * base - I, base_s = 1, k = 5, polarity = "negative"); chk(nrow(r2$events) == 25, "negative spikes with polarity = negative")
r3 = detect_spikes(t, base + rnorm(n, 0, noise), base_s = 1, k = 5); chk(nrow(r3$events) == 0, "pure noise: no events at k = 5")

# through the app: the trace saved as a CHI data-only i-t export (time, current, charge)
library(shiny)
f = file.path(tempdir(), "simulated impacts i-t.txt")
writeLines(sprintf("%.3f\t%.6e\t%.6e", t, I, cumsum(I) * dt), f)
app = shinyAppDir(".")
testServer(app, {
  session$setInputs(ca_files = data.frame(name = basename(f), size = file.size(f), type = "text/plain", datapath = f), ca_mode = "spikes",
                    ca_sfrom = 0, ca_sto = 1e9, ca_base = 1, ca_k = 5, ca_pol = "positive", ca_sep = 0, `ca_el-shape` = "preset", `ca_el-preset` = "2")
  s = ca_res(); print(s[, c("file", "events", "frequency_Hz", "noise_pA", "median_amplitude_pA", "median_charge_fC")], row.names = FALSE)
  chk(s$events == 25 && nrow(ca_events()) == 25, "app: 25 events from the i-t file")
  invisible(output$ca_plot); invisible(output$ca_events)
})
cat("spike checks passed\n")
