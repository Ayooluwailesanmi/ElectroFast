# Sensor performance functions on tables with known answers. Run from the project folder: Rscript tests/check_performance.R
source("R/methods.R"); source("R/steps.R"); source("R/performance.R")
chk = function(ok, what) { cat(if (isTRUE(ok)) "ok  " else "FAIL", what, "\n"); if (!isTRUE(ok)) quit(status = 1) }
set.seed(2)

# sensitivity: y = 2 + 3 c with noise SD 0.1, blanks with SD 0.1
cc = rep(c(0, 1, 2, 4, 8, 16), each = 3); y = 2 + 3 * cc + rnorm(length(cc), 0, 0.1)
s = perf_sensitivity(data.frame(`Concentration (mM)` = cc, `Signal (uA)` = y, check.names = FALSE), area_cm2 = 0.0314)
chk(abs(s$fit$slope - 3) < 0.02 && s$fit$r2 > 0.999, "sensitivity: slope 3 recovered, units in the headers ignored")
sb = sd(y[cc == 0]); chk(abs(s$fit$lod - 3.3 * sb / s$fit$slope) < 1e-12 && abs(s$fit$loq - 10 * sb / s$fit$slope) < 1e-12, "sensitivity: LOD 3.3 s/slope and LOQ 10 s/slope from the blanks")
chk(abs(s$fit$slope_per_area - s$fit$slope / 0.0314) < 1e-9, "sensitivity: slope per unit area")
s2 = perf_sensitivity(data.frame(concentration = cc[cc > 0], signal = y[cc > 0]), from = 1, to = 8)
chk(s2$fit$to == 8 && grepl("residual", s2$fit$lod_basis), "sensitivity: range limits and the residual SD without blanks")

# selectivity: one interferent changes the signal by +10%, the other by -2%
sel = perf_selectivity(data.frame(interferent = c("ascorbate", "ascorbate", "urea"), signal_analyte = c(10, 10, 10), signal_mixture = c(11, 11, 9.8)), 5)
chk(isTRUE(all.equal(sel$change_pct, c(10, -2))) && identical(sel$within_tolerance, c(FALSE, TRUE)) && sel$replicates[1] == 2, "selectivity: % change, replicates averaged, tolerance")

# specificity: cross-reactivity relative to the target
sp = perf_specificity(data.frame(species = c("glucose", "fructose", "galactose"), signal = c(50, 5, 1), target = c("yes", "no", "no")))
chk(isTRUE(all.equal(sp$cross_reactivity_pct, c(100, 10, 2))) && isTRUE(all.equal(sp$selectivity_ratio, c(1, 10, 50))), "specificity: cross-reactivity 100, 10 and 2%, ratios 1, 10 and 50")
sp0 = perf_specificity(data.frame(species = c("a", "b"), signal = c(4, 0), target = c("yes", "no"))); chk(is.na(sp0$selectivity_ratio[2]), "specificity: no ratio for a species with no signal")

# stability: 100, 95, 85% at days 0, 7, 14: falls below 90% at day 10.5
st = perf_stability(data.frame(time = c(0, 7, 14), signal = c(20, 19, 17)))
chk(isTRUE(all.equal(st$points$retention_pct, c(100, 95, 85))) && abs(st$summary$first_time_below_limit - 10.5) < 1e-9, "stability: retention and the day it falls below 90% (10.5)")

# speed: first-order rise with time constant 2 s after an addition at 10 s: t90 = 2 ln 10 = 4.61 s
t = seq(0, 60, by = 0.01); yy = ifelse(t < 10, 1, 1 + 5 * (1 - exp(-(t - 10) / 2))) + rnorm(length(t), 0, 0.01)
v = perf_speed(t, yy, t_add = 10)
chk(abs(v$t90_s - 2 * log(10)) < 0.1 && abs(v$t95_s - 2 * log(20)) < 0.1, sprintf("speed: t90 %.2f s (true 4.61), t95 %.2f s (true 5.99)", v$t90_s, v$t95_s))
vd = perf_speed(t, 2 - (yy - 1), t_add = 10); chk(abs(vd$t90_s - 2 * log(10)) < 0.1, "speed: a falling signal works the same way")

# precision: day SD 2, electrode SD 1, residual SD 0.5, 8 days x 6 electrodes x 3 replicates
g = expand.grid(rep = 1:3, electrode = paste0("E", 1:6), day = paste0("D", 1:8))
de = rnorm(8, 0, 2); ee = rnorm(6, 0, 1)
g$value = 100 + de[as.integer(factor(g$day))] + ee[as.integer(factor(g$electrode))] + rnorm(nrow(g), 0, 0.5)
p = perf_precision(g)
sdv = setNames(p$components$SD, p$components$component)
chk(abs(sdv[["Repeatability (within-run SD)"]] - 0.5) < 0.1, sprintf("precision: repeatability %.2f (true 0.5)", sdv[["Repeatability (within-run SD)"]]))
chk(isTRUE(all.equal(sdv[["Intermediate precision (within the laboratory)"]]^2,
                     sum(sdv[c("Repeatability (within-run SD)", "Between-day SD", "Between-electrode SD")]^2))), "precision: intermediate precision combines the within-laboratory components")
chk(!any(grepl("Reproducibility", names(sdv))), "precision: no reproducibility without a laboratory column")
g$laboratory = rep(c("Lab A", "Lab B"), length.out = nrow(g)); p2 = perf_precision(g)
chk(any(grepl("Reproducibility", p2$components$component)) && any(grepl("Only 2 laboratory", p2$notes)), "precision: reproducibility with two laboratories, flagged as indicative")
chk(length(perf_precision(data.frame(value = c(1, 2, 3)))$notes) == 1, "precision: values alone give the repeatability only, with a note")
e = tryCatch(perf_sensitivity(data.frame(conc = 1:3)), error = function(e) conditionMessage(e)); chk(grepl("template", e), "missing columns: a clear message")
cat("sensor performance checks passed\n")
