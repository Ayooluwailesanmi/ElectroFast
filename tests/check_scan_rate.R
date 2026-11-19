# Scan-rate fits for both peaks (scan_rate_fits, scan_rate_reading) on peak currents with a known answer.
# Run from the project folder: Rscript tests/check_scan_rate.R
source("R/methods.R")
chk = function(ok, what) { cat(if (isTRUE(ok)) "ok  " else "FAIL", what, "\n"); if (!isTRUE(ok)) quit(status = 1) }
v = c(0.01, 0.025, 0.05, 0.1, 0.2, 0.5)

# diffusing species: ip proportional to v^1/2, the return peak a little smaller and of opposite sign
f = scan_rate_fits(v, 8e-5 * sqrt(v), -7e-5 * sqrt(v))
ll = f$loglog; sq = f$fits[f$fits$against == "√v", ]
chk(all(abs(ll$loglog_slope - 0.5) < 1e-9), "diffusion: log-log slope 0.5 for both peaks")
chk(isTRUE(all.equal(sq$slope, c(8e-5, 7e-5))) && all(abs(sq$intercept_A) < 1e-15) && all(sq$r2 > 1 - 1e-12), "diffusion: line against v^1/2 through the origin, slopes 8e-5 and 7e-5")
chk(all(f$fits$r2[f$fits$against == "v"] < 0.99), "diffusion: the line against v fits worse")
chk(all(grepl("closer to 0.5", scan_rate_reading(f))) && all(grepl("close to zero", scan_rate_reading(f))), "diffusion: reading says 0.5 and intercept zero")

# surface-confined species: ip proportional to v
g = scan_rate_fits(v, 3e-4 * v, -3e-4 * v)
chk(all(abs(g$loglog$loglog_slope - 1) < 1e-9) && all(g$fits$r2[g$fits$against == "v"] > 1 - 1e-12), "adsorbed: log-log slope 1, line against v exact")
chk(all(grepl("closer to 1", scan_rate_reading(g))), "adsorbed: reading says 1")

# a constant background added to a diffusing peak: the intercept against v^1/2 is no longer zero
set.seed(1); h = scan_rate_fits(v, 8e-5 * sqrt(v) + 2e-6 + rnorm(6, 0, 1e-8), -7e-5 * sqrt(v))
chk(grepl("may contribute", scan_rate_reading(h)[1]) && grepl("close to zero", scan_rate_reading(h)[2]), "a 2 uA background (25% of the smallest peak) is flagged for that peak only")
chk(is.null(scan_rate_fits(v[1:2], c(1, 2), c(-1, -2))) && grepl("At least 3", scan_rate_reading(NULL)), "fewer than 3 scan rates: no fit, with a message")
cat("scan-rate checks passed\n")
