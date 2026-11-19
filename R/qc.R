# Quality control across sessions: a history of results (one row per session and quantity), control limits from a
# reference period, and the Westgard multi-rule checks (Westgard, Barry, Hunt and Groth, Clin. Chem. 1981, 27, 493):
#   1-3s  one point beyond 3 SD                      2-2s  two in a row beyond 2 SD on the same side
#   R-4s  two in a row more than 4 SD apart (applied here to consecutive sessions)   4-1s  four in a row beyond 1 SD on the same side
#   10x   ten in a row on the same side of the mean
# 1-2s (one point beyond 2 SD) is a warning only. shinyapps.io keeps nothing between visits, so the history lives
# in a CSV the user downloads and loads again next time.

qc_columns = c("date", "quantity", "value", "unit", "electrode", "note")
qc_template = function() data.frame(date = character(), quantity = character(), value = numeric(), unit = character(),
                                     electrode = character(), note = character())

qc_read = function(path) {
  h = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(h) = tolower(trimws(names(h)))
  miss = setdiff(c("date", "quantity", "value"), names(h))
  if (length(miss)) stop("The history needs the columns date, quantity and value (missing: ", paste(miss, collapse = ", "), ")")
  for (k in setdiff(qc_columns, names(h))) h[[k]] = ""
  h$value = suppressWarnings(as.numeric(h$value))
  h[, qc_columns]
}

# Limits from the first n_ref points of each quantity (all points when n_ref is NA)
qc_limits = function(h, n_ref = NA) do.call(rbind, lapply(split(h, h$quantity), function(d) {
  d = d[order(d$date), ]; r = if (is.na(n_ref)) d$value else utils::head(d$value, n_ref)
  data.frame(quantity = d$quantity[1], unit = d$unit[1], n_ref = length(r), mean = mean(r), sd = if (length(r) > 1) stats::sd(r) else NA)
}))

qc_rules = function(z) {
  n = length(z); f = vector("list", n)
  run = function(i, k, cond) i >= k && all(cond(z[(i - k + 1):i]))
  for (i in seq_len(n)) {
    r = character(0)
    if (abs(z[i]) > 3) r = c(r, "1-3s")
    if (run(i, 2, function(x) x > 2) || run(i, 2, function(x) x < -2)) r = c(r, "2-2s")
    if (i >= 2 && abs(z[i] - z[i - 1]) > 4) r = c(r, "R-4s")
    if (run(i, 4, function(x) x > 1) || run(i, 4, function(x) x < -1)) r = c(r, "4-1s")
    if (run(i, 10, function(x) x > 0) || run(i, 10, function(x) x < 0)) r = c(r, "10x")
    if (!length(r) && abs(z[i]) > 2) r = "1-2s (warning)"
    f[[i]] = if (length(r)) paste(r, collapse = ", ") else ""
  }
  unlist(f)
}

# Things to check when a rule fires, by the kind of quantity. These are prompts, not diagnoses.
qc_hint = function(quantity, rule) {
  if (!nzchar(rule)) return("")
  q = tolower(quantity)
  base = if (grepl("slope", q)) "standards (age, dilution, contamination), electrode surface or membrane, reference electrode, temperature"
    else if (grepl("background|b_bkg|intercept|e0", q)) "standards and blank (contamination), reference electrode drift, filling solution"
    else if (grepl("check|standard", q)) "the check standard itself, reference electrode, reading time"
    else if (grepl("ecsa|area|charge", q)) "electrode cleaning or polishing, film deposition, electrolyte, scan window"
    else "the electrode, the solutions and the instrument settings"
  # a point beyond 3 SD (or two beyond 2 SD) is read as a sudden change even when the jump also breaks R-4s
  kind = if (grepl("1-3s|2-2s|warning", rule)) "A sudden change: " else if (grepl("R-4s", rule)) "Poor repeatability: " else "A gradual shift: "
  paste0(kind, "check ", base, ".")
}

qc_table = function(h, n_ref = NA) {
  L = qc_limits(h, n_ref)
  do.call(rbind, lapply(split(h, h$quantity), function(d) {
    d = d[order(d$date), ]; l = L[L$quantity == d$quantity[1], ]
    z = if (is.finite(l$sd) && l$sd > 0) (d$value - l$mean) / l$sd else rep(0, nrow(d))
    d$z = z; d$rule = qc_rules(z); d$check = vapply(seq_len(nrow(d)), function(i) qc_hint(d$quantity[i], d$rule[i]), ""); d
  }))
}

qc_plot = function(h, n_ref = NA) {
  t = qc_table(h, n_ref); L = qc_limits(h, n_ref)
  lines = do.call(rbind, lapply(seq_len(nrow(L)), function(i) data.frame(quantity = L$quantity[i],
    y = L$mean[i] + c(-3, -2, 0, 2, 3) * ifelse(is.finite(L$sd[i]), L$sd[i], 0), kind = c("3 SD", "2 SD", "mean", "2 SD", "3 SD"))))
  t$flag = ifelse(grepl("1-3s|2-2s|R-4s|4-1s|10x", t$rule), "rule broken", ifelse(nzchar(t$rule), "warning", "in control"))
  t$date = as.Date(t$date)
  ggplot2::ggplot(t, ggplot2::aes(date, value)) +
    ggplot2::geom_hline(data = lines, ggplot2::aes(yintercept = y, linetype = kind), colour = "grey45") +
    ggplot2::geom_line(colour = "grey60") + ggplot2::geom_point(ggplot2::aes(colour = flag), size = 2.6) +
    ggplot2::scale_colour_manual(values = c(`in control` = "#1A6FDF", warning = "#E0A100", `rule broken` = "#C0392B"), drop = FALSE) +
    ggplot2::scale_linetype_manual(values = c(mean = "solid", `2 SD` = "dashed", `3 SD` = "dotted")) +
    ggplot2::facet_wrap(~ quantity, scales = "free_y", ncol = 1) + ggplot2::labs(x = NULL, y = NULL, colour = NULL, linetype = NULL)
}
