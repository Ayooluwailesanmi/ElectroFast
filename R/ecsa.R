# Electrochemically active surface area of gold from the oxide reduction peak
# in a CV recorded in H2SO4. ECSA = Q / q_ref, with Q the charge under the
# reduction peak above a linear baseline and q_ref = 390 uC/cm2 for a monolayer
# of gold oxide (Trasatti and Petrii 1991, Pure Appl. Chem. 63, 711-734).

# Sweeps where the potential falls are the cathodic ones; the analysis uses the last
# stable cycle, so by default the last cathodic sweep in the file.
cathodic_segments = function(data) {
  s = tapply(data$E_V, data$segment, function(e) mean(diff(e)))
  as.integer(names(s)[s < 0])
}

# Baseline limits: walk out from the peak on each side until the smoothed current
# changes by less than flat_frac of the peak height per 10 mV, i.e. until the
# current has returned to a flat region.
suggest_limits = function(E, I, peak_i, flat_frac = 0.001, smooth_n = 11) {
  Is = stats::filter(I, rep(1 / smooth_n, smooth_n), sides = 2)
  Is[is.na(Is)] = I[is.na(Is)]
  step = stats::median(abs(diff(E)))
  k = max(1, round(0.01 / step))
  tol = flat_frac * abs(I[peak_i])
  walk = function(dir) {
    j = peak_i
    while (j + dir * k >= 1 && j + dir * k <= length(E)) {
      if (abs(Is[j + dir * k] - Is[j]) < tol && abs(E[j] - E[peak_i]) > 0.05) break
      j = j + dir * k
    }
    j
  }
  sort(c(E[walk(-1)], E[walk(1)]))
}

ecsa_gold = function(x, q_ref_uC_cm2 = 390, segment = NULL, limits = NULL,
                     geometric_area_cm2 = NULL) {
  d = x$data
  v = x$params$scan_rate_V_s
  if (is.null(v) || is.na(v)) stop("Scan rate unknown for ", x$file)
  cath = cathodic_segments(d)
  if (is.null(segment)) segment = max(cath)
  s = d[d$segment == segment, ]

  # Peak search skips the first and last 15% of the sweep, where the current
  # rises towards the potential limits rather than forming the oxide peak.
  n = nrow(s)
  inner = seq(ceiling(0.15 * n), floor(0.85 * n))
  peak_i = inner[which.min(s$I_A[inner])]
  # On the cathodic sweep the oxide reduction peak is a negative current. A positive peak much larger than any negative
  # current means the sign is inverted (leads or current polarity reversed): the limits would then span the whole sweep.
  inverted = max(s$I_A[inner]) > 2 * abs(min(s$I_A[inner]))

  if (is.null(limits)) limits = suggest_limits(s$E_V, s$I_A, peak_i)
  w = s$E_V >= limits[1] & s$E_V <= limits[2]
  E = s$E_V[w]
  I = s$I_A[w]
  base = stats::approx(c(E[1], E[length(E)]), c(I[1], I[length(I)]), xout = E)$y
  net = I - base

  # Charge from the potential axis: dt = dE / v, trapezoid rule.
  Q_C = abs(sum(diff(E) * (net[-1] + net[-length(net)]) / 2)) / v
  ecsa = Q_C / (q_ref_uC_cm2 * 1e-6)

  # Stability of the oxide charge over the cathodic sweeps, same limits.
  per_cycle = sapply(cath, function(g) {
    sg = d[d$segment == g & d$E_V >= limits[1] & d$E_V <= limits[2], ]
    b = stats::approx(c(sg$E_V[1], sg$E_V[nrow(sg)]), c(sg$I_A[1], sg$I_A[nrow(sg)]), xout = sg$E_V)$y
    nt = sg$I_A - b
    abs(sum(diff(sg$E_V) * (nt[-1] + nt[-length(nt)]) / 2)) / v
  })

  res = data.frame(
    file = x$file, segment = segment, scan_rate_V_s = v,
    peak_E_V = s$E_V[peak_i], peak_I_A = s$I_A[peak_i],
    lower_limit_V = limits[1], upper_limit_V = limits[2],
    Q_uC = Q_C * 1e6, q_ref_uC_cm2 = q_ref_uC_cm2, ECSA_cm2 = ecsa,
    change_last_cycle_pct = if (length(per_cycle) > 1)
      100 * (per_cycle[length(per_cycle)] / per_cycle[length(per_cycle) - 1] - 1) else NA
  )
  if (!is.null(geometric_area_cm2)) res$roughness_factor = ecsa / geometric_area_cm2
  res$note = if (inverted) "The current sign looks inverted: the reduction peak on the cathodic sweep is positive. Check the leads or the current polarity setting. This area is not valid." else ""

  chi = x$chi_peaks
  if (!is.null(chi) && any(chi$segment == segment & chi$ip_A < 0)) {
    r = chi[chi$segment == segment & chi$ip_A < 0, ]
    r = r[which.min(r$ip_A), ]
    res$chi_Ah_uC = abs(r$Ah_C) * 1e6
  }

  list(result = res, sweep = s, window = data.frame(E_V = E, I_A = I, base_A = base),
       per_cycle = data.frame(segment = cath, Q_uC = per_cycle * 1e6))
}

# Origin 2026 default line colour increment list, decoded from
# Themes/IncrementList/Line Color List.oth in the Origin installation.
origin_line_colours = c("#515151", "#F14040", "#1A6FDF", "#37AD6B", "#B177DE", "#CC9900",
                        "#00CBCC", "#7D4E4E", "#8E8E00", "#FB6501", "#6699CC", "#6FB802")

# Origin 2026 "Candy" list (Themes/Color/Candy.oth), the largest of Origin's
# lists of distinct colours; used once the default line list runs out.
origin_candy_colours = c("#EF0000", "#336699", "#FEC211", "#3BC371", "#666699", "#999999", "#FF6666",
                         "#6699CC", "#CC6600", "#009999", "#6B67BC", "#99867A", "#CC3333", "#669999",
                         "#FF9999")

# Distinct colours for the cycles other than the one used: Origin's line list
# without its dark grey (too close to the black used cycle), then Candy, with
# duplicates and the shading colour removed.
cycle_palette = function(area_fill) {
  p = unique(toupper(c(origin_line_colours[-1], origin_candy_colours)))
  p[p != toupper(area_fill)]
}

# Full CV view with the integrated oxide reduction area shaded above its
# baseline. The cycle used for ECSA is black and every other cycle gets its own
# Origin colour, listed in the legend (inside the plot for up to 8 cycles, to the
# right beyond that). Only when there are more cycles than distinct Origin
# colours are the earlier cycles coloured along the viridis scale (also an
# Origin colour map) with a colour bar, since colours would otherwise repeat.
# area_fill defaults to gold, which is kept out of the cycle colours.
plot_ecsa_cv = function(fit, x, title = NULL, reference = NULL, area_fill = "#CC9900",
                        max_inside = 8) {
  r = fit$result
  d = x$data
  w = fit$window
  ua = 1e6

  # A cycle is an anodic sweep followed by its cathodic sweep
  d$cycle = ceiling(d$segment / 2)
  used = ceiling(r$segment / 2)
  n_cyc = max(d$cycle)
  pal = cycle_palette(area_fill)
  listed = n_cyc - 1 <= length(pal)
  inside = n_cyc <= max_inside
  earlier = d[d$cycle != used, ]
  last = d[d$cycle == used, ]
  used_lab = sprintf("Cycle %d (used for ECSA)", used)

  # Direction arrows on the cycle used, away from the peaks
  an = d[d$segment == r$segment - 1, ]
  ca = d[d$segment == r$segment, ]
  arrow_at = function(s, e0, e1) {
    i0 = which.min(abs(s$E_V - e0)); i1 = which.min(abs(s$E_V - e1))
    data.frame(x = s$E_V[i0], y = s$I_A[i0] * ua, xend = s$E_V[i1], yend = s$I_A[i1] * ua)
  }
  arrows = rbind(arrow_at(an, 0.25, 0.40), arrow_at(ca, 0.45, 0.30))

  ox_form = an[an$E_V > 1.0 & an$E_V < 1.45, ]
  ox_i = which.max(ox_form$I_A)

  lab = sprintf("Q = %.2f µC\nECSA = %.4f cm²", r$Q_uC, r$ECSA_cm2)
  if (!is.null(r$roughness_factor)) lab = sprintf("%s\nRoughness factor = %.2f", lab, r$roughness_factor)

  keep = function(z) z[abs(z$E_V - r$peak_E_V) <= 0.03 & z$segment %% 2 == 0 & z$I_A <= r$peak_I_A * 0.80, ]
  tip_e = keep(earlier)
  tip_l = keep(last)

  # Earlier cycles as a listed legend or a colour bar; the inset follows the same scheme
  if (listed) {
    cyc_ids = sort(unique(earlier$cycle))
    cols = setNames(pal[seq_along(cyc_ids)], paste("Cycle", cyc_ids))
    earlier$cycle_f = factor(paste("Cycle", earlier$cycle), levels = names(cols))
    tip_e$cycle_f = factor(paste("Cycle", tip_e$cycle), levels = names(cols))
    colour_layers = list(
      ggplot2::geom_path(data = earlier, ggplot2::aes(E_V, I_A * ua, colour = cycle_f, group = cycle),
                         linewidth = 0.55),
      ggplot2::scale_colour_manual(values = cols, name = NULL,
                                   guide = ggplot2::guide_legend(order = 1, ncol = if (length(cols) > 14) 2 else 1)))
    inset_layers = list(
      ggplot2::geom_path(data = tip_e, ggplot2::aes(E_V, I_A * ua, colour = cycle_f, group = cycle),
                         linewidth = 0.6),
      ggplot2::scale_colour_manual(values = cols, guide = "none"))
  } else {
    colour_layers = list(
      ggplot2::geom_path(data = earlier, ggplot2::aes(E_V, I_A * ua, colour = cycle, group = cycle),
                         linewidth = 0.45),
      ggplot2::scale_colour_viridis_c(name = "Earlier cycles", end = 0.9,
                                      breaks = unique(round(seq(1, n_cyc - 1, length.out = 5))),
                                      guide = ggplot2::guide_colourbar(order = 1, barheight = ggplot2::unit(30, "mm"))))
    inset_layers = list(
      ggplot2::geom_path(data = tip_e, ggplot2::aes(E_V, I_A * ua, colour = cycle, group = cycle),
                         linewidth = 0.5),
      ggplot2::scale_colour_viridis_c(end = 0.9, guide = "none"))
  }

  xlab = if (is.null(reference)) "Potential (V)" else sprintf("Potential (V vs. %s)", reference)

  main = ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = w, ggplot2::aes(E_V, ymin = I_A * ua, ymax = base_A * ua),
                         fill = area_fill, alpha = 0.40) +
    colour_layers +
    # The used cycle gets its own legend entry through linetype, so it stays
    # black whichever colour scheme the earlier cycles use
    ggplot2::geom_path(data = last, ggplot2::aes(E_V, I_A * ua, linetype = used_lab), colour = "black",
                       linewidth = 0.75) +
    ggplot2::scale_linetype_manual(values = setNames("solid", used_lab), name = NULL,
                                   guide = ggplot2::guide_legend(order = 2)) +
    ggplot2::scale_x_continuous(breaks = seq(floor(min(d$E_V) * 10) / 10, ceiling(max(d$E_V) * 10) / 10, by = 0.1)) +
    ggplot2::scale_y_continuous(breaks = seq(floor(min(d$I_A) * ua / 10) * 10, ceiling(max(d$I_A) * ua / 10) * 10, by = 10)) +
    ggplot2::geom_line(data = w, ggplot2::aes(E_V, base_A * ua), linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_segment(data = arrows, ggplot2::aes(x, y, xend = xend, yend = yend),
                          arrow = ggplot2::arrow(length = ggplot2::unit(2.2, "mm"), type = "closed"),
                          linewidth = 0.5) +
    ggplot2::annotate("segment", x = r$peak_E_V + 0.03, xend = r$peak_E_V + 0.22,
                      y = r$peak_I_A * ua * 0.40, yend = r$peak_I_A * ua * 0.40, linewidth = 0.3) +
    ggplot2::annotate("text", x = r$peak_E_V + 0.23, y = r$peak_I_A * ua * 0.40, hjust = 0, size = 3.2,
                      label = sprintf("Au oxide reduction peak\nEp = %.3f V\nShaded: area used for ECSA", r$peak_E_V)) +
    ggplot2::annotate("text", x = ox_form$E_V[ox_i], y = ox_form$I_A[ox_i] * ua + 4, size = 3.2,
                      label = "Au oxide formation") +
    # Results sit top left, in the empty space above the double-layer region,
    # so the legend can run down the left side
    ggplot2::annotate("text", x = min(d$E_V) + 0.02, y = max(d$I_A) * ua, hjust = 0, vjust = 1, size = 3.2,
                      label = lab) +
    ggplot2::labs(x = xlab, y = expression("Current ("*mu*"A)"),
                  title = if (is.null(title)) "Au electrode in 0.5 M H₂SO₄" else title,
                  subtitle = sprintf("%d cycles at %.1f V/s, ECSA from cycle %d with %g µC/cm²",
                                     n_cyc, r$scan_rate_V_s, used, r$q_ref_uC_cm2)) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11),
                   plot.subtitle = ggplot2::element_text(size = 9),
                   legend.position = if (inside) c(0.02, 0.60) else "right",
                   legend.justification = if (inside) c(0, 1) else c(0, 0.5),
                   legend.key.width = ggplot2::unit(8, "mm"),
                   legend.spacing.y = ggplot2::unit(1, "mm"),
                   legend.title = ggplot2::element_text(size = 9),
                   legend.background = ggplot2::element_blank())

  # The cycles overlap almost everywhere, so an inset around the reduction peak
  # tip shows how each cycle differs.
  inset = ggplot2::ggplot() +
    inset_layers +
    ggplot2::geom_path(data = tip_l, ggplot2::aes(E_V, I_A * ua), colour = "black", linewidth = 0.7) +
    ggplot2::labs(x = xlab, y = expression("Current ("*mu*"A)"), title = "Peak tip, all cycles") +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), plot.title = ggplot2::element_text(size = 8),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA))

  main + patchwork::inset_element(inset, left = 0.64, bottom = 0.03, right = 0.99, top = 0.36)
}
