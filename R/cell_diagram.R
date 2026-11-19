# Cell diagram: a clean schematic of an electrochemical cell, drawn only from what the user states. The shapes are
# generic symbols (sheathed disk, reference tube with frit, coiled wire), not a picture of the user's glassware, and
# every label is the user's own text. A blank label shows the role alone ("Working electrode"), never a guessed material.
# Layouts: one compartment with three electrodes; two compartments with a separator (potentiometric, one electrode
# per chamber, or a potentiostat with WE and RE in chamber A and CE in chamber B); a rotating disk electrode.

cell_layouts = c("One compartment, three electrodes" = "single", "Two compartments with a separator" = "two",
                 "Two electrodes (patch-clamp amplifier)" = "pc",
                 "Rotating disk electrode" = "rde")
cell_we_kinds = c("Disk electrode (sheathed)" = "disk", "Plate or foil" = "plate", "Wire" = "wire")
cell_ce_kinds = c("Coiled wire" = "coil", "Rod" = "rod", "Plate or mesh" = "plate", "Wire" = "wire")
cell_re_kinds = c("Single-junction tube with frit" = "reference", "Double-junction tube (outer salt-bridge sleeve)" = "double",
                  "Pseudo-reference wire" = "pseudo")
cell_any_kinds = c(cell_re_kinds, cell_we_kinds, "Coiled wire" = "coil", "Rod" = "rod")

cd_col = list(line = "#2f3a45", metal = "#8a939c", dark = "#3d4650", glass = "#ffffff", sol_a = "#d6e6f5",
              sol_b = "#e2efe0", sep = "#c9822f", box = "#f2f4f6", frit = "#b9c0c7")

# Drawing pieces are collected in an environment and turned into data frames at the end
cd_new = function() { e = new.env(); e$poly = list(); e$path = list(); e$arrow = list(); e$text = list(); e }
cd_poly = function(e, x, y, fill, col = cd_col$line, lw = 0.45) { i = length(e$poly) + 1
  e$poly[[i]] = data.frame(x = x, y = y, fill = fill, col = col, lw = lw, id = i) }
cd_rect = function(e, x0, x1, y0, y1, ...) cd_poly(e, c(x0, x1, x1, x0), c(y0, y0, y1, y1), ...)
cd_path = function(e, x, y, col = cd_col$line, lw = 0.5, arrow = FALSE) { k = if (arrow) "arrow" else "path"; i = length(e[[k]]) + 1
  e[[k]][[i]] = data.frame(x = x, y = y, col = col, lw = lw, id = i) }
cd_text = function(e, x, y, label, size = 3.6, hjust = 0.5, vjust = 0.5, face = "plain")
  e$text[[length(e$text) + 1]] = data.frame(x = x, y = y, label = label, size = size, hjust = hjust, vjust = vjust, face = face)
cd_arc = function(cx, cy, r, a0, a1, n = 16) { a = seq(a0, a1, length.out = n); list(x = cx + r * cos(a), y = cy + r * sin(a)) }
cd_circle = function(e, cx, cy, r, ...) { a = cd_arc(cx, cy, r, 0, 2 * pi, 40); cd_poly(e, a$x, a$y, ...) }

# Open vessel with rounded bottom corners, and the solution inside it up to the level
cd_vessel_outline = function(x0, x1, top, r = 0.7) {
  a1 = cd_arc(x0 + r, r, r, pi, 1.5 * pi); a2 = cd_arc(x1 - r, r, r, 1.5 * pi, 2 * pi)
  list(x = c(x0, a1$x, a2$x, x1), y = c(top, a1$y, a2$y, top)) }
cd_solution = function(e, x0, x1, level, fill, round_left = TRUE, round_right = TRUE, r = 0.7) {
  left = if (round_left) cd_arc(x0 + r, r, r, pi, 1.5 * pi) else list(x = x0, y = 0)
  right = if (round_right) cd_arc(x1 - r, r, r, 1.5 * pi, 2 * pi) else list(x = x1, y = 0)
  cd_poly(e, c(x0, left$x, right$x, x1), c(level, left$y, right$y, level), fill = fill, col = NA, lw = 0) }

# One electrode at x = cx, from its tip up to top; the shape is the symbol chosen by the user
cd_electrode = function(e, kind, cx, tip, top) {
  if (cd_is3d(e) && kind %in% c("disk", "rde", "reference", "double", "rod", "plate")) return(cd_electrode3d(e, kind, cx, tip, top))
  switch(kind,
    disk = { cd_rect(e, cx - 0.35, cx + 0.35, tip, top, fill = cd_col$glass); cd_rect(e, cx - 0.24, cx + 0.24, tip, tip + 0.16, fill = cd_col$dark) },
    plate = { cd_path(e, c(cx, cx), c(tip + 2.2, top), lw = 0.6); cd_rect(e, cx - 0.5, cx + 0.5, tip, tip + 2.2, fill = cd_col$metal) },
    wire = cd_path(e, c(cx, cx), c(tip, top), col = cd_col$dark, lw = 0.9),
    rod = cd_rect(e, cx - 0.2, cx + 0.2, tip, top, fill = cd_col$dark),
    coil = { cd_path(e, c(cx, cx), c(tip + 2, top), col = cd_col$dark, lw = 0.7)
             t = seq(0, 6 * 2 * pi, length.out = 240)
             cd_path(e, cx + 0.38 * sin(t), tip + 2 - 2 * t / max(t), col = cd_col$dark, lw = 0.7) },
    reference = { cd_rect(e, cx - 0.4, cx + 0.4, tip, top, fill = "#f7fafc")
                  cd_rect(e, cx - 0.4, cx + 0.4, tip, tip + 0.35, fill = cd_col$frit)
                  cd_path(e, c(cx, cx), c(tip + 1.6, top), col = cd_col$dark, lw = 0.5)
                  cd_rect(e, cx - 0.09, cx + 0.09, tip + 0.8, tip + 1.6, fill = cd_col$dark, lw = 0.3) },
    double = { cd_rect(e, cx - 0.55, cx + 0.55, tip, tip + 4.2, fill = "#eef3f8")
               cd_rect(e, cx - 0.55, cx + 0.55, tip, tip + 0.35, fill = cd_col$frit)
               cd_rect(e, cx - 0.3, cx + 0.3, tip + 0.9, top, fill = "#f7fafc")
               cd_rect(e, cx - 0.3, cx + 0.3, tip + 0.9, tip + 1.2, fill = cd_col$frit)
               cd_path(e, c(cx, cx), c(tip + 2.4, top), col = cd_col$dark, lw = 0.5)
               cd_rect(e, cx - 0.08, cx + 0.08, tip + 1.6, tip + 2.4, fill = cd_col$dark, lw = 0.3) },
    pseudo = cd_path(e, c(cx, cx), c(tip, top), col = cd_col$dark, lw = 0.9),
    ume = { cd_rect(e, cx - 0.3, cx + 0.3, tip + 1.4, top, fill = cd_col$glass)
            cd_poly(e, c(cx - 0.3, cx - 0.1, cx + 0.1, cx + 0.3), c(tip + 1.4, tip, tip, tip + 1.4), fill = cd_col$glass, lw = 0.45)
            cd_path(e, c(cx, cx), c(tip + 0.06, top), col = cd_col$dark, lw = 0.4)
            cd_rect(e, cx - 0.05, cx + 0.05, tip, tip + 0.06, fill = cd_col$dark, lw = 0.2) },
    rde = { cd_rect(e, cx - 0.45, cx + 0.45, tip, top, fill = cd_col$glass); cd_rect(e, cx - 0.3, cx + 0.3, tip, tip + 0.16, fill = cd_col$dark)
            a = cd_arc(cx, top - 1.2, 1, 0.75 * pi, 2.25 * pi, 30); cd_path(e, a$x, top - 1.2 + (a$y - (top - 1.2)) * 0.3, lw = 0.5, arrow = TRUE) },
    stop("unknown electrode symbol: ", kind))
}

# 3D style (optional): restrained perspective. Elliptical rims, caps and liquid surface, a light highlight and a
# faint shade on each cylinder, and a box in oblique view for the two-compartment cell. No glossy gradients; the
# labels and their positions are the same as in the flat style.
cd_ry = 0.6; cd_dx = 1.2; cd_dy = 0.75
cd_ell = function(xc, yc, rx, ry, a0, a1, n = 48) { a = seq(a0, a1, length.out = n); list(x = xc + rx * cos(a), y = yc + ry * sin(a)) }
cd_is3d = function(e) identical(e$style, "3d")

# Beaker as a cylinder: curved base, rim ellipse (back half lighter), and an elliptical liquid surface
cd_vessel3d = function(e, x0, x1, top) {
  xc = (x0 + x1) / 2; rx = (x1 - x0) / 2; f = cd_ell(xc, cd_ry, rx, cd_ry, pi, 2 * pi)
  cd_path(e, c(x0, f$x, x1), c(top, f$y, top), lw = 0.9)
  b = cd_ell(xc, top, rx, cd_ry, 0, pi); cd_path(e, b$x, b$y, col = "#8a96a3", lw = 0.6)
  r = cd_ell(xc, top, rx, cd_ry, pi, 2 * pi); cd_path(e, r$x, r$y, lw = 0.9)
}
cd_solution3d = function(e, x0, x1, level, fill, surface) {
  xc = (x0 + x1) / 2; rx = (x1 - x0) / 2; f = cd_ell(xc, cd_ry, rx, cd_ry, pi, 2 * pi)
  cd_poly(e, c(x0, f$x, x1), c(level, f$y, level), fill = fill, col = NA, lw = 0)
  s = cd_ell(xc, level, rx, cd_ry, 0, 2 * pi); cd_poly(e, s$x, s$y, fill = surface, col = "#7fa7cf", lw = 0.4)
}

# Two-compartment cell as a box: liquid tops and the right face of chamber B (drawn before the electrodes), then
# the box edges and the separator set in the dividing wall (drawn after)
cd_box3d_liquid = function(e, x0, xm, x1, level) {
  cd_poly(e, c(x0, xm, xm + cd_dx, x0 + cd_dx), c(level, level, level + cd_dy, level + cd_dy), fill = "#e6f0f9", col = "#7fa7cf", lw = 0.4)
  cd_poly(e, c(xm, x1, x1 + cd_dx, xm + cd_dx), c(level, level, level + cd_dy, level + cd_dy), fill = "#edf5ec", col = "#8fbf8a", lw = 0.4)
  cd_poly(e, c(x1, x1 + cd_dx, x1 + cd_dx, x1), c(0, cd_dy, level + cd_dy, level), fill = "#d3e6d0", col = NA, lw = 0)
}
cd_box3d_frame = function(e, x0, xm, x1, top, o1, o2) {
  cd_path(e, c(x1, x1 + cd_dx, x1 + cd_dx, x1), c(0, cd_dy, top + cd_dy, top), lw = 0.9)
  cd_path(e, c(x0, x0 + cd_dx, x1 + cd_dx), c(top, top + cd_dy, top + cd_dy), col = "#8a96a3", lw = 0.6)
  cd_path(e, c(xm, xm + cd_dx), c(top, top + cd_dy), col = "#8a96a3", lw = 0.6)
  cd_poly(e, c(xm, xm + cd_dx, xm + cd_dx, xm), c(o1, o1 + cd_dy, o2 + cd_dy, o2), fill = "#c9822fb3", col = cd_col$sep, lw = 0.5)
}

# Cylinder look for a straight electrode body from y0 to y1: highlight, shade, rounded bottom edge and a top cap
cd_cyl = function(e, cx, hw, y0, y1, body, bottom = body) {
  ry = hw * 0.32
  cd_rect(e, cx - hw * 0.62, cx - hw * 0.28, y0, y1, fill = "#FFFFFFA6", col = NA, lw = 0)
  cd_rect(e, cx + hw * 0.58, cx + hw - 0.03, y0, y1, fill = "#0000001A", col = NA, lw = 0)
  b = cd_ell(cx, y0, hw, ry, pi, 2 * pi)
  cd_poly(e, c(cx - hw, b$x, cx + hw), c(y0 + 0.04, b$y, y0 + 0.04), fill = bottom, col = NA, lw = 0)
  cd_path(e, c(cx - hw, b$x, cx + hw), c(y0, b$y, y0), lw = 0.45)
  t = cd_ell(cx, y1, hw, ry, 0, 2 * pi); cd_poly(e, t$x, t$y, fill = body, lw = 0.45)
}
cd_electrode3d = function(e, kind, cx, tip, top) {
  switch(kind,
    disk = { cd_rect(e, cx - 0.35, cx + 0.35, tip, top, fill = cd_col$glass); cd_cyl(e, cx, 0.35, tip, top, cd_col$glass)
             d = cd_ell(cx, tip - 0.02, 0.26, 0.085, 0, 2 * pi); cd_poly(e, d$x, d$y, fill = cd_col$dark, lw = 0.3) },
    rde = { cd_rect(e, cx - 0.45, cx + 0.45, tip, top, fill = cd_col$glass); cd_cyl(e, cx, 0.45, tip, top, cd_col$glass)
            d = cd_ell(cx, tip - 0.02, 0.32, 0.1, 0, 2 * pi); cd_poly(e, d$x, d$y, fill = cd_col$dark, lw = 0.3)
            a = cd_arc(cx, top - 1.2, 1, 0.75 * pi, 2.25 * pi, 30); cd_path(e, a$x, top - 1.2 + (a$y - (top - 1.2)) * 0.3, lw = 0.5, arrow = TRUE) },
    reference = { cd_rect(e, cx - 0.4, cx + 0.4, tip, top, fill = "#f7fafc"); cd_rect(e, cx - 0.4, cx + 0.4, tip, tip + 0.35, fill = cd_col$frit)
                  cd_path(e, c(cx, cx), c(tip + 1.6, top), col = cd_col$dark, lw = 0.5)
                  cd_rect(e, cx - 0.09, cx + 0.09, tip + 0.8, tip + 1.6, fill = cd_col$dark, lw = 0.3)
                  cd_cyl(e, cx, 0.4, tip, top, "#f7fafc", cd_col$frit) },
    double = { cd_rect(e, cx - 0.55, cx + 0.55, tip, tip + 4.2, fill = "#eef3f8"); cd_rect(e, cx - 0.55, cx + 0.55, tip, tip + 0.35, fill = cd_col$frit)
               cd_cyl(e, cx, 0.55, tip, tip + 4.2, "#eef3f8", cd_col$frit)
               cd_rect(e, cx - 0.3, cx + 0.3, tip + 0.9, top, fill = "#f7fafc"); cd_rect(e, cx - 0.3, cx + 0.3, tip + 0.9, tip + 1.2, fill = cd_col$frit)
               cd_path(e, c(cx, cx), c(tip + 2.4, top), col = cd_col$dark, lw = 0.5)
               cd_rect(e, cx - 0.08, cx + 0.08, tip + 1.6, tip + 2.4, fill = cd_col$dark, lw = 0.3)
               cd_cyl(e, cx, 0.3, tip + 0.9, top, "#f7fafc", cd_col$frit) },
    rod = { cd_rect(e, cx - 0.2, cx + 0.2, tip, top, fill = cd_col$dark); cd_cyl(e, cx, 0.2, tip, top, cd_col$dark) },
    plate = { cd_path(e, c(cx, cx), c(tip + 2.2, top), lw = 0.6); cd_rect(e, cx - 0.5, cx + 0.5, tip, tip + 2.2, fill = cd_col$metal)
              cd_poly(e, c(cx + 0.5, cx + 0.64, cx + 0.64, cx + 0.5), c(tip, tip + 0.1, tip + 2.3, tip + 2.2), fill = "#6f7881", lw = 0.4)
              cd_poly(e, c(cx - 0.5, cx + 0.5, cx + 0.64, cx - 0.36), c(tip + 2.2, tip + 2.2, tip + 2.3, tip + 2.3), fill = "#aab2b9", lw = 0.4) },
    cd_electrode(e, kind, cx, tip, top))
}
cd_instrument3d = function(e, x0, x1, top) {
  t = top + 4.9; b = top + 2.6; dx = 0.7; dy = 0.45
  cd_poly(e, c(x0, x1, x1 + dx, x0 + dx), c(t, t, t + dy, t + dy), fill = "#e3e7eb", lw = 0.45)
  cd_poly(e, c(x1, x1 + dx, x1 + dx, x1), c(b, b + dy, t + dy, t), fill = "#d5dadf", lw = 0.45)
}

# Instrument box across the top, with one straight lead per electrode and terminal names where they apply
cd_instrument = function(e, label, xs, terminals, x0, x1, top) {
  cd_rect(e, x0, x1, top + 2.6, top + 4.9, fill = cd_col$box)
  if (cd_is3d(e)) cd_instrument3d(e, x0, x1, top)
  cd_text(e, (x0 + x1) / 2, top + 4.25, label, size = 4.8, face = "bold")
  for (i in seq_along(xs)) {
    cd_path(e, c(xs[i], xs[i]), c(top, top + 2.6), lw = 0.45)
    cd_circle(e, xs[i], top + 2.6, 0.13, fill = cd_col$line, lw = 0.3)
    if (nzchar(terminals[i])) cd_text(e, xs[i], top + 3.2, terminals[i], size = 3.9, face = "bold")
  }
}

# Digits that follow an element symbol or a closing bracket become subscripts: H2SO4 to H₂SO₄, N2 to N₂, Mg(ClO4)2 to
# Mg(ClO₄)₂. Used only for solutions, electrolyte and gas (not electrode or instrument names such as CHI760E), and only
# when the user leaves the option on. Concentrations such as "0.5 M" are untouched: their digits do not follow a letter.
cd_subscripts = function(x) {
  sub_digit = function(d) chartr("0123456789", "₀₁₂₃₄₅₆₇₈₉", d)
  m = gregexpr("(?<=[A-Za-z)\\]])[0-9]+", x, perl = TRUE)
  regmatches(x, m) = lapply(regmatches(x, m), sub_digit)
  x
}

# One entry as text for the report: the role, then the user's text if any
cd_entry = function(role, value) { v = trimws(value %||% ""); if (nzchar(v)) paste0(role, ": ", v) else role }

# Half-width of an electrode symbol at height y, so a leader starts at the edge of the part and not inside it
cd_halfwidth = function(kind, tip, y) switch(kind, disk = 0.35, ume = if (y < tip + 1.4) 0.2 else 0.3, rde = 0.45, reference = 0.4, double = if (y < tip + 4.2) 0.55 else 0.3,
  plate = if (y < tip + 2.2) 0.5 else 0, rod = 0.2, coil = if (y < tip + 2) 0.38 else 0, 0)

# Labels beside the drawing, each joined to its part by a thin leader. The anchor (ax, ay) is on the part; the leader
# runs level from the anchor to just past the vessel wall, so ay is a height where nothing else is in the way, then
# bends to its label. Labels keep the top-to-bottom order of their anchors, so leaders on one side never cross.
cd_callouts = function(e, items, side, wall, floor = -1.2) {
  if (!length(items)) return(0)
  items = items[order(-vapply(items, `[[`, 0, "ay"))]
  # role only in "full" labels; in "short" labels the role is NULL and only the description is shown (all text bold)
  lines = lapply(items, function(it) c(it$role, unlist(lapply(it$values[nzchar(trimws(it$values))], strwrap, 26))))
  h = (lengths(lines) - 1) * 0.95; ly = vapply(items, `[[`, 0, "ay")
  for (i in seq_along(items)[-1]) ly[i] = min(ly[i], ly[i - 1] - h[i - 1] / 2 - h[i] / 2 - 1.25)
  low = min(ly - h / 2); if (low < floor) ly = ly + (floor - low)
  s = if (side == "left") -1 else 1; xa = wall + s * 0.6; xb = wall + s * 1.4; xt = wall + s * 1.65
  for (i in seq_along(items)) { it = items[[i]]
    cd_path(e, c(it$ax, xa, xb), c(it$ay, it$ay, ly[i]), col = "#6b7682", lw = 0.35)
    cd_circle(e, it$ax, it$ay, 0.09, fill = "#6b7682", col = NA, lw = 0)
    y1 = ly[i] + h[i] / 2
    for (j in seq_along(lines[[i]])) cd_text(e, xt, y1 - (j - 1) * 0.95, lines[[i]][j], size = 4.2,
                                             hjust = if (side == "left") 1 else 0, face = "bold")
  }
  max(nchar(unlist(lines)))
}
cd_item = function(ax, ay, role, ...) list(ax = ax, ay = ay, role = role, values = c(...))

# spec: a list with layout and the user's entries (see the app's Cell diagram tab). Returns the plot and the entries.
cell_diagram = function(spec) {
  g = function(k, d = "") { v = spec[[k]]; if (is.null(v) || (length(v) == 1 && is.na(v))) d else v }
  f = function(k) if (isTRUE(spec$subscripts)) cd_subscripts(g(k)) else g(k)
  e = cd_new(); e$style = g("style", "flat"); is3d = cd_is3d(e); L = list(); R = list(); entries = character(0)
  # each label is added to one side and to the report entries; an optional part (analyte, gas, rotation) only when entered or ticked
  # labels: "short" shows only the user's description (RE, WE and CE at the instrument already name the electrodes),
  # "full" adds the role, "none" draws no outside labels. The report entries are the same in every case.
  labels = g("labels", "short")
  lab = function(side, ax, ay, role, value = "", extra = NULL) {
    desc = c(value, extra); has_desc = any(nzchar(trimws(desc)))
    if (labels == "full" || (labels == "short" && has_desc)) {
      it = cd_item(ax, ay, if (labels == "full") role else NULL, value, extra)
      if (side == "left") L[[length(L) + 1]] <<- it else R[[length(R) + 1]] <<- it
    }
    entries <<- c(entries, cd_entry(role, value), extra)
  }
  layout = g("layout", "single"); top = 12.5
  analyte = f("analyte"); has_analyte = nzchar(trimws(analyte))
  # electrolyte and analyte are written inside the liquid (no leaders), in a band at the bottom of each compartment that
  # the electrodes do not reach; chamber names the compartment: "single", "A" or "B"
  inside = list()
  content = function(chamber, role, value) {
    entries <<- c(entries, cd_entry(role, value))
    txt = if (labels == "full") cd_entry(role, value) else if (labels == "short") trimws(value) else ""
    if (nzchar(txt)) inside[[length(inside) + 1]] <<- list(chamber = chamber, text = txt, analyte = startsWith(role, "Analyte"))
  }
  if (layout %in% c("single", "rde")) {
    x0 = 0; x1 = 10; level = 8; xs = c(RE = 2.3, WE = 4.6, CE = 8.2)
    we_kind = if (layout == "rde") "rde" else g("we_kind", "disk"); re_kind = g("re_kind", "reference"); ce_kind = g("ce_kind", "coil")
    we_tip = if (layout == "rde") 4 else 4.2; re_tip = 4.9; ce_tip = 4
    if (is3d) cd_solution3d(e, x0, x1, level, cd_col$sol_a, "#e6f0f9") else cd_solution(e, x0, x1, level, cd_col$sol_a)
    stir = isTRUE(spec$stir) && layout == "single"; purge = isTRUE(spec$purge)
    # stir bar as a white capsule: most stir bars are PTFE-coated
    if (stir) { r = 0.19; a1 = cd_arc(4.3, 0.45, r, 0.5 * pi, 1.5 * pi); a2 = cd_arc(5.7, 0.45, r, -0.5 * pi, 0.5 * pi)
      cd_poly(e, c(a1$x, a2$x), c(a1$y, a2$y), fill = "#ffffff", lw = 0.5) }
    cd_electrode(e, we_kind, xs[["WE"]], we_tip, top); cd_electrode(e, re_kind, xs[["RE"]], re_tip, top); cd_electrode(e, ce_kind, xs[["CE"]], ce_tip, top)
    # the gas tube sits between WE and CE, so the CE leader runs clear to the right wall
    if (purge) { cd_path(e, c(6.4, 6.4), c(3.6, top), col = cd_col$dark, lw = 0.6)
      for (b in seq(4.2, 7.4, by = 0.8)) cd_circle(e, 6.4 + 0.18 * sin(b * 3), b, 0.14, fill = "#ffffff", lw = 0.35) }
    if (is3d) cd_vessel3d(e, x0, x1, 10) else { o = cd_vessel_outline(x0, x1, 10); cd_path(e, o$x, o$y, lw = 0.9)
      cd_path(e, c(x0 + 0.05, x1 - 0.05), c(level, level), col = "#7fa7cf", lw = 0.4) }
    cd_instrument(e, if (nzchar(trimws(g("instrument")))) g("instrument") else "Potentiostat", xs, c("RE", "WE", "CE"), x0 - 0.5, x1 + 0.5, top)
    lab("left", xs[["RE"]] - cd_halfwidth(re_kind, re_tip, 9), 9, "Reference electrode", g("re"))
    # the WE leader runs under the tip of the reference electrode (between 4.2 and 4.9)
    rot = trimws(g("rotation"))
    lab("left", xs[["WE"]] - cd_halfwidth(we_kind, we_tip, 4.55), 4.55, if (layout == "rde") "Rotating disk electrode" else "Working electrode",
        g("we"), if (layout == "rde" && nzchar(rot)) paste("Rotation:", rot))
    lab("right", xs[["CE"]] + cd_halfwidth(ce_kind, ce_tip, 9), 9, "Counter electrode", g("ce"))
    # the gas leader runs from the end of the tube, under the counter electrode
    if (purge) lab("right", 6.4, 3.65, "Gas purge", f("gas"))
    # the stir bar sits at the centre of the bottom; its leader runs under the gas tube and the counter electrode
    if (stir) lab("right", 5.9, 0.45, "Stir bar")
    content("single", "Electrolyte", f("electrolyte")); if (has_analyte) content("single", "Analyte", analyte)
    band = list(single = c((x0 + x1) / 2, if (stir) 1.25 else 0.75, if (purge) 3.3 else min(we_tip, ce_tip) - 0.45))
  } else if (layout == "pc") {
    # two electrodes on a patch-clamp amplifier, no counter electrode: the ultramicroelectrode on the headstage input and
    # the Ag/AgCl electrode on ground (Axon convention: potential at the headstage input with respect to ground)
    x0 = 0; x1 = 10; level = 8; xs = c(UME = 3.2, AG = 7.0)
    re_kind = g("re_kind", "reference"); ume_tip = 4.2; re_tip = 4.9
    if (is3d) cd_solution3d(e, x0, x1, level, cd_col$sol_a, "#e6f0f9") else cd_solution(e, x0, x1, level, cd_col$sol_a)
    # particles as small dots on a fixed, evenly scattered pattern (the R2 low-discrepancy sequence, no random numbers),
    # clear of the electrodes and the text
    if (isTRUE(spec$particles)) for (i in 1:60) { px = 0.6 + 8.8 * ((0.5 + i * 0.7548776662) %% 1); py = 4.0 + 3.6 * ((0.5 + i * 0.5698402910) %% 1)
      if (abs(px - xs[["UME"]]) > 0.6 && abs(px - xs[["AG"]]) > 0.8) cd_circle(e, px, py, 0.07, fill = cd_col$dark, col = NA, lw = 0) }
    cd_electrode(e, "ume", xs[["UME"]], ume_tip, top); cd_electrode(e, re_kind, xs[["AG"]], re_tip, top)
    if (is3d) cd_vessel3d(e, x0, x1, 10) else { o = cd_vessel_outline(x0, x1, 10); cd_path(e, o$x, o$y, lw = 0.9)
      cd_path(e, c(x0 + 0.05, x1 - 0.05), c(level, level), col = "#7fa7cf", lw = 0.4) }
    cd_instrument(e, if (nzchar(trimws(g("instrument")))) g("instrument") else "Patch-clamp amplifier", xs, c("Headstage", "Ground"), x0 - 0.5, x1 + 0.5, top)
    lab("left", xs[["UME"]] - cd_halfwidth("ume", ume_tip, 9), 9, "Working electrode", g("we"))
    lab("right", xs[["AG"]] + cd_halfwidth(re_kind, re_tip, 9), 9, "Reference and counter electrode", g("re"))
    content("single", "Electrolyte", f("electrolyte")); if (has_analyte) content("single", "Analyte", analyte)
    band = list(single = c((x0 + x1) / 2, 0.75, min(ume_tip, re_tip) - 0.45))
  } else if (layout == "two") {
    x0 = 0; xm = 10; x1 = 20; level = 8.5; o1 = 3; o2 = 6
    rr = if (is3d) 0.01 else 0.7
    cd_solution(e, x0, xm, level, cd_col$sol_a, round_right = FALSE, r = rr); cd_solution(e, xm, x1, level, cd_col$sol_b, round_left = FALSE, r = rr)
    if (is3d) cd_box3d_liquid(e, x0, xm, x1, level)
    potentiometric = g("mode", "potentiometric") == "potentiometric"
    if (potentiometric) {
      ka = g("a_kind", "reference"); kb = g("b_kind", "reference")
      xs = c(5, 15); cd_electrode(e, ka, 5, 6.4, top); cd_electrode(e, kb, 15, 6.4, top); terms = c("", ""); instr_default = "Potentiometer"
    } else {
      we_kind = g("we_kind", "disk"); re_kind = g("re_kind", "reference"); ce_kind = g("ce_kind", "coil")
      xs = c(3.5, 6.5, 15); cd_electrode(e, re_kind, 3.5, 6.9, top); cd_electrode(e, we_kind, 6.5, 6.2, top); cd_electrode(e, ce_kind, 15, 5.8, top)
      terms = c("RE", "WE", "CE"); instr_default = "Potentiostat"
    }
    o = cd_vessel_outline(x0, x1, 10, r = rr); cd_path(e, o$x, o$y, lw = 0.9)
    if (is3d) cd_box3d_frame(e, x0, xm, x1, 10, o1, o2)
    cd_path(e, c(xm, xm), c(0, o1), lw = 0.9); cd_path(e, c(xm, xm), c(o2, 10), lw = 0.9)
    cd_path(e, c(xm, xm), c(o1, o2), col = cd_col$sep, lw = 3)
    cd_path(e, c(x0 + 0.05, xm - 0.05), c(level, level), col = "#7fa7cf", lw = 0.4)
    cd_path(e, c(xm + 0.05, x1 - 0.05), c(level, level), col = "#8fbf8a", lw = 0.4)
    cd_instrument(e, if (nzchar(trimws(g("instrument")))) g("instrument") else instr_default, xs, terms, x0 - 0.5, x1 + 0.5, top)
    if (potentiometric) {
      lab("left", 5 - cd_halfwidth(ka, 6.4, 9.2), 9.2, "Electrode in chamber A", g("a"))
      lab("right", 15 + cd_halfwidth(kb, 6.4, 9.2), 9.2, "Electrode in chamber B", g("b"))
    } else {
      lab("left", 3.5 - cd_halfwidth(re_kind, 6.9, 9.2), 9.2, "Reference electrode", g("re"))
      # the WE leader runs under the tip of the reference electrode
      lab("left", 6.5 - cd_halfwidth(we_kind, 6.2, 6.55), 6.55, "Working electrode", g("we"))
      lab("right", 15 + cd_halfwidth(ce_kind, 5.8, 9.2), 9.2, "Counter electrode", g("ce"))
    }
    # the separator leader runs under every electrode tip in chamber A and above the text band
    lab("left", xm, 4.5, "Separator", g("separator"))
    content("A", "Chamber A", f("sol_a")); content("B", "Chamber B", f("sol_b"))
    in_b = g("analyte_ch", "B") == "B"
    if (has_analyte) content(if (in_b) "B" else "A", paste("Analyte, chamber", if (in_b) "B" else "A"), analyte)
    band = list(A = c((x0 + xm) / 2, 0.75, 4.05), B = c((xm + x1) / 2, 0.75, 4.05))
  } else stop("unknown layout: ", layout)

  # text inside the liquid: bold, centred in the compartment, lines stacked upward from the bottom of the band and
  # wrapped to the compartment width
  for (ch in names(band)) {
    its = Filter(function(it) it$chamber == ch, inside)
    txt = vapply(its, `[[`, "", "text")
    if (!length(txt)) next
    # short labels: a centred "+" on its own line between the electrolyte and the analyte, as methods sections write
    # "0.1 M KOH + 50 mM ethanol"; each entry stays whole and wraps cleanly
    is_an = vapply(its, `[[`, TRUE, "analyte")
    plus = labels == "short" && any(is_an) && any(!is_an)
    # largest size at which the lines fit between the bottom of the band and its top limit; smaller text wraps wider
    for (sz in c(5, 4.4, 3.8, 3.3)) {
      wr = round(22 * 5 / sz); sp = 0.9 * sz / 5
      w = if (plus) c(unlist(lapply(txt[!is_an], strwrap, wr)), "+", unlist(lapply(txt[is_an], strwrap, wr))) else unlist(lapply(txt, strwrap, wr))
      if (band[[ch]][2] + (length(w) - 1) * sp + 0.25 <= band[[ch]][3]) break
    }
    for (j in seq_along(w)) cd_text(e, band[[ch]][1], band[[ch]][2] + (length(w) - j) * sp, w[j], size = sz, face = "bold")
  }

  wl = cd_callouts(e, L, "left", x0); wr = cd_callouts(e, R, "right", x1 + if (is3d && layout == "two") cd_dx else 0)
  if (is3d && layout == "two") x1 = x1 + cd_dx
  # about 0.33 plot units per bold character at this text size, so the labels fit inside the plot
  xl = x0 - 1.65 - wl * 0.33 - 0.3; xr = x1 + 1.65 + wr * 0.33 + 0.3
  title = trimws(g("title")); if (nzchar(title)) cd_text(e, min(xl, x0 - 0.5), top + 6.2, title, size = 5.4, hjust = 0, face = "bold")
  ys = c(unlist(lapply(e$text, `[[`, "y")), 0)

  P = do.call(rbind, e$poly); Ln = do.call(rbind, e$path); A = if (length(e$arrow)) do.call(rbind, e$arrow); Tx = do.call(rbind, e$text)
  p = ggplot2::ggplot() +
    ggplot2::geom_polygon(data = P, ggplot2::aes(x, y, group = id, fill = fill, colour = col, linewidth = lw)) +
    ggplot2::geom_path(data = Ln, ggplot2::aes(x, y, group = id, colour = col, linewidth = lw), lineend = "round") +
    (if (!is.null(A)) ggplot2::geom_path(data = A, ggplot2::aes(x, y, group = id, colour = col, linewidth = lw),
                                         arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed"))) +
    ggplot2::geom_text(data = Tx, ggplot2::aes(x, y, label = label, size = size, hjust = hjust, vjust = vjust, fontface = face), colour = cd_col$line) +
    ggplot2::scale_fill_identity() + ggplot2::scale_colour_identity() + ggplot2::scale_linewidth_identity() + ggplot2::scale_size_identity() +
    ggplot2::coord_fixed(clip = "off", xlim = c(min(xl, x0 - 1), xr), ylim = c(min(-0.8, min(ys) - 0.6), top + 6.9)) +
    ggplot2::theme_void() + ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  list(plot = p, key = data.frame(Part = entries, check.names = FALSE))
}

# The Cell diagram tab. Text fields start empty: a blank field shows the role only in the drawing.
panel_cell = function() nav_panel("Cell diagram",
  layout_sidebar(
    sidebar = sidebar(width = 340,
      selectInput("cd_layout", "Cell", cell_layouts),
      conditionalPanel("input.cd_layout == 'two'",
        radioButtons("cd_mode", "Set-up", c("Potentiometric: one electrode in each chamber" = "potentiometric",
                                            "Potentiostat: working and reference in A, counter in B" = "potentiostat"))),
      conditionalPanel("input.cd_layout == 'two' && input.cd_mode == 'potentiometric'",
        selectInput("cd_a_kind", "Electrode in chamber A", cell_any_kinds), textInput("cd_a", "Its label", "", placeholder = "e.g. Ag/AgCl"),
        selectInput("cd_b_kind", "Electrode in chamber B", cell_any_kinds), textInput("cd_b", "Its label", "", placeholder = "e.g. Ag/AgCl")),
      conditionalPanel("input.cd_layout != 'two' || input.cd_mode == 'potentiostat'",
        conditionalPanel("['rde', 'pc'].indexOf(input.cd_layout) < 0", selectInput("cd_we_kind", "Working electrode", cell_we_kinds)),
        textInput("cd_we", "Working electrode label", "", placeholder = "e.g. Au disk, 2 mm"),
        selectInput("cd_re_kind", "Reference electrode", cell_re_kinds),
        textInput("cd_re", "Reference electrode label", "", placeholder = "e.g. Ag/AgCl (1 M KCl)"),
        conditionalPanel("input.cd_layout != 'pc'",
          selectInput("cd_ce_kind", "Counter electrode", cell_ce_kinds),
          textInput("cd_ce", "Counter electrode label", "", placeholder = "e.g. Pt wire"))),
      conditionalPanel("input.cd_layout == 'pc'", checkboxInput("cd_particles", "Show particles in the solution", FALSE)),
      textInput("cd_analyte", "Analyte (optional)", "", placeholder = "e.g. 50 mM ethanol"),
      conditionalPanel("input.cd_layout == 'two'", radioButtons("cd_analyte_ch", "Analyte in chamber", c("A", "B"), selected = "B", inline = TRUE)),
      conditionalPanel("input.cd_layout != 'two'",
        textInput("cd_electrolyte", "Electrolyte", "", placeholder = "e.g. 0.5 M H2SO4")),
      conditionalPanel("input.cd_layout == 'two'",
        textInput("cd_sol_a", "Chamber A solution", "", placeholder = "e.g. 1 mM LiClO4"),
        textInput("cd_sol_b", "Chamber B solution", "", placeholder = "e.g. LiClO4 standards"),
        textInput("cd_sep", "Separator", "", placeholder = "e.g. membrane, frit or glass-ceramic sheet")),
      conditionalPanel("input.cd_layout == 'rde'", textInput("cd_rot", "Rotation", "", placeholder = "e.g. 400 to 2500 rpm")),
      conditionalPanel("['single', 'rde'].indexOf(input.cd_layout) >= 0",
        checkboxInput("cd_purge", "Gas purge tube", FALSE),
        conditionalPanel("input.cd_purge", textInput("cd_gas", "Gas", "", placeholder = "e.g. N2"))),
      conditionalPanel("input.cd_layout == 'single'", checkboxInput("cd_stir", "Stir bar", FALSE)),
      textInput("cd_instr", "Instrument", "", placeholder = "Potentiostat or potentiometer"),
      textInput("cd_title", "Title (optional)", ""),
      checkboxInput("cd_sub", "Show formula numbers as subscripts: H2SO4 becomes H₂SO₄ (electrolyte, solution and gas fields only)", TRUE),
      radioButtons("cd_style", "Style", c("Flat" = "flat", "3D-style" = "3d"), inline = TRUE),
      radioButtons("cd_labels", "Labels beside the drawing", c("Short: your descriptions only" = "short",
                   "Full: part name and description" = "full", "None: drawing only" = "none")),
      checkboxInput("cd_report", "Add the diagram to the report", FALSE),
      fileInput("cd_photo", "Photo of the cell (optional, shown for comparison)", accept = c(".png", ".jpg", ".jpeg")),
      helpText("The drawing shows only what is entered here, with generic symbols for each part. It is a schematic, not a",
               "copy of the photo: check it against the cell before using it."),
      downloadButton("cd_png", "PNG", class = "btn-sm"), downloadButton("cd_pdf", "PDF", class = "btn-sm"),
      downloadButton("cd_svg", "SVG", class = "btn-sm")
    ),
    card(card_header("Schematic"), plotOutput("cd_plot", height = "520px")),
    conditionalPanel("output.cd_has_photo", card(card_header("Photo"), imageOutput("cd_photo_img", height = "auto")))
  ))
