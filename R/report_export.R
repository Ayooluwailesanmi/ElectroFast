# Slide deck (PowerPoint, built with the officer package) and PDF of an ElectroFast analysis: a title slide, the
# settings, then for every tab with data its summary, plots (one per slide) and tables (split over slides
# when long). The same parts feed the HTML report, so both always show the same results.

# The deck follows the laboratory's own report decks (EtOH fall report 2, 21 Oct 2026): 16:9 white slides in Calibri
# with dark text, a title slide listing the cell, a main-points slide, then one slide per plot with the plot large
# and centred and the result as bullets under it. Every slide is a blank layout with text boxes at fixed positions,
# so the deck looks the same in PowerPoint, WPS and LibreOffice.
pptx_ink = "#111111"
pptx_par = function(txt, size, bold = FALSE, after = 6)
  officer::fpar(officer::ftext(txt, officer::fp_text(font.size = size, bold = bold, color = pptx_ink, font.family = "Calibri")),
                fp_p = officer::fp_par(padding.bottom = after))
pptx_box = function(doc, paras, left, top, width, height)
  officer::ph_with(doc, do.call(officer::block_list, paras), location = officer::ph_location(left = left, top = top, width = width, height = height))
pptx_bullets = function(lines, size) lapply(as.character(lines), function(l) pptx_par(paste0("•  ", l), size))
# rough number of printed lines, to decide what fits: Calibri averages about 0.5 em per character
pptx_lines = function(lines, size, width = 12) sum(ceiling(nchar(as.character(lines)) / (width * 72 / (0.5 * size))))
pptx_new = function(doc, title) {
  doc = officer::add_slide(doc, layout = "Blank", master = "Office Theme")
  pptx_box(doc, list(pptx_par(title, 24, bold = TRUE, after = 0)), 0.6, 0.25, 12.1, 0.8)
}
pptx_text = function(doc, title, lines) {
  # 15 pt bullets, as on the text slides of the report decks; long lists continue on a second slide
  per = 16; groups = list(character(0)); used = 0
  for (l in as.character(lines)) { n = pptx_lines(l, 15)
    if (used + n > per && length(groups[[length(groups)]])) { groups[[length(groups) + 1]] = character(0); used = 0 }
    groups[[length(groups)]] = c(groups[[length(groups)]], l); used = used + n }
  for (k in seq_along(groups)) {
    doc = pptx_new(doc, if (length(groups) > 1) sprintf("%s (%d of %d)", title, k, length(groups)) else title)
    doc = pptx_box(doc, pptx_bullets(groups[[k]], 15), 0.7, 1.2, 12, 5.8)
  }
  doc
}
pptx_plot = function(doc, title, g, caption = NULL) {
  f = tempfile(fileext = ".png"); ggplot2::ggsave(f, g, width = 10, height = 5, dpi = 200)
  doc = pptx_new(doc, title)
  doc = officer::ph_with(doc, officer::external_img(f, width = 10, height = 5), location = officer::ph_location(left = 1.67, top = 1.2, width = 10, height = 5))
  if (length(caption)) doc = pptx_box(doc, pptx_bullets(caption, 13), 0.6, 6.3, 12.1, 1.1)
  doc
}
pptx_table = function(doc, title, d, rows_per_slide = 12, max_cols = 8, caption = NULL) {
  d = as.data.frame(d); if (exists("nice_names")) names(d) = nice_names(names(d))
  if (ncol(d) > max_cols) d = d[, seq_len(max_cols)]
  d[] = lapply(d, function(v) if (is.numeric(v)) format(signif(v, 4)) else as.character(v))
  chunks = split(seq_len(nrow(d)), ceiling(seq_len(nrow(d)) / rows_per_slide))
  for (k in seq_along(chunks)) {
    doc = pptx_new(doc, if (length(chunks) > 1) sprintf("%s (%d of %d)", title, k, length(chunks)) else title)
    doc = officer::ph_with(doc, d[chunks[[k]], , drop = FALSE], alignment = rep("l", ncol(d)), location = officer::ph_location(left = 0.6, top = 1.2, width = 12.1, height = 0.3 * (length(chunks[[k]]) + 1)))
    # the result goes under the table on the first slide, as in the report decks
    if (k == 1 && length(caption)) doc = pptx_box(doc, pptx_bullets(caption, 15), 0.7, 1.2 + 0.3 * (length(chunks[[k]]) + 1) + 0.4, 12, 7.1 - (1.6 + 0.3 * (length(chunks[[k]]) + 1)))
  }
  doc
}
# A section: list(title, text, plots = named list of ggplots, tables = named list of data frames). Its text goes under
# the first plot (or the first table when there is no plot) when it fits there, and on a slide of its own otherwise.
pptx_section = function(doc, s) {
  txt = as.character(s$text); first_tab = if (length(s$tables)) nrow(as.data.frame(s$tables[[1]])) else Inf
  under_plot = length(txt) && length(s$plots) && pptx_lines(txt, 13) <= 4
  under_tab = length(txt) && !length(s$plots) && first_tab <= 8 && pptx_lines(txt, 15) <= 14 - first_tab
  if (length(txt) && !under_plot && !under_tab) doc = pptx_text(doc, s$title, txt)
  for (i in seq_along(s$plots)) doc = pptx_plot(doc, paste0(s$title, ": ", names(s$plots)[i]), s$plots[[i]], if (i == 1 && under_plot) txt)
  for (i in seq_along(s$tables)) doc = pptx_table(doc, paste0(s$title, ": ", names(s$tables)[i]), s$tables[[i]], caption = if (i == 1 && under_tab) txt)
  doc
}
# Main points: the summary lines of every section, under the section's name
pptx_main_points = function(doc, sections) {
  s = Filter(function(x) length(x$text) > 0, sections); if (!length(s)) return(doc)
  pages = list(list()); used = 0
  for (x in s) { n = 1 + pptx_lines(x$text, 15)
    if (used + n > 16 && length(pages[[length(pages)]])) { pages[[length(pages) + 1]] = list(); used = 0 }
    pages[[length(pages)]] = c(pages[[length(pages)]], list(pptx_par(x$title, 15, bold = TRUE, after = 2)), pptx_bullets(x$text, 15)); used = used + n }
  for (k in seq_along(pages)) {
    doc = pptx_new(doc, if (length(pages) > 1) sprintf("Main points (%d of %d)", k, length(pages)) else "Main points")
    doc = pptx_box(doc, pages[[k]], 0.7, 1.2, 12, 5.8)
  }
  doc
}
# details: named character vector for the title slide (for example Electrolyte, Working electrode); empty entries are left out
build_pptx = function(file, title, author, settings, sections, details = NULL) {
  tpl = "report/electrofast_slides_16x9.pptx"
  doc = if (file.exists(tpl)) officer::read_pptx(tpl) else officer::read_pptx()
  details = details[nzchar(trimws(details))]
  doc = officer::add_slide(doc, layout = "Blank", master = "Office Theme")
  lines = c(if (nzchar(author)) author, if (length(details)) paste0(names(details), ": ", details), paste("Date:", format(Sys.Date(), "%m/%d/%Y")), "Made with ElectroFast")
  doc = pptx_box(doc, c(list(pptx_par(title, 26, bold = TRUE)), lapply(lines, pptx_par, size = 26)), 0.8, 2.0, 11.7, 4.1)
  doc = pptx_main_points(doc, sections)
  for (s in sections) doc = pptx_section(doc, s)
  if (length(settings)) doc = pptx_table(doc, "Settings and files used", data.frame(Setting = names(settings), Value = vapply(settings, function(v) paste(v, collapse = ", "), "")))
  print(doc, target = file)
  pptx_table_font(file)
  invisible(file)
}
# officer writes table runs without run properties, so they take PowerPoint's 18 pt default. Give each table run
# 13 pt Calibri, as in the report decks (the header row stays bold through the table style).
pptx_table_font = function(file) {
  file = normalizePath(file); d = tempfile(); utils::unzip(file, exdir = d)
  for (f in list.files(file.path(d, "ppt", "slides"), "[.]xml$", full.names = TRUE)) {
    x = paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    if (!grepl("<a:tbl>", x, fixed = TRUE)) next
    m = gregexpr("(?s)<a:tbl>.*?</a:tbl>", x, perl = TRUE)
    regmatches(x, m) = list(gsub("<a:r>([[:space:]]*)<a:t(/?)>", "<a:r><a:rPr lang=\"en-US\" sz=\"1300\" dirty=\"0\"><a:latin typeface=\"Calibri\"/></a:rPr>\\1<a:t\\2>",
                                 regmatches(x, m)[[1]], perl = TRUE))
    writeLines(enc2utf8(x), f, useBytes = TRUE)
  }
  unlink(file); old = setwd(d); on.exit(setwd(old))
  zip::zipr(file, list.files(".", all.files = TRUE, no.. = TRUE))
}

# The same sections as a PDF (landscape Letter pages, base PDF device): a title page, a page of text for each
# section, one plot per page, and tables split over pages. No LaTeX or browser is needed (Cairo PDF device), so it also works
# on shinyapps.io.
pdf_page_text = function(title, lines, size = 11) {
  grid::grid.newpage()
  grid::grid.text(title, x = 0.05, y = 0.94, just = "left", gp = grid::gpar(fontsize = 18, fontface = "bold"))
  if (length(lines)) {
    w = unlist(lapply(as.character(lines), function(l) c(strwrap(l, 115, prefix = "", initial = "• ", exdent = 2), "")))
    grid::grid.text(paste(w, collapse = "\n"), x = 0.05, y = 0.86, just = c("left", "top"), gp = grid::gpar(fontsize = size, lineheight = 1.15))
  }
}
pdf_page_plot = function(title, g) {
  grid::grid.newpage()
  grid::grid.text(title, x = 0.05, y = 0.95, just = "left", gp = grid::gpar(fontsize = 15, fontface = "bold"))
  print(g, vp = grid::viewport(x = 0.5, y = 0.46, width = 0.92, height = 0.86))
}
pdf_page_table = function(title, d, rows_per_page = 16, max_cols = 9) {
  d = as.data.frame(d); if (exists("nice_names")) names(d) = nice_names(names(d)); if (ncol(d) > max_cols) d = d[, seq_len(max_cols)]
  d[] = lapply(d, function(v) if (is.numeric(v)) format(signif(v, 4)) else substr(as.character(v), 1, 45))
  chunks = split(seq_len(nrow(d)), ceiling(seq_len(nrow(d)) / rows_per_page))
  for (k in seq_along(chunks)) {
    grid::grid.newpage()
    grid::grid.text(if (length(chunks) > 1) sprintf("%s (%d of %d)", title, k, length(chunks)) else title, x = 0.05, y = 0.95, just = "left", gp = grid::gpar(fontsize = 15, fontface = "bold"))
    tg = gridExtra::tableGrob(d[chunks[[k]], , drop = FALSE], rows = NULL, theme = gridExtra::ttheme_minimal(base_size = 10))
    h_in = sum(grid::convertHeight(tg$heights, "in", valueOnly = TRUE)); w_in = sum(grid::convertWidth(tg$widths, "in", valueOnly = TRUE))
    grid::pushViewport(grid::viewport(x = 0.5, y = grid::unit(0.88, "npc"), just = c("centre", "top"), height = grid::unit(h_in, "in"), width = grid::unit(w_in, "in")))
    grid::grid.draw(tg); grid::popViewport()
  }
}
# PDF device for every PDF the app writes. On Windows the default sans font (Arial) has no subscript digits, so H₂SO₄
# printed with boxes; Segoe UI has them. Elsewhere the system default font is kept (fontconfig finds a font that has them).
echem_pdf = function(...) if (.Platform$OS.type == "windows") grDevices::cairo_pdf(..., family = "Segoe UI") else grDevices::cairo_pdf(...)

build_pdf = function(file, title, author, settings, sections) {
  echem_pdf(file, width = 11, height = 8.5, onefile = TRUE)
  on.exit(grDevices::dev.off())
  grid::grid.newpage()
  grid::grid.text(title, y = 0.6, gp = grid::gpar(fontsize = 26, fontface = "bold"))
  grid::grid.text(paste(c(if (nzchar(author)) author, format(Sys.time(), "%d %B %Y"), "Made with ElectroFast"), collapse = "\n"), y = 0.42, gp = grid::gpar(fontsize = 13, col = "grey30"))
  if (length(settings)) pdf_page_table("Settings", data.frame(Setting = names(settings), Value = vapply(settings, function(v) paste(v, collapse = ", "), "")))
  for (s in sections) {
    pdf_page_text(s$title, s$text)
    for (n in names(s$plots)) pdf_page_plot(paste0(s$title, ": ", n), s$plots[[n]])
    for (n in names(s$tables)) pdf_page_table(paste0(s$title, ": ", n), s$tables[[n]])
  }
  invisible(file)
}

# The same sections as one HTML file built directly in R. Used where Pandoc is not available (the browser version
# of the app, Shinylive), so rmarkdown cannot render: plots are embedded as PNG images and tables as HTML tables.
html_table = function(d) {
  d = as.data.frame(d); if (exists("nice_names")) names(d) = nice_names(names(d)); d[] = lapply(d, function(v) if (is.numeric(v)) format(signif(v, 4)) else as.character(v))
  htmltools::tags$table(htmltools::tags$thead(htmltools::tags$tr(lapply(names(d), htmltools::tags$th))),
    htmltools::tags$tbody(lapply(seq_len(nrow(d)), function(i) htmltools::tags$tr(lapply(d[i, ], htmltools::tags$td)))))
}
html_plot = function(g) {
  f = tempfile(fileext = ".png"); ggplot2::ggsave(f, g, width = 9, height = 5.2, dpi = 150)
  htmltools::tags$img(src = paste0("data:image/png;base64,", jsonlite::base64_enc(readBin(f, "raw", file.size(f)))), style = "max-width:100%")
}
build_html = function(file, title, author, settings, sections) {
  # the section list sits on the left and stays in view while scrolling, as in the Rmd report; on a narrow screen it goes on top
  css = "body{font-family:'Source Sans 3',Arial,sans-serif;max-width:1280px;margin:24px auto;padding:0 16px;color:#1f2933;display:grid;grid-template-columns:230px minmax(0,1fr);gap:32px}
    nav{position:sticky;top:16px;align-self:start;max-height:calc(100vh - 32px);overflow-y:auto;border:1px solid #d0d7de;border-radius:4px;padding:8px 0;font-size:14px}
    nav ul{list-style:none;margin:0;padding:0} nav a{display:block;padding:5px 14px;color:#1f2933;text-decoration:none} nav a:hover{background:#f5f7fa;color:#0b5cad}
    main{min-width:0} h1{margin:0 0 4px} h2{border-top:1px solid #d0d7de;padding-top:16px;margin-top:32px} .meta{color:#52606d}
    table{border-collapse:collapse;margin:8px 0 20px;font-size:14px;display:block;overflow-x:auto} th,td{border-bottom:1px solid #e4e7eb;padding:4px 10px;text-align:left}
    th{background:#f5f7fa} img{max-width:100%}
    @media (max-width:800px){body{display:block} nav{position:static;max-height:none;margin-bottom:16px}}"
  ids = sprintf("s%d", seq_along(sections))
  # the head is written as text: htmltools drops a tags$head() placed inside tags$html() when it renders
  nav_items = c(if (length(settings)) list(c("settings", "Settings")), lapply(seq_along(sections), function(i) c(ids[i], sections[[i]]$title)))
  body = htmltools::tags$body(
    htmltools::tags$nav(htmltools::tags$ul(lapply(nav_items, function(x) htmltools::tags$li(htmltools::tags$a(href = paste0("#", x[1]), x[2]))))),
    htmltools::tags$main(
      htmltools::tags$h1(title),
      htmltools::tags$p(class = "meta", paste(c(if (nzchar(author)) author, format(Sys.time(), "%d %B %Y"), "Made with ElectroFast"), collapse = " | ")),
      if (length(settings)) list(htmltools::tags$h2(id = "settings", "Settings"), html_table(data.frame(Setting = names(settings), Value = vapply(settings, function(v) paste(v, collapse = ", "), ""), check.names = FALSE))),
      lapply(seq_along(sections), function(i) { s = sections[[i]]
        list(htmltools::tags$h2(id = ids[i], s$title),
             if (length(s$text)) htmltools::tags$ul(lapply(as.character(s$text), htmltools::tags$li)),
             lapply(names(s$plots), function(n) list(htmltools::tags$h3(n), html_plot(s$plots[[n]]))),
             lapply(names(s$tables), function(n) list(htmltools::tags$h3(n), html_table(s$tables[[n]])))) })))
  writeLines(c("<!DOCTYPE html>", "<html lang=\"en\">", "<head>", "<meta charset=\"utf-8\">", paste0("<title>", htmltools::htmlEscape(title), "</title>"),
               paste0("<style>", css, "</style>"), "</head>", as.character(body), "</html>"), file, useBytes = TRUE)
  invisible(file)
}
