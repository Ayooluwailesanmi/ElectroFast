# The invisible software note in downloaded figures: a PNG gets "Software" and "Creation Time" text chunks with valid
# CRCs (checked against the known CRC-32 of "IEND", AE426082), an SVG a <metadata> element, and both still open.
# Run from the project folder: Rscript tests/check_fig_stamp.R [folder to keep the test figures]
source("R/provenance.R"); source("R/app_core.R")
chk = function(ok, what) { cat(if (ok) "ok  " else "FAIL", what, "\n"); if (!ok) quit(status = 1) }
out = commandArgs(TRUE)[1]; if (is.na(out)) out = tempdir()
crc = as.numeric(png_crc32(charToRaw("IEND"))) %% 2^32
chk(sprintf("%04X%04X", as.integer(crc %/% 65536), as.integer(crc %% 65536)) == "AE426082", "CRC-32 of IEND is AE426082")
g = ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
png_f = file.path(out, "stamp_test.png"); svg_f = file.path(out, "stamp_test.svg")
fig_save(png_f, g, width = 3, height = 2, dpi = 100); fig_save(svg_f, g, width = 3, height = 2)
b = readBin(png_f, "raw", file.size(png_f)); s = rawToChar(b[b != as.raw(0)])
chk(grepl(paste0("tEXtSoftware", fig_credit_text()), s, fixed = TRUE, useBytes = TRUE), paste("PNG carries Software:", fig_credit_text()))
chk(grepl("tEXtCreation Time", s, fixed = TRUE, useBytes = TRUE), "PNG carries Creation Time")
chk(is.array(png::readPNG(png_f)), "PNG still opens (png package)")
x = readLines(svg_f, warn = FALSE)
chk(any(grepl(sprintf("<metadata>Made with %s</metadata>", fig_credit_text()), x, fixed = TRUE)), "SVG carries the metadata element")
chk(inherits(xml2::read_xml(svg_f), "xml_document"), "SVG is still valid XML")
cat("figure stamp checks passed\n")
