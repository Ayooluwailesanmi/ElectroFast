# ElectroFast

ElectroFast is an R/Shiny application for electrochemical data from any potentiostat and any electrode. It covers
surface area, CV kinetics, Tafel analysis, pulse voltammetry, OCP, chronoamperometry (including single-entity spikes),
chronocoulometry, chronopotentiometry, RDE/RRDE, impedance fitting, ISE calibration, selectivity, precision and QC
control charts, with HTML, PDF and PowerPoint reports. Version 1.0.0.

## Running it

Install R (tested with R 4.5.3) and the packages:

```r
install.packages(c("shiny", "bslib", "ggplot2", "DT", "lme4", "jsonlite", "rmarkdown", "knitr", "htmltools", "zip",
                   "patchwork", "officer", "xml2", "readABF", "gridExtra", "svglite"))
```

Then, from this folder, run `shiny::runApp()`, or open `app.R` in RStudio and click Run App. The Guide tab in the app
says what each tab needs, what it gives and what to check. For a walkthrough, unzip the example data and follow the
eight tests in its README.txt. Each test has a known answer.

## Files it reads

| Instrument | Files | Reader |
|---|---|---|
| CH Instruments | text exports with or without the header, every technique | `R/read_chi_txt.R` |
| PalmSens | PSTrace sessions (.pssession) | `R/read_pssession.R` |
| BioLogic | EC-Lab text (.mpt) | `R/read_other.R` |
| Gamry | Framework files (.DTA) | `R/read_other.R` |
| Pine | AfterMath text exports | `R/read_pine_txt.R` |
| Lawson Labs | EMF16 recordings | `R/read_emf16_grf.R` |
| Axon | patch-clamp recordings (.abf) | `R/read_abf.R` |
| Any other | text or CSV with a header row (columns found by name and unit) | `R/read_other.R` |

**BioLogic users:** ElectroFast does not read EC-Lab binary files (.mpr). In EC-Lab, export each run as text (.mpt)
and load the .mpt file. If a .mpr file is loaded, ElectroFast says so and asks for the .mpt export.

## Example data

`www/examples/ElectroFast_example_data_simulated.zip` holds simulated data with known answers: cyclic voltammograms
(made with the eChem package), a Randles impedance spectrum, a Cottrell transient, current spikes, current steps,
precision studies and a QC history. Its `README.txt` describes each file. No measured data are included.

## Windows version

A Windows installer is attached to each release. It includes its own copy of R and the packages, so nothing
else needs to be installed. The installer is not code-signed, so Windows may show a SmartScreen warning when it has
been downloaded. Click "More info" and then "Run anyway".

## Browser version

ElectroFast also runs in the browser, with nothing to install, at
https://ayooluwailesanmi.github.io/ElectroFast/. R runs inside the browser, so the files you load stay on your
computer. It works in Chrome and Edge, not in Safari. The first start takes about a minute.

## Tests

Each file in `tests/` is run from this folder with `Rscript`. `tests/check_other_readers.R` takes a folder of public
example files as an argument: the test files of yadg (GPL-3) for BioLogic .mpt and of gamry-parser (MIT) for Gamry. The other tests make their own simulated data or use the example data above.

## Licence and terms of use

Copyright (C) 2026 Ayooluwa Ilesanmi.

ElectroFast is publicly available source code, free for academic, educational and noncommercial research use, under
the PolyForm Noncommercial License 1.0.0 (see `LICENSE`). Commercial use requires permission. For a commercial
licence, contact ailesanmi@missouri.edu.

The name ElectroFast refers to this software. The licence does not give permission to use the name for other or
modified software.

These terms apply to the code of ElectroFast only. R and the R packages that ElectroFast uses keep their own
licences. The Windows installer includes `THIRD_PARTY_PACKAGES.csv`, which lists the bundled packages with their versions and
licences.
