# BioLogic EC-Lab binary files (.mpr) are not read in this version. A .mpr file is recognised from its first bytes, so
# that the user gets a clear message instead of a failed read: export the run from EC-Lab as text (.mpt), which the
# BioLogic reader in R/read_other.R reads.

is_mpr = function(path) { b = readBin(path, "raw", 22); identical(rawToChar(b[b != as.raw(0)]), "BIO-LOGIC MODULAR FILE") }

read_biologic_mpr = function(path, name = basename(path)) {
  stop(name, " is a BioLogic EC-Lab binary file (.mpr), which ElectroFast does not read. In EC-Lab, export the run as ",
       "text (.mpt) and load that file instead.", call. = FALSE)
}
