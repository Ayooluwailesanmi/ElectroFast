# Axon Binary Format (.abf) files from pCLAMP / Clampex. Patch-clamp amplifiers are also used as low-noise
# potentiostats, for example for single-particle impact chronoamperometry at an ultramicroelectrode: a constant
# potential and the current recorded against time. The data are read with the readABF package (GPL-3, pure R).
# The holding potential, the amplifier filter and the protocol name, which readABF does not return, are read from
# the ABF2 header (section map at byte 76, ADC and DAC records and the strings section, as in the ABF2 layout that
# pyABF also follows). The string list is aligned on the channel units that readABF returns, so a wrong offset stops
# the read of those settings instead of giving wrong ones.

is_abf = function(path) {
  s = tryCatch(readBin(path, "raw", 4), error = function(e) raw(0))
  length(s) == 4 && rawToChar(s) %in% c("ABF2", "ABF ")
}

# current units to amperes
abf_current_scale = c(A = 1, mA = 1e-3, uA = 1e-6, "µA" = 1e-6, nA = 1e-9, pA = 1e-12, fA = 1e-15)

# Settings from an ABF2 header: telegraph filter, holding potential of DAC 0 with its units, protocol file name
abf2_settings = function(path, adc_units) {
  h = readBin(path, "raw", 65536)
  i32 = function(o) readBin(h[(o + 1):(o + 4)], "integer", size = 4, endian = "little")
  f32 = function(o) readBin(h[(o + 1):(o + 4)], "numeric", size = 4, endian = "little")
  sec = function(i) c(offset = i32(76 + 16 * i) * 512, bytes = i32(76 + 16 * i + 4), n = i32(76 + 16 * i + 8))
  adc = sec(1); dac = sec(2); str = sec(9)
  out = list(filter_Hz = f32(adc[["offset"]] + 10), lowpass_Hz = f32(adc[["offset"]] + 56), gain = f32(adc[["offset"]] + 6),
             holding = f32(dac[["offset"]] + 12), holding_units = NA_character_, protocol = NA_character_)
  # strings: the block split at zero bytes; index 1 is aligned so that the ADC units index points at the units readABF gave
  if (str[["offset"]] > 0 && str[["offset"]] + str[["bytes"]] <= length(h)) {
    b = h[(str[["offset"]] + 1):(str[["offset"]] + str[["bytes"]])]
    # split the bytes at zeros directly (the block starts with binary header bytes, which are not valid text)
    z = which(b == as.raw(0)); st0 = c(1, z + 1); en0 = c(z - 1, length(b))
    s = vapply(seq_along(st0), function(k) if (en0[k] >= st0[k]) iconv(rawToChar(b[st0[k]:en0[k]]), "latin1", "UTF-8") else "", "")
    u_idx = i32(adc[["offset"]] + 78)          # lADCUnitsIndex of ADC 0
    pos = which(s == adc_units)
    if (length(pos)) {
      off = pos[1] - u_idx
      get = function(k) if (k > 0 && off + k <= length(s)) s[off + k] else NA_character_
      out$holding_units = get(i32(dac[["offset"]] + 28))      # lDACChannelUnitsIndex of DAC 0
      p = get(i32(72))                                          # uProtocolPathIndex in the file header (byte 72)
      if (!is.na(p)) out$protocol = basename(gsub("\\\\", "/", p))
    }
  }
  out
}

# One file to the same structure as the other readers: data with t_s and I_A (and sweep when there are several),
# settings in params. The current channel is the first channel recorded in a current unit.
read_abf = function(path, name = basename(path)) {
  a = readABF::readABF(path)
  units = trimws(a$channelUnits)
  ch = which(units %in% names(abf_current_scale))[1]
  if (is.na(ch)) stop(name, ": no channel recorded in a current unit (channels: ", paste(a$channelNames, units, collapse = ", "), ")")
  dt = a$samplingIntervalInSec; sc = abf_current_scale[[units[ch]]]
  d = if (length(a$data) == 1) data.frame(t_s = (seq_len(nrow(a$data[[1]])) - 1) * dt, I_A = a$data[[1]][, ch] * sc) else
    do.call(rbind, lapply(seq_along(a$data), function(k) data.frame(sweep = k, t_s = (seq_len(nrow(a$data[[k]])) - 1) * dt, I_A = a$data[[k]][, ch] * sc)))
  hd = a$header
  mode = c("1" = "event-driven", "2" = "oscilloscope", "3" = "gap-free", "4" = "high-speed oscilloscope", "5" = "episodic")[as.character(hd$nOperationMode)]
  st = if (isTRUE(hd$fFileVersionNumber >= 2)) tryCatch(abf2_settings(path, units[ch]), error = function(e) NULL)
  E = if (!is.null(st) && identical(st$holding_units, "mV")) st$holding / 1000 else if (!is.null(st) && identical(st$holding_units, "V")) st$holding else NA_real_
  date = as.character(hd$uFileStartDate); ms = hd$uFileStartTimeMS
  start = if (!is.null(ms) && nchar(date) == 8) sprintf("%s-%s-%s %02d:%02d:%02d", substr(date, 1, 4), substr(date, 5, 6), substr(date, 7, 8),
                                                         ms %/% 3600000, (ms %/% 60000) %% 60, (ms %/% 1000) %% 60) else NA_character_
  list(file = name, kind = "ca", data = d,
       params = list(technique = paste("Chronoamperometry,", if (is.na(mode)) "ABF" else mode, "recording (Axon ABF)"),
                     instrument = "Axon pCLAMP", scan_rate_V_s = NA_real_, E_V = E,
                     sampling_Hz = round(1 / dt, 6), sweeps = length(a$data), channel = a$channelNames[ch], channel_units = units[ch],
                     filter_Hz = if (!is.null(st)) st$filter_Hz else NA_real_, protocol = if (!is.null(st)) st$protocol else NA_character_,
                     recorded = start))
}
