## =============================================================================
## MethaneSignalProcessor (MSP) -- R port
## -----------------------------------------------------------------------------
## Port of https://github.com/ACBonet/MethaneSignalProcessor (main.py + inc/functions.py)
## Cardona, A., Butturini, A., & Fonollosa, J. (2026). MethaneSignalProcessor (MSP):
## Automated discrimination of diffusive and ebullitive methane fluxes at the water–air
## interface from time-series data. Ecological Informatics, 95, 103781.
## https://doi.org/10.1016/j.ecoinf.2026.103781
##
## Aim: reproduce the Python outputs (scipy / pandas) as closely as possible so
## MSP results can be compared with goAquaFlux on the same chamber data
## (compare_goAquaFlux_FluxSeparator_MSP.R). The port was checked against the
## Python implementation with check_msp_port_parity.R.
##
## Dependencies: base R (processing) + ggplot2 (plots only).
##   The 'signal' package is deliberately NOT used: signal::filtfilt() zero-pads
##   and has no initial conditions, so it does not reproduce scipy.signal.filtfilt.
##   Butterworth design, filtfilt, find_peaks and linregress are re-implemented
##   below following the scipy algorithms.
##
## Python quirks that are reproduced on purpose (for parity) are flagged
## "## PY-QUIRK:". Additions that do not exist in the Python tool, and do not
## affect its results, are flagged "## R-ONLY:".
## Indices that mirror Python are kept 0-based (variables ending in "0").
##
## Usage
##   * RStudio: edit CONFIG and Source the file.
##   * Terminal: Rscript R/MethaneSignalProcessor.R --dir "path/to/folder" [--file name] [--window_peaks 5]
##     where "path/to/folder" contains a "Raw data/" sub-folder with MSP input
##     files (tab-separated .txt, as for the Python tool); results are written
##     to "path/to/folder/Processed data/".
##   * From your own code (e.g. goFlux data already in memory):
##       options(msp.source_only = TRUE); source("MethaneSignalProcessor.R")
##       res <- msp_process_signal(time_s, ch4_ppm, temp_C, pres_mmHg, cfg = CONFIG)
## =============================================================================

## ---- CONFIG -----------------------------------------------------------------
CONFIG <- list(
  base_dir     = file.path("results", "msp_parity", "MSP_R"),   # folder that contains "Raw data/" (Python: script folder)
  raw_dir      = NULL,          # NULL -> <base_dir>/Raw data
  output_dir   = NULL,          # NULL -> <base_dir>/Processed data
  file_name    = NULL,          # NULL -> all .txt in raw_dir; or e.g. "data" (no extension)
  window_peaks = 5,             # half-window (samples) for ebullition peak analysis

  # Input format (Python: pd.read_csv(sep='\t', header=1) -> first line ignored)
  skip_lines   = 1,
  col_time     = "time(s)",
  col_ch4      = "CH4(ppm)",
  col_temp     = c("Temp", "temp"),       # first match is used
  col_pres     = "Pressure(Hg_mm)",
  temp_C_fallback        = NULL,          ## R-ONLY: used only if no temperature column (Python crashes)
  pressure_mmHg_fallback = NULL,          ## R-ONLY: used only if no pressure column (Python crashes)

  # Chamber geometry (hard-coded in Python process_file)
  volume_m3    = 57.5*0.001,
  area_m2      = 1.436544,

  # Processing constants (hard-coded in Python)
  bandpass_hz  = c(0.01, 0.15), # Butterworth band-pass corners (Hz)
  butter_order = 4,
  seg_window   = 10,            # samples excluded on each side of a peak for diffusive fits
  r2_min       = 0.7,
  only_positive = TRUE,

  make_plots   = TRUE,
  write_segments_csv = TRUE     ## R-ONLY: table of all diffusive segments (for goFlux comparison)
)

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- Python-semantics helpers ---------------------------------------------------

# R indices for the Python slice x[start:end] (0-based, end exclusive, negatives wrap)
py_slice <- function(n, start, end) {
  if (start < 0) start <- max(0, n + start)
  if (end < 0)   end   <- max(0, n + end)
  start <- min(start, n); end <- min(end, n)
  if (end <= start) integer(0) else seq.int(start + 1L, end)
}

# R index for the Python element x[i0] (negative i0 wraps)
py_at <- function(n, i0) {
  if (i0 < 0) i0 <- n + i0
  if (i0 < 0 || i0 >= n) stop("Python index ", i0, " out of range for length ", n)
  i0 + 1L
}

pop_sd <- function(x) sqrt(mean((x - mean(x))^2))   # np.std (ddof = 0)

# pandas: Series.rolling(window = w, center = TRUE, min_periods = 1).mean()
msp_rolling_mean <- function(x, w) {
  w <- as.integer(w)
  if (is.na(w) || w < 1)
    stop("Rolling window < 1. MSP needs at least 100 samples per file ",
         "(windows are n/50, n/25 and n/100); the Python version fails too.")
  n <- length(x)
  left  <- w %/% 2L
  right <- w - left - 1L
  ok  <- !is.na(x)
  cs  <- c(0, cumsum(ifelse(ok, x, 0)))
  cnt <- c(0, cumsum(ok))
  i  <- seq_len(n)
  lo <- pmax(1L, i - left)
  hi <- pmin(n, i + right)
  k  <- cnt[hi + 1] - cnt[lo]
  out <- (cs[hi + 1] - cs[lo]) / k
  out[k == 0] <- NA_real_
  out
}

# Python fill_nan_with_local_mean()
msp_fill_nan_local_mean <- function(x, window = 5) {
  filled <- x
  n <- length(x)
  for (i in which(is.na(x))) {
    left  <- if (i > 1) x[max(1, i - window):(i - 1)] else numeric(0)
    right <- if (i < n) x[(i + 1):min(n, i + window)] else numeric(0)
    nb <- c(left, right)
    nb <- nb[!is.na(nb)]
    if (length(nb)) filled[i] <- mean(nb)
  }
  filled
}

## ---- scipy.signal re-implementations ----------------------------------------------

poly_from_roots <- function(r) {            # numpy.poly
  a <- 1 + 0i
  for (ri in r) a <- c(a, 0) - c(0, ri * a)
  a
}

# scipy.signal.butter(N, Wn = c(low, high), btype = "band"), Wn normalised to Nyquist
msp_butter_bandpass <- function(order, low, high) {
  fs <- 2
  warped <- 2 * fs * tan(pi * c(low, high) / fs)
  m  <- seq(-order + 1, order - 1, by = 2)
  p  <- -exp(1i * pi * m / (2 * order))          # buttap poles, k = 1, no zeros
  bw <- warped[2] - warped[1]
  wo <- sqrt(warped[1] * warped[2])
  p_lp <- p * bw / 2                              # lp2bp_zpk
  p_bp <- c(p_lp + sqrt(p_lp^2 - wo^2), p_lp - sqrt(p_lp^2 - wo^2))
  z_bp <- rep(0 + 0i, order)
  k_bp <- bw^order
  fs2 <- 2 * fs                                   # bilinear_zpk
  z_d <- c((fs2 + z_bp) / (fs2 - z_bp), rep(-1 + 0i, length(p_bp) - length(z_bp)))
  p_d <- (fs2 + p_bp) / (fs2 - p_bp)
  k_d <- Re(k_bp * prod(fs2 - z_bp) / prod(fs2 - p_bp))
  list(b = Re(k_d * poly_from_roots(z_d)), a = Re(poly_from_roots(p_d)))
}

msp_lfilter_zi <- function(b, a) {              # scipy.signal.lfilter_zi
  b <- b / a[1]; a <- a / a[1]
  n <- max(length(a), length(b))
  a <- c(a, rep(0, n - length(a))); b <- c(b, rep(0, n - length(b)))
  m <- n - 1
  C <- matrix(0, m, m)
  C[1, ] <- -a[-1]
  if (m > 1) C[cbind(2:m, 1:(m - 1))] <- 1
  solve(diag(m) - t(C), b[-1] - a[-1] * b[1])
}

msp_lfilter <- function(b, a, x, zi) {          # direct form II transposed
  ns <- length(zi); z <- zi; y <- numeric(length(x))
  for (i in seq_along(x)) {
    xi <- x[i]
    yi <- b[1] * xi + z[1]
    if (ns > 1) z[1:(ns - 1)] <- b[2:ns] * xi + z[2:ns] - a[2:ns] * yi
    z[ns] <- b[ns + 1] * xi - a[ns + 1] * yi
    y[i] <- yi
  }
  y
}

# scipy.signal.filtfilt(b, a, x) with defaults: padtype = "odd", padlen = 3 * max(len(a), len(b))
msp_filtfilt <- function(b, a, x) {
  b <- b / a[1]; a <- a / a[1]
  padlen <- 3L * max(length(a), length(b))
  n <- length(x)
  if (n <= padlen) stop("Signal too short for filtfilt: need > ", padlen, " samples.")
  ext <- c(2 * x[1] - x[(padlen + 1):2], x, 2 * x[n] - x[(n - 1):(n - padlen)])
  zi <- msp_lfilter_zi(b, a)
  y <- msp_lfilter(b, a, ext, zi * ext[1])
  y <- rev(msp_lfilter(b, a, rev(y), zi * y[length(y)]))
  y[(padlen + 1):(padlen + n)]
}

# scipy.signal.find_peaks(x, height, distance); returns 0-based indices
msp_find_peaks <- function(x, height = NULL, distance = NULL) {
  n <- length(x); peaks0 <- integer(0); i <- 2L
  while (i < n) {                                  # _local_maxima_1d (plateau midpoints)
    if (x[i - 1] < x[i]) {
      ia <- i + 1L
      while (ia < n && x[ia] == x[i]) ia <- ia + 1L
      if (x[ia] < x[i]) {
        peaks0 <- c(peaks0, ((i - 1L) + (ia - 2L)) %/% 2L)
        i <- ia
      }
    }
    i <- i + 1L
  }
  if (!is.null(height)) peaks0 <- peaks0[which(x[peaks0 + 1L] >= height)]
  if (!is.null(distance) && length(peaks0) > 1) {  # _select_by_peak_distance
    distance <- ceiling(distance)
    np_ <- length(peaks0); keep <- rep(TRUE, np_)
    ord <- order(x[peaks0 + 1L], method = "radix")
    for (ii in np_:1) {
      j <- ord[ii]
      if (!keep[j]) next
      k <- j - 1L
      while (k >= 1 && peaks0[j] - peaks0[k] < distance) { keep[k] <- FALSE; k <- k - 1L }
      k <- j + 1L
      while (k <= np_ && peaks0[k] - peaks0[j] < distance) { keep[k] <- FALSE; k <- k + 1L }
    }
    peaks0 <- peaks0[keep]
  }
  peaks0
}

msp_linregress <- function(x, y) {             # scipy.stats.linregress (slope, intercept, r)
  xm <- mean(x); ym <- mean(y)
  ssxm <- mean((x - xm)^2); ssym <- mean((y - ym)^2); ssxym <- mean((x - xm) * (y - ym))
  r <- if (ssxm == 0 || ssym == 0) 0 else max(-1, min(1, ssxym / sqrt(ssxm * ssym)))
  slope <- ssxym / ssxm
  list(slope = slope, intercept = ym - slope * xm, r = r)
}

## ---- MSP functions ----------------------------------------------------------------

msp_ppm_per_s_to_umol_per_m2h <- function(pressure_mmHg, ppm_per_s, volume_m3, temperature_C, area_m2) {
  R  <- 8.314
  mol_total <- (pressure_mmHg * 133.322 * volume_m3) / (R * (temperature_C + 273.15))
  (ppm_per_s / 1e6) * mol_total * 1e6 * 3600 / area_m2
}

# Segments between peaks, 0-based (start, end) as in Python
msp_segments <- function(peaks0, n, window) {
  if (length(peaks0) == 0) return(list(c(0L, n - 1L)))
  seg <- list()
  if (peaks0[1] > window) seg[[length(seg) + 1]] <- c(0L, peaks0[1] - window)
  for (i in seq_len(max(0, length(peaks0) - 1))) {
    s <- peaks0[i] + window; e <- peaks0[i + 1] - window
    if (e > s) seg[[length(seg) + 1]] <- c(s, e)
  }
  last <- peaks0[length(peaks0)]
  if (last + window < n - 1) seg[[length(seg) + 1]] <- c(last + window, n - 1L)
  seg
}

# All segments with regression stats; 'valid' = what Python reports as a diffusive flux
msp_diffusive_segments <- function(time, y, peaks0, temp_C, pres_mmHg, cfg) {
  n <- length(time)
  rows <- lapply(msp_segments(peaks0, n, cfg$seg_window), function(se) {
    idx <- py_slice(n, se[1], se[2])   ## PY-QUIRK: x[start:end] drops the 'end' sample
    if (length(idx) < 2) return(NULL)
    lr <- msp_linregress(time[idx], y[idx])
    avg_t <- mean(temp_C[idx]); avg_p <- mean(pres_mmHg[idx])
    data.frame(start_idx0 = se[1], end_idx0 = se[2],
               t_start_s = time[idx[1]], t_end_s = time[idx[length(idx)]], n_points = length(idx),
               slope_ppm_s = lr$slope, intercept = lr$intercept, r2 = lr$r^2,
               temp_C = avg_t, pressure_mmHg = avg_p,
               flux_umol_m2_h = msp_ppm_per_s_to_umol_per_m2h(avg_p, lr$slope, cfg$volume_m3, avg_t, cfg$area_m2),
               plotted = lr$r^2 > cfg$r2_min,
               valid = lr$r^2 > cfg$r2_min && (!cfg$only_positive || lr$slope > 0))
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

# Python get_describe()
msp_ebullition_summary <- function(ch4, time, peaks0, window) {
  n <- length(ch4)
  total_adj <- 0; adj_is_float <- FALSE; t_bub <- 0
  for (p in peaks0) {
    s <- max(p - window, 0); e <- min(p + window, n - 1)
    adj <- max(ch4[(s + 1):(e + 1)]) - ch4[s + 1]     # .loc[start:end] is inclusive
    if (adj > 0) { total_adj <- total_adj + adj; adj_is_float <- TRUE }
    t_bub <- t_bub + max(time[e + 1] - time[s + 1], 0)
  }
  n_b <- length(peaks0)
  t_h <- t_bub / 3600
  total_conc <- max(ch4)                               ## PY-QUIRK: labelled "Final" but it is max(CH4)
  out <- list(
    list(v = window * 2, type = "int"),
    list(v = round(total_adj, 2), type = if (adj_is_float) "float" else "int"),
    list(v = round(total_conc, 2), type = "float"),
    list(v = if (total_conc > 0) round(total_adj / total_conc * 100, 2) else 0,
         type = if (total_conc > 0) "float" else "int"),
    list(v = n_b, type = "int"),
    list(v = peaks0, type = "list"),
    list(v = round(t_h, 3), type = "float"),
    list(v = if (t_h > 0) round(n_b / t_h, 2) else 0, type = if (t_h > 0) "float" else "int"))
  names(out) <- c("Peak Analysis Interval",
                  "Total Adjusted CH\u2084 Concentration (ppm)",
                  "Final CH\u2084 Concentration (ppm)",
                  "Contribution of boiling to the total (%)",
                  "Number of Bubbles", "Index of Bubbles",
                  "Total Bubble Time (h)", "Bubbles per Hour")
  out
}

# Core algorithm on vectors (no file I/O) -- the part to call on goFlux data
msp_process_signal <- function(time, ch4, temp_C, pres_mmHg, window_peaks = 5, cfg = CONFIG) {
  signal <- suppressWarnings(as.numeric(ch4)); time <- suppressWarnings(as.numeric(time))
  n <- length(signal)
  if (length(time) != n) stop("time and ch4 have different lengths (", length(time), " vs ", n, ").")
  if (anyNA(signal) || anyNA(time))
    stop(sum(is.na(signal) | is.na(time)), " NA/non-numeric value(s) in time or CH4: remove them first.")
  if (length(unique(round(diff(time), 6))) > 1)
    warning("Irregular time step: MSP uses only the first step (", time[2] - time[1], " s) to set the filter.")
  if (n < 100) stop("MSP needs at least 100 samples (got ", n, ").")
  if (length(temp_C) == 1)    temp_C    <- rep(temp_C, n)
  if (length(pres_mmHg) == 1) pres_mmHg <- rep(pres_mmHg, n)

  # 1. Butterworth band-pass + zero-phase filtering
  fs  <- 1 / (time[2] - time[1])                        ## PY-QUIRK: fs from the first time step only
  nyq <- 0.5 * fs
  bf  <- msp_butter_bandpass(cfg$butter_order, cfg$bandpass_hz[1] / nyq, cfg$bandpass_hz[2] / nyq)
  signal_filtered <- msp_filtfilt(bf$b, bf$a, signal)

  # 2. Smoothing and adaptive peak detection
  window <- n %/% 50
  signal_smoothed <- msp_rolling_mean(signal_filtered, window)
  all_peaks0 <- msp_find_peaks(signal_smoothed)
  pv <- signal_smoothed[all_peaks0 + 1]
  threshold <- mean(pv) + pop_sd(pv) / 3
  above0 <- all_peaks0[which(signal_smoothed[all_peaks0 + 1] > threshold)]
  min_distance <- 1
  if (length(above0) > 1) {
    iv <- diff(above0)
    min_distance <- max(1, trunc(mean(iv) - pop_sd(iv)))
  }
  valid0 <- msp_find_peaks(signal_smoothed, height = threshold, distance = min_distance)

  # 3. Peak correction (manual and automatic branches are identical in Python)
  correct <- function() {
    s <- signal
    for (i in seq_len(max(0, length(valid0) - 1))) {   ## PY-QUIRK: last peak never corrected
      idx <- py_slice(n, valid0[i] - window, n)          ## PY-QUIRK: negative start wraps to the end
      s[idx] <- msp_rolling_mean(s[idx], window) - signal[py_at(n, valid0[i])] +
        signal[py_at(n, valid0[i] - window)]
    }
    s
  }
  corrected <- correct()

  ms <- msp_fill_nan_local_mean(msp_rolling_mean(corrected, n %/% 50))
  ms <- msp_fill_nan_local_mean(msp_rolling_mean(ms, n %/% 25))
  ms <- msp_fill_nan_local_mean(msp_rolling_mean(ms, n %/% 50))
  ms <- msp_fill_nan_local_mean(msp_rolling_mean(ms, n %/% 100))
  ms[seq_len(window)] <- msp_rolling_mean(signal[seq_len(window)], 2)

  auto <- msp_fill_nan_local_mean((signal + corrected) / 2)
  auto <- msp_rolling_mean(auto, n %/% 25)

  final_signal <- (ms + auto) / 2

  # 4. Diffusive fluxes (seg_window = 10 fixed) and ebullition summary (window_peaks)
  segs <- msp_diffusive_segments(time, final_signal, valid0, temp_C, pres_mmHg, cfg)
  ebul <- msp_ebullition_summary(signal, time, valid0, window_peaks)

  # Step-like representation
  step <- numeric(n)
  for (i in seq_len(max(0, length(valid0) - 1)))
    step[py_slice(n, valid0[i], valid0[i + 1])] <- signal[valid0[i] + 1]
  if (length(valid0)) step[py_slice(n, valid0[length(valid0)], n)] <- signal[valid0[length(valid0)] + 1]

  list(final_signal = final_signal, peaks_idx0 = valid0, segments = segs, ebullition = ebul,
       step_signal = step, threshold = threshold, min_distance = min_distance,
       filtered = signal_filtered, smoothed = signal_smoothed, butter = bf)
}

## ---- Output formatting (identical text layout to Python) ------------------------------

py_float_str <- function(x) {                  # Python str(float) for rounded values
  if (is.na(x)) return("nan")
  if (x == trunc(x) && abs(x) < 1e16) return(sprintf("%.1f", x))
  format(x, digits = 15, scientific = FALSE)
}

py_value_str <- function(item) {
  switch(item$type,
         int   = format(as.integer(round(item$v))),
         float = py_float_str(item$v),
         list  = paste0("[", paste(item$v, collapse = ", "), "]"))
}

# str(pd.Series(fluxes).describe().round(2)) with float_format '{:,.2f}'
pandas_describe_str <- function(x) {
  x <- x[!is.na(x)]
  cnt <- length(x)
  vals <- if (cnt == 0) c(0, rep(NaN, 7)) else
    c(cnt, mean(x), if (cnt > 1) sd(x) else NaN, min(x),
      stats::quantile(x, c(.25, .5, .75), names = FALSE, type = 7), max(x))
  vals <- round(vals, 2)
  s <- ifelse(is.nan(vals), "NaN", formatC(vals, format = "f", digits = 2, big.mark = ","))
  lab <- c("count", "mean", "std", "min", "25%", "50%", "75%", "max")
  w <- max(nchar(s)) + 3
  paste(c(sprintf("%-5s%*s", lab, w, s), "dtype: float64"), collapse = "\n")
}

msp_write_results_txt <- function(path, base_name, res) {
  seg <- res$segments
  v <- if (nrow(seg)) seg[seg$valid, ] else seg
  lines <- if (nrow(v)) sprintf(
    "- Slope: %.4f ppm/s | r\u00b2: %.3f | T: %.1f\u00b0C | P: %.1f mmHg | Diffusive Flux: %.2f \u00b5mol/m\u00b2\u00b7h",
    v$slope_ppm_s, v$r2, v$temp_C, v$pressure_mmHg, v$flux_umol_m2_h) else character(0)
  out <- c(paste0("# Source File: ", base_name), "",
           "--- Diffusive Flux Segments ---", lines, "",
           "--- Summary Statistics of Diffusive Fluxes (\u00b5mol/m\u00b2\u00b7h) ---",
           pandas_describe_str(if (nrow(v)) v$flux_umol_m2_h else numeric(0)), "",
           "--- Summary of Ebullitive Events ---",
           paste0(names(res$ebullition), ": ", vapply(res$ebullition, py_value_str, "")))
  ## FIX: explicit UTF-8 (Python's default cp1252 on Windows cannot encode "CH\u2084")
  con <- file(path, open = "w")
  on.exit(close(con))
  writeLines(enc2utf8(out), con, useBytes = TRUE)
}

## ---- Plots (ggplot2) ------------------------------------------------------------

msp_plots <- function(time, signal, res, base_name, plot_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not installed: plots skipped."); return(invisible(NULL))
  }
  library(ggplot2)
  ylab <- "CH\u2084 (ppm)"
  p0 <- res$peaks_idx0 + 1
  dirs <- file.path(plot_dir, c("with_peaks", "slopes", "steps"))
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  d1 <- rbind(data.frame(t = time, y = signal, what = "Original Signal"),
              data.frame(t = time, y = res$final_signal, what = "Processed Signal"))
  g1 <- ggplot(d1, aes(t, y, colour = what)) + geom_line() +
    geom_point(data = data.frame(t = time[p0], y = signal[p0]), aes(t, y),
               inherit.aes = FALSE, colour = "red", size = 2) +
    scale_colour_manual(values = c("#1f77b4", "#ff7f0e"), name = NULL) +
    labs(x = "Time (s)", y = ylab, title = "Original vs Processed Signals with Peaks") +
    theme_bw() + theme(legend.position = "top")
  ggsave(file.path(dirs[1], paste0(base_name, "_peaks_comparison.png")), g1, width = 12, height = 5, dpi = 100)

  seg <- res$segments
  sp <- if (nrow(seg)) seg[seg$plotted, ] else seg      # Python plots any r2 > 0.7 (any sign)
  g2 <- ggplot(data.frame(t = time, y = res$final_signal), aes(t, y)) +
    geom_line(colour = "blue", alpha = .6) +
    labs(x = "Time (s)", y = "CH\u2084 Concentration (ppm)", title = "Slopes & r\u00b2 of Signal Segments",
         subtitle = "Processed signal (CH4_final)") + theme_bw()
  if (nrow(sp)) {
    fit <- do.call(rbind, lapply(seq_len(nrow(sp)), function(i) {
      tt <- c(sp$t_start_s[i], sp$t_end_s[i])
      data.frame(id = i, t = tt, y = sp$slope_ppm_s[i] * tt + sp$intercept[i])
    }))
    rng <- diff(range(res$final_signal))
    lab <- data.frame(t = (sp$t_start_s + sp$t_end_s) / 2)
    lab$y <- sp$slope_ppm_s * lab$t + sp$intercept + 0.02 * rng
    lab$txt <- sprintf("Slope: %.2f\nr\u00b2: %.2f", sp$slope_ppm_s, sp$r2)
    g2 <- g2 + geom_line(data = fit, aes(t, y, group = id), colour = "red", linewidth = 1) +
      geom_label(data = lab, aes(t, y, label = txt), size = 2.8, colour = "darkred",
                 label.size = 0, fill = "white", alpha = .7, vjust = 0)
  }
  ggsave(file.path(dirs[2], paste0(base_name, "_slopes_on_signal.png")), g2, width = 12, height = 6, dpi = 100)

  g3 <- ggplot(data.frame(t = time, y = res$step_signal), aes(t, y)) +
    geom_line(colour = "orange") +
    labs(x = "Time (s)", y = ylab, title = "Analog-like response of valid peaks") + theme_bw()
  ggsave(file.path(dirs[3], paste0(base_name, "_peak_steps.png")), g3, width = 10, height = 5, dpi = 100)
  invisible(NULL)
}

## ---- File wrapper (equivalent of Python process_file) ---------------------------------

msp_read_input <- function(filepath, cfg) {
  df <- tryCatch(
    utils::read.delim(filepath, skip = cfg$skip_lines, check.names = FALSE,
                      stringsAsFactors = FALSE, strip.white = TRUE),
    error = function(e) stop("Error loading ", filepath, ": ", conditionMessage(e), call. = FALSE))
  miss <- setdiff(c(cfg$col_time, cfg$col_ch4), names(df))
  if (length(miss))
    stop(basename(filepath), ": missing column(s) ", paste(shQuote(miss), collapse = ", "),
         ". Found: ", paste(shQuote(names(df)), collapse = ", "),
         ". Check CONFIG$col_* and CONFIG$skip_lines (Python ignores the first line).", call. = FALSE)
  df[[cfg$col_time]] <- suppressWarnings(as.numeric(df[[cfg$col_time]]))
  df[[cfg$col_ch4]]  <- suppressWarnings(as.numeric(df[[cfg$col_ch4]]))
  keep <- !is.na(df[[cfg$col_time]]) & !is.na(df[[cfg$col_ch4]])
  if (any(!keep)) message(basename(filepath), ": dropped ", sum(!keep), " row(s) with missing time/CH4.",
                          " Note: Python uses row labels after dropna, so its results can differ here.")
  df[keep, , drop = FALSE]
}

msp_get_aux <- function(df, cols, fallback, what, fname) {
  hit <- intersect(cols, names(df))
  if (length(hit)) return(as.numeric(df[[hit[1]]]))
  if (is.null(fallback))
    stop(fname, ": no ", what, " column (", paste(cols, collapse = "/"), "). Python requires it too. ",
         "Set CONFIG$", if (what == "temperature") "temp_C_fallback" else "pressure_mmHg_fallback",
         " to use a constant.", call. = FALSE)
  warning(fname, ": no ", what, " column, using constant ", fallback, call. = FALSE)
  rep(fallback, nrow(df))
}

msp_process_file <- function(filepath, output_dir, cfg = CONFIG, window_peaks = cfg$window_peaks) {
  base_name <- tools::file_path_sans_ext(basename(filepath))
  df <- msp_read_input(filepath, cfg)
  temp <- msp_get_aux(df, cfg$col_temp, cfg$temp_C_fallback, "temperature", base_name)
  pres <- msp_get_aux(df, cfg$col_pres, cfg$pressure_mmHg_fallback, "pressure", base_name)

  res <- msp_process_signal(df[[cfg$col_time]], df[[cfg$col_ch4]], temp, pres,
                            window_peaks = window_peaks, cfg = cfg)
  df[["CH4_final (ppm)"]] <- res$final_signal

  dirs <- file.path(output_dir, c("data", "plots", "results"))
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  utils::write.csv(df, file.path(dirs[1], paste0(base_name, "_processed.csv")), row.names = FALSE)
  msp_write_results_txt(file.path(dirs[3], paste0(base_name, "_results.txt")), base_name, res)
  if (isTRUE(cfg$write_segments_csv) && nrow(res$segments))       ## R-ONLY
    utils::write.csv(res$segments, file.path(dirs[3], paste0(base_name, "_segments.csv")), row.names = FALSE)
  if (isTRUE(cfg$make_plots)) msp_plots(df[[cfg$col_time]], df[[cfg$col_ch4]], res, base_name, dirs[2])

  message("Processed and saved: ", base_name)
  invisible(c(list(file = filepath, data = df), res))
}

## ---- Main --------------------------------------------------------------------------

msp_parse_cli <- function(cfg, args = commandArgs(trailingOnly = TRUE)) {
  get <- function(flag) { i <- match(flag, args); if (!is.na(i) && i < length(args)) args[i + 1] else NULL }
  cfg$base_dir     <- get("--dir") %||% cfg$base_dir
  cfg$file_name    <- get("--file") %||% cfg$file_name
  wp <- get("--window_peaks"); if (!is.null(wp)) cfg$window_peaks <- as.integer(wp)
  cfg
}

msp_main <- function(cfg = CONFIG) {
  raw_dir    <- cfg$raw_dir    %||% file.path(cfg$base_dir, "Raw data")
  output_dir <- cfg$output_dir %||% file.path(cfg$base_dir, "Processed data")
  if (!dir.exists(raw_dir)) stop("Raw data directory not found: ", normalizePath(raw_dir, mustWork = FALSE))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  files <- if (is.null(cfg$file_name)) {
    list.files(raw_dir, pattern = "\\.txt$", full.names = TRUE)
  } else file.path(raw_dir, paste0(sub("\\.txt$", "", cfg$file_name), ".txt"))
  files <- files[file.exists(files)]
  if (!length(files)) stop("No .txt input files found in ", raw_dir)

  out <- lapply(files, function(f) tryCatch(
    msp_process_file(f, output_dir, cfg),
    error = function(e) { message("FAILED ", basename(f), ": ", conditionMessage(e)); NULL }))
  names(out) <- tools::file_path_sans_ext(basename(files))
  invisible(out)
}

# options(msp.source_only = TRUE) before source() = define functions only, run nothing
if (!isTRUE(getOption("msp.source_only"))) msp_results <- msp_main(msp_parse_cli(CONFIG))
