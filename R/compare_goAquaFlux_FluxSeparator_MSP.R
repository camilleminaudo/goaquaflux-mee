###############################################################################
## compare_goAquaFlux_FluxSeparator_MSP.R
##
## Objective comparison of three approaches that separate diffusive and
## ebullitive CH4 fluxes in chamber incubations:
##
##   goAquaFlux     goFlux extension (find.bubbles + goFlux model selection)
##   FluxSeparator  Sø et al. (2024, JGR-B) R package, v2.0.0
##   MSP            MethaneSignalProcessor, Cardona et al. (2026, Ecol. Inform.)
##                  https://doi.org/10.1016/j.ecoinf.2026.103781 (R port)
##
## Plus two references that do not separate pathways:
##   goFlux_noSep   goFlux + best.flux on the whole incubation
##   endpoint       (C_end - C_start) / duration: model-free total flux
##
## PART A  Capability table (what each approach can do)            -> CSV
## PART B  Synthetic incubations with known truth                    -> CSV, PDF
##         flux accuracy (total, diffusive, ebullitive) and bubble
##         identification (number, timing, magnitude)
## PART C  Real incubations (goFlux-formatted data, autoID)          -> CSV, PDF
##         agreement between methods, total flux vs endpoint
##
## All fluxes are expressed in nmol m-2 s-1 using the same goFlux flux.term
## (Vtot, Pcham, Area, Tcham, H2O) for every method, so that differences
## reflect the algorithms and not the unit conversion. For FluxSeparator and
## MSP, concentration rates (ppm or ppb per s) are taken from their outputs
## and multiplied by that flux.term (MSP's own conversion uses a hard-coded
## chamber geometry; FluxSeparator's ppm_to_umol is not used).
##
## Where a method does not report a quantity, it is derived and flagged here:
##   FluxSeparator total  = diffusive + ebullitive (both are FS outputs)
##   MSP diffusive        = mean of MSP's valid segment fluxes (as in its
##                          summary statistics)
##   MSP ebullitive       = sum of positive per-peak "adjusted" CH4 rises
##                          (MSP's own definition) / incubation duration
##   MSP total            = diffusive + ebullitive
## Bubble times: goAquaFlux t.bubble; MSP peak time; FluxSeparator time of
## the largest increment inside each detected event (FS gives no time).
##
## Dependencies: dplyr, tidyr, purrr, ggplot2, zoo, TTR, broom, lubridate
## + goFlux (with goAquaFlux), FluxSeparator (installed or its R/ folder),
## MethaneSignalProcessor.R (R port).
###############################################################################


## ============================================================================
## 0. CONFIGURATION
## ============================================================================


repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(repo_root)

CFG <- list(

  run_synthetic = FALSE,
  run_real      = TRUE,
  recompute     = TRUE,     # FALSE: reload saved results and only redo stats/figures

  out_dir = "C:/Projects/myGit/goaquaflux-mee/results/method_comparison_real_incubations",

  ## --- code sources ---------------------------------------------------------
  goflux_dir    = "C:/Projects/myGit/goFlux",   # devtools::load_all(); NULL -> library(goFlux)
  fluxsep_r_dir = NULL,     # path to FluxSeparator/R (v2.0.0); NULL -> library(FluxSeparator)
  msp_file      = "C:/Projects/myGit/goaquaflux-mee/R/MethaneSignalProcessor.R",   # R port of MSP

  ## --- methods to run -------------------------------------------------------
  methods = c("goAquaFlux", "goFlux_noSep", #"goAquaFlux_w30",
              "FS_default","MSP"),


  ## goAquaFlux variants: name -> extra arguments passed to goAquaFlux()
  goaquaflux_variants = list(
    goAquaFlux     = list(),                          # package defaults
    goAquaFlux_w30 = list(bubble.window.size = 30)),  # find.bubbles default window

  ## FluxSeparator "as shipped": data passed in ppm, package default cutoffs
  fs_default = list(runvar_cutoff = 0.5, IndexSpan = 30,
                    concentration_diffusion_cutoff = 1, top_selection = "last",
                    remove_observations_prior = 200,
                    number_of_observations_used = 400,
                    number_of_observations_required = 50),
  ## FluxSeparator tuned from the data of each incubation:
  ## noise sigma = MAD(first differences)/sqrt(2), trend s = median increment
  ##   runvar_cutoff                  = 2.5 s^2 + noise_mult * sigma^2
  ##   concentration_diffusion_cutoff = s * 2 * IndexSpan + mag_noise_mult * sigma
  ## and a diffusive window adapted to short chamber incubations.
  fs_tuned = list(noise_mult = 6, mag_noise_mult = 5, IndexSpan = 30,
                  top_selection = "last",
                  remove_observations_prior = 1,
                  number_of_observations_used = Inf,
                  number_of_observations_required = 30),

  msp_window_peaks = 5,     # MSP default half-window for peak magnitude

  ## --- evaluation ------------------------------------------------------------
  match_tol_s      = 15,    # max |t_est - t_true| to match two bubbles
  endpoint_window_s = 10,   # s averaged at start and end for the endpoint flux

  ## --- synthetic incubations ---------------------------------------------------
  syn = list(
    n_rep = 30, seed = 20260928,
    duration = 600, dt = 1, C0 = 2000,
    ramp = 3,         # s from bubble onset to transient peak at the analyser
    tau = 15,         # s, decay of the overshoot
    edge = 60, min_spacing = 20, size_sdlog = 0.5, ar1 = 0,
    noise_sd    = c(2, 10, 50),         # ppb
    slope       = c(0.1, 1),            # ppb s-1 (diffusive accumulation)
    bubble_mean = c(50, 300, 1500),     # ppb, mean settled step
    n_bubbles   = c(0, 1, 4, 10),
    overshoot   = c(0, 1),              # transient peak excess / settled step
    chamber = list(Vtot = 20, Area = 1000, Pcham = 101.325, Tcham = 20)),  # L, cm2, kPa, degC

  ## --- real incubations (same loading logic as the goFlux vs goAquaFlux script)
  real = list(
    data_file  = "C:/Projects/myGit/goaquaflux-mee/data/many_incubations.RData",
    data_object = "mydata_all", aux_object = "myauxfile",
    aux_duration_col = "duration",   # copied to obs.length for autoID()
    shoulder = 0,
    time_tolerance_s = 5,            # skip if aux start.time precedes data by more
    ids = NULL,                      # optional subset of UniqueIDs
    n_max = NULL,                    # optional random subset size (quick tests)
    prec_default = NULL,             # CH4 precision (ppb) if no CH4_prec column
    reference_method = "goAquaFlux", # reference for between-method agreement
    n_example_plots = 12)
)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a


## ============================================================================
## 1. LOAD CODE
## ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(ggplot2)
  library(zoo); library(TTR); library(broom); library(lubridate)
})

## goFlux / goAquaFlux

devtools::load_all("C:/Projects/myGit/goFlux")   # path to the package root
stopifnot(exists("goAquaFlux"), exists("goFlux"), exists("best.flux"),
          exists("flux.term"))

## FluxSeparator: exported functions in an environment FS
FS <- new.env()
if (!is.null(CFG$fluxsep_r_dir)) {
  for (f in c("utils-bubble-detection.R", "utils-smoothing.R",
              "ebullitive_flux.R", "diffusive_flux.R")) {
    p <- file.path(CFG$fluxsep_r_dir, f)
    if (file.exists(p)) sys.source(p, envir = FS)
  }
} else {
  FS$ebullitive_flux <- FluxSeparator::ebullitive_flux
  FS$diffusive_flux  <- FluxSeparator::diffusive_flux
}
stopifnot(exists("ebullitive_flux", envir = FS), exists("diffusive_flux", envir = FS))

## MSP (functions only, nothing is run)
MSP <- new.env()
options(msp.source_only = TRUE)
source(CFG$msp_file, local = MSP)
stopifnot(exists("msp_process_signal", envir = MSP))

dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)


## ============================================================================
## PART A. CAPABILITY TABLE
## ============================================================================

capabilities <- tribble(
  ~Feature, ~goAquaFlux, ~FluxSeparator, ~MSP,
  "Implementation",
  "R (goFlux extension)", "R package; Shiny app in v2.0.0", "Python (R port used here)",
  "Designed for",
  "Discrete floating-chamber incubations", "Automated DIY chambers, long pump cycles", "Chamber CH4 time series (one file per incubation)",
  "Gases",
  "Any goFlux gas; bubbles detected on CH4, other gases cut only at abrupt slope change", "Any column; ebullition meant for CH4, CO2 diffusive with look_for_bubbles = FALSE", "CH4 only",
  "Sampling requirements",
  ">= 30 obs; irregular sampling OK (interpolated for detection)", "Regular ~1 Hz implied (windows in samples, gaps in s); > 100 obs per cycle", "Regular (fs from first step); >= 100 samples; fs > 0.3 Hz (band-pass 0.01-0.15 Hz)",
  "Detection statistic",
  "Rolling variance of increments (or of concentration) on robust-standardised series", "Running variance of concentration (5 obs)", "Butterworth band-pass + rolling mean + peak finding",
  "Detection threshold",
  "Adaptive per incubation (quantile / median + k MAD) + max/median ratio guard", "Absolute, user-set (default 0.5 conc. units^2)", "Relative per incubation (mean + SD/3 of local maxima): always >= 1 peak if the smoothed signal has a local maximum",
  "Bubble magnitude",
  "Step regression with shared slope, optional overshoot / re-equilibration term", "Last (or max) minus min in a window padded by +/-30 obs", "Max within +/-5 samples of the peak minus value 5 samples before",
  "Bubble timing reported",
  "Yes (event start/end, step time, peak time)", "No (only per-cycle sums and counts)", "Yes (peak index)",
  "Diffusive flux",
  "goFlux LM/HM model selection before first bubble (bubble gas)", "Linear regression on first bubble-free segment (skip 200 obs, use 400); optional HMR", "Linear regressions on inter-peak segments of a smoothed signal; kept if r2 > 0.7 and slope > 0; several values per incubation",
  "Ebullitive flux",
  "Sum of steps / incubation time x flux term", "Sum of rises per hour (conc. units h-1)", "Not as a flux: total 'adjusted' CH4 rise (ppm), bubbles per hour, % of max CH4",
  "Total flux",
  "Diffusive + ebullitive, checked against endpoint estimate", "Not reported (sum of the two outputs)", "Not reported",
  "Unit conversion",
  "flux.term (V, P, T, A, H2O dilution)", "ppm_to_umol (user supplies V, A, T)", "Ideal gas with hard-coded chamber V and A; P, T from file",
  "Uncertainty",
  "SE for diffusive, each bubble, ebullitive and total", "SE of diffusive slope only", "None",
  "Quality flags",
  "MDF, model diagnostics, inconsistency / suspicious-total flags", "r2, p-value of diffusive fit", "r2 filter on segments",
  "User parameters",
  "Many (detection, magnitude model, diffusive criteria); adaptive defaults", "~4 main (cutoffs, IndexSpan, diffusive window)", "window_peaks; other constants hard-coded",
  "Diagnostics",
  "flux.plot.aqua (fitted event model redrawn)", "ggplot panels; Shiny app with manual selection", "PNG plots per file (peaks, segments, steps)"
)
write.csv(capabilities, file.path(CFG$out_dir, "A_capabilities.csv"), row.names = FALSE)


## ============================================================================
## 2. SHARED HELPERS
## ============================================================================

## Run an expression silently, keeping its warnings as text
quietly <- function(expr) {
  warns <- character(0)
  val <- NULL
  utils::capture.output(
    val <- withCallingHandlers(suppressMessages(expr), warning = function(w) {
      warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }),
    type = "output")
  list(value = val, warnings = unique(warns))
}

empty_events <- function() data.frame(t_event = numeric(0), t_start = numeric(0),
                                      t_end = numeric(0), magnitude_ppb = numeric(0))

## Put an irregular series on a regular grid (FS and MSP assume regular sampling)
regularize <- function(t, y) {
  dt <- stats::median(diff(t))
  if (all(abs(diff(t) - dt) <= 0.01 * dt))
    return(list(t = t, y = y, dt = dt, regularized = FALSE))
  tg <- seq(0, max(t), by = dt)
  list(t = tg, y = stats::approx(t, y, xout = tg, rule = 2)$y, dt = dt, regularized = TRUE)
}

## goFlux flux.term exactly as goAquaFlux computes it
incubation_flux_term <- function(df) {
  first_ok <- function(x, default = NA) { x <- x[!is.na(x)]; if (length(x)) x[1] else default }
  Vtot <- if ("Vtot" %in% names(df) && any(!is.na(df$Vtot))) first_ok(df$Vtot) else
    first_ok(df$Vcham) + first_ok(df$Area) * first_ok(df$offset) / 1000
  P <- first_ok(df$Pcham, 101.325); Tc <- first_ok(df$Tcham, 15)
  H2O <- if ("H2O_ppm" %in% names(df)) first_ok(df$H2O_ppm, 0) / 1e6 else 0
  list(flux_term = flux.term(Vtot, P, first_ok(df$Area), Tc, H2O), T_C = Tc, P_kPa = P)
}

## One incubation = the goFlux data frame (flag == 1, Etime from 0) + metadata
make_incubation <- function(id, df, prec_default = NULL) {
  df <- as.data.frame(df)
  df <- df[df$flag == 1 & !is.na(df$CH4dry_ppb) & !is.na(df$Etime), , drop = FALSE]
  df <- df[order(df$Etime), , drop = FALSE]
  df <- df[!duplicated(df$Etime), , drop = FALSE]
  if (nrow(df) < 30) stop("fewer than 30 CH4 observations")
  df$Etime <- df$Etime - df$Etime[1]
  df$UniqueID <- as.character(id)
  if (!"H2O_ppm" %in% names(df)) df$H2O_ppm <- 0
  if (!"CH4_prec" %in% names(df) || all(is.na(df$CH4_prec))) {
    if (is.null(prec_default)) stop("no CH4_prec column: set CFG$real$prec_default")
    df$CH4_prec <- prec_default
  }
  ft <- incubation_flux_term(df)
  list(id = as.character(id), df = df, t = df$Etime, ch4 = df$CH4dry_ppb,
       duration = max(df$Etime), flux_term = ft$flux_term,
       T_C = ft$T_C, P_kPa = ft$P_kPa)
}

## Model-free total flux: mean over the last minus the first window
endpoint_total <- function(inc, w = CFG$endpoint_window_s) {
  C0 <- mean(inc$ch4[inc$t <= w]); Cf <- mean(inc$ch4[inc$t >= inc$duration - w])
  (Cf - C0) / inc$duration * inc$flux_term
}

## Time of the largest increment inside [a, b] (event time for FS)
largest_increment_time <- function(t, y, a, b) {
  i <- which(t >= a & t <= b)
  if (length(i) < 2) return((a + b) / 2)
  j <- i[-1][which.max(diff(y[i]))]
  t[j]
}


## ============================================================================
## 3. METHOD RUNNERS
## Each returns list(flux = list(total, diffusive, ebullition),
##                   events = data.frame(t_event, t_start, t_end, magnitude_ppb),
##                   note = character, check = numeric)
## ============================================================================

## --- goAquaFlux ------------------------------------------------------------
run_goAquaFlux <- function(inc, args = list()) {
  out <- quietly(do.call(goAquaFlux, c(list(dataframe = inc$df, gastype = "CH4dry_ppb",
                                            H2O_col = "H2O_ppm"), args)))
  r <- out$value
  s <- r$flux_summary
  b <- r$bubbles
  ev <- if (is.null(b) || nrow(b) == 0) empty_events() else
    data.frame(t_event = if ("t.bubble" %in% names(b)) b$t.bubble else b$start,
               t_start = b$start, t_end = b$end, magnitude_ppb = b$magnitude)
  list(flux = list(total = s$flux_total[1], diffusive = s$flux_diffusive[1],
                   ebullition = s$flux_ebullition[1]),
       events = ev, note = paste(out$warnings, collapse = " | "))
}

## --- goFlux without separation ----------------------------------------------
run_goFlux_noSep <- function(inc) {
  out <- quietly(best.flux(goFlux(inc$df, "CH4dry_ppb", H2O_col = "H2O_ppm")))
  list(flux = list(total = out$value$best.flux[1], diffusive = NA, ebullition = NA),
       events = NULL, note = paste(out$warnings, collapse = " | "))
}

## --- FluxSeparator: event-level re-implementation --------------------------
## Reproduces detect_bubbles() + the classification in ebullitive_flux()
## (v2.0.0) for one cycle; quirks kept. The per-cycle sum is compared with
## ebullitive_flux() for every incubation (column `check`).
fs_events <- function(time, conc, runvar_cutoff, IndexSpan,
                      concentration_diffusion_cutoff, top_selection = "last",
                      runvar_window = 5, time_gap_seconds = 30, min_obs = 100) {
  if (length(conc) <= min_obs) return(empty_events())
  rv <- TTR::runVar(conc, n = runvar_window)
  flagged <- which(rv > runvar_cutoff)
  flagged <- flagged[-1]                         # as in FS (drop_na on lag)
  if (!length(flagged)) return(empty_events())
  ids <- unique(unlist(lapply(flagged, function(x) {
    v <- (x - IndexSpan):(x + IndexSpan); v[v > 0] })))
  ids <- sort(ids[ids <= length(conc)])
  tt <- time[ids]; cc <- conc[ids]
  td <- c(NA, diff(tt)); keep <- !is.na(td)      # as in FS
  tt <- tt[keep]; cc <- cc[keep]; td <- td[keep]
  if (!length(tt)) return(empty_events())
  grp <- 1 + cumsum(td > time_gap_seconds)
  ev <- lapply(split(seq_along(tt), grp), function(i) {
    x <- cc[i]; ti <- tt[i]
    if (!(x[1] < x[length(x)])) return(NULL)
    top <- if (top_selection == "last") x[length(x)] else max(x)
    d <- top - min(x)
    if (!(d > concentration_diffusion_cutoff && ti[which.min(x)] < ti[which.max(x)])) return(NULL)
    data.frame(t_event = largest_increment_time(time, conc, min(ti), max(ti)),
               t_start = min(ti), t_end = max(ti), magnitude_ppb = d)
  })
  ev <- do.call(rbind, ev)
  if (is.null(ev)) empty_events() else ev
}

run_fluxseparator <- function(inc, tuned = FALSE) {
  reg <- regularize(inc$t, inc$ch4)
  if (tuned) {
    a <- CFG$fs_tuned
    x <- reg$y; k <- 1                             # ppb
    d <- diff(x)
    sig <- max(stats::mad(d) / sqrt(2), 1e-6)
    s <- stats::median(d)
    p <- list(runvar_cutoff = 2.5 * s^2 + a$noise_mult * sig^2,
              IndexSpan = a$IndexSpan,
              concentration_diffusion_cutoff = max(s, 0) * 2 * a$IndexSpan + a$mag_noise_mult * sig,
              top_selection = a$top_selection,
              remove_observations_prior = a$remove_observations_prior,
              number_of_observations_used = a$number_of_observations_used,
              number_of_observations_required = a$number_of_observations_required)
  } else {
    p <- CFG$fs_default
    x <- reg$y / 1000; k <- 1000                   # ppm, as the defaults assume
  }
  dfs <- data.frame(datetime = as.POSIXct("2000-01-01", tz = "UTC") + reg$t,
                    PumpCycle = 1, station = inc$id, tempC = inc$T_C, sensor = 1,
                    CH4conc = x)
  notes <- character(0)

  eb <- tryCatch({
    o <- quietly(FS$ebullitive_flux(dfs, concentration_values = "CH4conc",
                                    runvar_cutoff = p$runvar_cutoff, IndexSpan = p$IndexSpan,
                                    concentration_diffusion_cutoff = p$concentration_diffusion_cutoff,
                                    top_selection = p$top_selection, show_plots = FALSE))
    v <- o$value; if (!is.data.frame(v) && is.list(v)) v <- v[[1]]
    as.data.frame(v)
  }, error = function(e) { notes <<- c(notes, paste("ebullitive_flux:", conditionMessage(e))); NULL })

  di <- tryCatch({
    o <- quietly(FS$diffusive_flux(dfs, concentration_values = "CH4conc",
                                   runvar_cutoff = p$runvar_cutoff, IndexSpan = p$IndexSpan,
                                   remove_observations_prior = p$remove_observations_prior,
                                   number_of_observations_used = p$number_of_observations_used,
                                   number_of_observations_required = p$number_of_observations_required,
                                   look_for_bubbles = TRUE, show_plots = FALSE))
    v <- o$value; if (!is.data.frame(v) && is.list(v)) v <- v[[1]]
    as.data.frame(v)
  }, error = function(e) { notes <<- c(notes, paste("diffusive_flux:", conditionMessage(e))); NULL })

  ev <- fs_events(reg$t, x, p$runvar_cutoff, p$IndexSpan,
                  p$concentration_diffusion_cutoff, p$top_selection)
  ev$magnitude_ppb <- ev$magnitude_ppb * k

  ebul <- if (!is.null(eb) && nrow(eb)) eb$concentration_per_time[1] * k / 3600 * inc$flux_term else NA
  if (!is.null(eb) && nrow(eb) && is.na(eb$concentration_per_time[1])) ebul <- 0
  dif  <- if (!is.null(di) && nrow(di)) di$slope_concentration_hr[1] * k / 3600 * inc$flux_term else NA
  if (is.null(di) || !nrow(di)) notes <- c(notes, "no diffusive segment accepted by FS")
  chk <- if (!is.null(eb) && nrow(eb))
    abs(sum(ev$magnitude_ppb) - k * dplyr::coalesce(eb$sum_bubbles_concentration[1], 0)) else NA
  if (reg$regularized) notes <- c(notes, "series regularized")

  list(flux = list(total = dif + ebul, diffusive = dif, ebullition = ebul),
       events = ev, note = paste(notes, collapse = " | "), check = chk)
}

## --- MSP ---------------------------------------------------------------------
run_msp <- function(inc, window_peaks = CFG$msp_window_peaks) {
  reg <- regularize(inc$t, inc$ch4)
  if (1 / reg$dt <= 2 * MSP$CONFIG$bandpass_hz[2])
    stop("sampling too slow for MSP band-pass (needs > ", 2 * MSP$CONFIG$bandpass_hz[2], " Hz)")
  sig_ppm <- reg$y / 1000
  out <- quietly(MSP$msp_process_signal(reg$t, sig_ppm, inc$T_C, inc$P_kPa * 7.50062,
                                        window_peaks = window_peaks, cfg = MSP$CONFIG))
  r <- out$value
  seg <- r$segments
  v <- if (nrow(seg)) seg[seg$valid, , drop = FALSE] else seg
  dif <- if (nrow(v)) mean(v$slope_ppm_s) * 1000 * inc$flux_term else NA

  ## per-peak magnitude exactly as MSP get_describe()
  n <- length(sig_ppm); p0 <- r$peaks_idx0
  ev <- if (!length(p0)) empty_events() else do.call(rbind, lapply(p0, function(p) {
    s <- max(p - window_peaks, 0); e <- min(p + window_peaks, n - 1)
    adj <- max(sig_ppm[(s + 1):(e + 1)]) - sig_ppm[s + 1]
    data.frame(t_event = reg$t[p + 1], t_start = reg$t[s + 1], t_end = reg$t[e + 1],
               magnitude_ppb = max(adj, 0) * 1000)
  }))
  ebul <- sum(ev$magnitude_ppb) / max(reg$t) * inc$flux_term
  notes <- c(out$warnings, if (!nrow(v)) "no valid diffusive segment",
             if (reg$regularized) "series regularized")
  list(flux = list(total = dif + ebul, diffusive = dif, ebullition = ebul),
       events = ev, note = paste(notes, collapse = " | "))
}

## --- registry --------------------------------------------------------------------
RUNNERS <- list()
for (m in CFG$methods) {
  RUNNERS[[m]] <- local({
    mm <- m
    if (mm %in% names(CFG$goaquaflux_variants)) {
      function(inc) run_goAquaFlux(inc, CFG$goaquaflux_variants[[mm]])
    } else switch(mm,
                  goFlux_noSep = function(inc) run_goFlux_noSep(inc),
                  FS_default   = function(inc) run_fluxseparator(inc, tuned = FALSE),
                  # FS_tuned     = function(inc) run_fluxseparator(inc, tuned = TRUE),
                  MSP          = function(inc) run_msp(inc),
                  stop("Unknown method: ", mm))
  })
}
METHOD_LEVELS <- c(CFG$methods, "endpoint")

## Run all methods on one incubation
process_incubation <- function(inc) {
  fl <- list(); ev <- list()
  for (m in names(RUNNERS)) {
    t0 <- proc.time()[["elapsed"]]
    r <- tryCatch(RUNNERS[[m]](inc), error = function(e) list(error = conditionMessage(e)))
    el <- proc.time()[["elapsed"]] - t0
    if (!is.null(r$error)) {
      fl[[m]] <- data.frame(id = inc$id, method = m, total = NA, diffusive = NA,
                            ebullition = NA, n_events = NA, status = "error",
                            note = r$error, check = NA, seconds = el)
      next
    }
    ne <- if (is.null(r$events)) NA else nrow(r$events)
    fl[[m]] <- data.frame(id = inc$id, method = m,
                          total = r$flux$total %||% NA, diffusive = r$flux$diffusive %||% NA,
                          ebullition = r$flux$ebullition %||% NA, n_events = ne,
                          status = "ok", note = r$note %||% "", check = r$check %||% NA,
                          seconds = el)
    if (!is.null(r$events) && nrow(r$events))
      ev[[m]] <- cbind(id = inc$id, method = m, r$events)
  }
  fl$endpoint <- data.frame(id = inc$id, method = "endpoint", total = endpoint_total(inc),
                            diffusive = NA, ebullition = NA, n_events = NA,
                            status = "ok", note = "", check = NA, seconds = 0)
  list(flux = bind_rows(fl), events = bind_rows(ev))
}

run_all <- function(incs, label) {
  cat("\nRunning", length(incs), label, "incubations x", length(RUNNERS), "methods\n")
  pb <- txtProgressBar(0, length(incs), style = 3)
  res <- lapply(seq_along(incs), function(k) {
    setTxtProgressBar(pb, k)
    tryCatch(process_incubation(incs[[k]]),
             error = function(e) { message("\n", incs[[k]]$id, ": ", conditionMessage(e)); NULL })
  })
  close(pb)
  list(flux   = bind_rows(lapply(res, `[[`, "flux")),
       events = bind_rows(lapply(res, `[[`, "events")))
}


## ============================================================================
## 4. EVALUATION HELPERS
## ============================================================================

## Greedy one-to-one matching of event times within a tolerance
match_events <- function(t_ref, t_est, tol = CFG$match_tol_s) {
  if (!length(t_ref) || !length(t_est)) return(matrix(integer(0), 0, 2))
  d <- abs(outer(t_ref, t_est, "-"))
  cand <- which(d <= tol, arr.ind = TRUE)
  if (!nrow(cand)) return(matrix(integer(0), 0, 2))
  cand <- cand[order(d[cand]), , drop = FALSE]
  used_r <- used_e <- integer(0); keep <- logical(nrow(cand))
  for (k in seq_len(nrow(cand))) {
    if (!(cand[k, 1] %in% used_r) && !(cand[k, 2] %in% used_e)) {
      keep[k] <- TRUE; used_r <- c(used_r, cand[k, 1]); used_e <- c(used_e, cand[k, 2])
    }
  }
  cand[keep, , drop = FALSE]
}

## Compare events of every method against reference events (truth or a method)
evaluate_events <- function(events, ref_events, ids, methods) {
  per_inc <- list(); pairs <- list()
  for (i in ids) {
    r <- ref_events[ref_events$id == i, , drop = FALSE]
    for (m in methods) {
      e <- events[events$id == i & events$method == m, , drop = FALSE]
      mt <- match_events(r$t_event, e$t_event)
      per_inc[[length(per_inc) + 1]] <- data.frame(
        id = i, method = m, n_ref = nrow(r), n_est = nrow(e), TP = nrow(mt),
        FP = nrow(e) - nrow(mt), FN = nrow(r) - nrow(mt))
      if (nrow(mt)) pairs[[length(pairs) + 1]] <- data.frame(
        id = i, method = m,
        t_ref = r$t_event[mt[, 1]], t_est = e$t_event[mt[, 2]],
        mag_ref = r$magnitude_ppb[mt[, 1]], mag_est = e$magnitude_ppb[mt[, 2]])
    }
  }
  list(per_inc = bind_rows(per_inc), pairs = bind_rows(pairs))
}

## Lin's concordance correlation coefficient
ccc <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 3) return(NA_real_)
  2 * stats::cov(x, y) / (stats::var(x) + stats::var(y) + (mean(x) - mean(y))^2)
}

## Trace with the events of every method, one panel per method
plot_incubation <- function(inc, events, truth = NULL, title = inc$id) {
  lv <- setdiff(names(RUNNERS), "goFlux_noSep")
  tr <- tidyr::crossing(data.frame(t = inc$t, ch4 = inc$ch4), method = lv)
  tr$method <- factor(tr$method, levels = lv)
  e <- events[events$id == inc$id & events$method %in% lv, , drop = FALSE]
  p <- ggplot(tr, aes(t, ch4))
  if (nrow(e)) {
    e$method <- factor(e$method, levels = lv)
    p <- p + geom_rect(data = e, inherit.aes = FALSE,
                       aes(xmin = t_start, xmax = t_end, ymin = -Inf, ymax = Inf),
                       fill = "steelblue", alpha = 0.2) +
      geom_vline(data = e, aes(xintercept = t_event), colour = "steelblue4", linewidth = 0.4)
  }
  if (!is.null(truth) && nrow(truth))
    p <- p + geom_vline(xintercept = truth$t_event, colour = "red", linetype = 2, linewidth = 0.4)
  p + geom_line(linewidth = 0.3) +
    facet_wrap(~ method, ncol = 1, scales = "free_y") +
    labs(title = title, x = "Etime (s)", y = "CH4 (ppb)",
         subtitle = paste0("blue = detected events (line = event time)",
                           if (!is.null(truth)) "; red dashed = true bubble onsets")) +
    theme_bw(base_size = 9)
}

component_levels <- c("total", "diffusive", "ebullition")


## ============================================================================
## PART B. SYNTHETIC INCUBATIONS
## ============================================================================

bubble_shape <- function(u, M, A, tau, ramp) {
  out <- numeric(length(u))
  if (ramp > 0) { r <- u >= 0 & u < ramp; out[r] <- (M + A) * u[r] / ramp }
  p <- u >= ramp
  out[p] <- M + A * exp(-(u[p] - ramp) / tau)
  out
}

place_bubbles <- function(n, lower, upper, min_spacing, max_try = 5000) {
  if (n == 0) return(numeric(0))
  for (k in seq_len(max_try)) {
    tb <- sort(round(runif(n, lower, upper)))
    if (n == 1 || all(diff(tb) >= min_spacing)) return(tb)
  }
  stop("Could not place ", n, " bubbles")
}

simulate_series <- function(sc, s = CFG$syn) {
  t  <- seq(0, s$duration - s$dt, by = s$dt)
  tb <- place_bubbles(sc$n_bubbles, s$edge, s$duration - s$edge, s$min_spacing)
  M  <- if (sc$n_bubbles > 0) sc$bubble_mean *
    exp(rnorm(sc$n_bubbles, -s$size_sdlog^2 / 2, s$size_sdlog)) else numeric(0)
  y <- s$C0 + sc$slope * t
  for (i in seq_along(tb)) y <- y + bubble_shape(t - tb[i], M[i], sc$overshoot * M[i], s$tau, s$ramp)
  eps <- if (s$ar1 > 0) as.numeric(arima.sim(list(ar = s$ar1), n = length(t),
                                             sd = sc$noise_sd * sqrt(1 - s$ar1^2))) else
                                               rnorm(length(t), 0, sc$noise_sd)
  list(t = t, y = y + eps, tb = tb, M = M)
}

build_synthetic <- function(s = CFG$syn) {
  g1 <- expand.grid(noise_sd = s$noise_sd, slope = s$slope, bubble_mean = s$bubble_mean,
                    n_bubbles = s$n_bubbles[s$n_bubbles > 0], overshoot = s$overshoot)
  g0 <- if (0 %in% s$n_bubbles) expand.grid(noise_sd = s$noise_sd, slope = s$slope,
                                            bubble_mean = NA, n_bubbles = 0, overshoot = NA) else NULL
  sc_grid <- rbind(g1, g0); sc_grid$scenario <- seq_len(nrow(sc_grid))
  jobs <- expand.grid(scenario = sc_grid$scenario, rep = seq_len(s$n_rep))
  incs <- vector("list", nrow(jobs)); truth_flux <- truth_ev <- vector("list", nrow(jobs))
  ch <- s$chamber
  for (j in seq_len(nrow(jobs))) {
    sc <- sc_grid[jobs$scenario[j], ]
    set.seed(s$seed + j)
    sim <- simulate_series(sc, s)
    id <- sprintf("syn_%05d", j)
    df <- data.frame(UniqueID = id, Etime = sim$t, flag = 1, CH4dry_ppb = sim$y,
                     CH4_prec = sc$noise_sd, H2O_ppm = 0, Vtot = ch$Vtot, Area = ch$Area,
                     Pcham = ch$Pcham, Tcham = ch$Tcham)
    inc <- make_incubation(id, df)
    incs[[j]] <- inc
    ft <- inc$flux_term
    truth_flux[[j]] <- cbind(id = id, sc, rep = jobs$rep[j],
                             true_diffusive = sc$slope * ft,
                             true_ebullition = sum(sim$M) / inc$duration * ft,
                             true_total = sc$slope * ft + sum(sim$M) / inc$duration * ft)
    if (length(sim$tb)) truth_ev[[j]] <- data.frame(id = id, t_event = sim$tb,
                                                    magnitude_ppb = sim$M)
  }
  list(incs = incs, truth_flux = bind_rows(truth_flux), truth_ev = bind_rows(truth_ev),
       grid = sc_grid)
}

if (CFG$run_synthetic) {

  f_syn <- file.path(CFG$out_dir, "B_synthetic_results.rds")
  if (CFG$recompute || !file.exists(f_syn)) {
    syn <- build_synthetic()
    t_start <- Sys.time()
    res_syn <- run_all(syn$incs, "synthetic")
    cat("Run time:", format(round(Sys.time() - t_start, 1)), "\n")
    saveRDS(list(syn = syn, res = res_syn), f_syn)
  } else {
    tmp <- readRDS(f_syn); syn <- tmp$syn; res_syn <- tmp$res
  }

  truth_flux <- syn$truth_flux
  truth_ev   <- if (nrow(syn$truth_ev)) syn$truth_ev else
    data.frame(id = character(0), t_event = numeric(0), magnitude_ppb = numeric(0))
  fl <- res_syn$flux %>% mutate(method = factor(method, levels = METHOD_LEVELS))
  ev <- res_syn$events

  ## --- FS fidelity check (event re-implementation vs ebullitive_flux) --------
  fs_chk <- fl %>% filter(grepl("^FS_", method), !is.na(check))
  if (nrow(fs_chk)) cat("\nFluxSeparator event re-implementation vs ebullitive_flux():",
                        "max |difference| =", signif(max(fs_chk$check), 3), "ppb over",
                        nrow(fs_chk), "runs\n")

  ## --- B1. Flux accuracy ---------------------------------------------------------
  flux_long <- fl %>%
    select(id, method, status, total, diffusive, ebullition) %>%
    pivot_longer(c(total, diffusive, ebullition), names_to = "component", values_to = "est") %>%
    left_join(truth_flux %>% select(id, noise_sd, slope, bubble_mean, n_bubbles, overshoot,
                                    true_total, true_diffusive, true_ebullition), by = "id") %>%
    mutate(true = case_when(component == "total" ~ true_total,
                            component == "diffusive" ~ true_diffusive,
                            TRUE ~ true_ebullition),
           component = factor(component, levels = component_levels)) %>%
    select(-true_total, -true_ebullition) %>%
    ## drop components a method does not provide at all
    group_by(method, component) %>% filter(any(!is.na(est))) %>% ungroup() %>%
    mutate(rel_err = ifelse(true > 0, (est - true) / true, NA_real_))

  flux_summary <- flux_long %>%
    group_by(method, component) %>%
    summarise(n = n(), fail_rate = mean(is.na(est)),
              med_rel_err = median(rel_err, na.rm = TRUE),
              med_abs_rel_err = median(abs(rel_err), na.rm = TRUE),
              q90_abs_rel_err = quantile(abs(rel_err), 0.9, na.rm = TRUE),
              .groups = "drop")
  flux_by_factor <- flux_long %>%
    group_by(method, component, noise_sd, n_bubbles, overshoot) %>%
    summarise(n = n(), fail_rate = mean(is.na(est)),
              med_rel_err = median(rel_err, na.rm = TRUE),
              med_abs_rel_err = median(abs(rel_err), na.rm = TRUE), .groups = "drop")
  false_ebul <- flux_long %>%
    filter(n_bubbles == 0, component == "ebullition") %>%
    group_by(method) %>%
    summarise(n = n(), share_nonzero = mean(est > 0, na.rm = TRUE),
              median_false_flux = median(est, na.rm = TRUE),
              median_ratio_to_true_diffusive = median(est / true_diffusive, na.rm = TRUE),
              .groups = "drop")

  ## --- B2. Bubble identification ---------------------------------------------------
  ev_methods <- setdiff(CFG$methods, "goFlux_noSep")
  ev_eval <- evaluate_events(ev, truth_ev, truth_flux$id, ev_methods)
  det_per_inc <- ev_eval$per_inc %>%
    left_join(truth_flux %>% select(id, noise_sd, slope, bubble_mean, n_bubbles, overshoot),
              by = "id") %>%
    mutate(method = factor(method, levels = METHOD_LEVELS))
  det_summary <- det_per_inc %>%
    group_by(method) %>%
    summarise(recall = sum(TP) / sum(n_ref),
              precision = if (sum(n_est) > 0) sum(TP) / sum(n_est) else NA_real_,
              F1 = 2 * sum(TP) / (sum(n_ref) + sum(n_est)),
              count_ratio = sum(n_est) / sum(n_ref),
              false_events_no_bubble = mean(n_est[n_bubbles == 0]),
              share_null_with_events = mean(n_est[n_bubbles == 0] > 0), .groups = "drop")
  det_by_snr <- det_per_inc %>%
    filter(n_bubbles > 0) %>%
    group_by(method, noise_sd, bubble_mean, n_bubbles) %>%
    summarise(recall = sum(TP) / sum(n_ref),
              precision = if (sum(n_est) > 0) sum(TP) / sum(n_est) else NA_real_,
              count_ratio = sum(n_est) / sum(n_ref), .groups = "drop") %>%
    mutate(snr = bubble_mean / noise_sd)
  pairs <- ev_eval$pairs %>%
    left_join(truth_flux %>% select(id, noise_sd, overshoot, n_bubbles), by = "id") %>%
    mutate(method = factor(method, levels = METHOD_LEVELS),
           timing_err = t_est - t_ref, mag_rel_err = (mag_est - mag_ref) / mag_ref)
  pair_summary <- pairs %>%
    group_by(method) %>%
    summarise(n_matched = n(),
              med_timing_err_s = median(timing_err), mad_timing_err_s = mad(timing_err),
              med_mag_rel_err = median(mag_rel_err),
              med_abs_mag_rel_err = median(abs(mag_rel_err)), .groups = "drop")
  pair_by_overshoot <- pairs %>%
    group_by(method, overshoot) %>%
    summarise(n_matched = n(), med_mag_rel_err = median(mag_rel_err), .groups = "drop")

  ## --- save tables --------------------------------------------------------------------
  write.csv(fl, file.path(CFG$out_dir, "B_flux_per_incubation.csv"), row.names = FALSE)
  write.csv(ev, file.path(CFG$out_dir, "B_events_per_incubation.csv"), row.names = FALSE)
  write.csv(truth_flux, file.path(CFG$out_dir, "B_truth_flux.csv"), row.names = FALSE)
  write.csv(truth_ev, file.path(CFG$out_dir, "B_truth_events.csv"), row.names = FALSE)
  write.csv(flux_summary, file.path(CFG$out_dir, "B_flux_summary.csv"), row.names = FALSE)
  write.csv(flux_by_factor, file.path(CFG$out_dir, "B_flux_by_factor.csv"), row.names = FALSE)
  write.csv(false_ebul, file.path(CFG$out_dir, "B_false_ebullition_no_bubbles.csv"), row.names = FALSE)
  write.csv(det_summary, file.path(CFG$out_dir, "B_detection_summary.csv"), row.names = FALSE)
  write.csv(det_by_snr, file.path(CFG$out_dir, "B_detection_by_snr.csv"), row.names = FALSE)
  write.csv(pair_summary, file.path(CFG$out_dir, "B_timing_magnitude_summary.csv"), row.names = FALSE)
  write.csv(pairs, file.path(CFG$out_dir, "B_matched_events.csv"), row.names = FALSE)

  cat("\n==== B. Flux accuracy (synthetic) ====\n")
  print(as.data.frame(flux_summary %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
  cat("\n==== B. False ebullition in bubble-free incubations ====\n")
  print(as.data.frame(false_ebul %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
  cat("\n==== B. Bubble detection ====\n")
  print(as.data.frame(det_summary %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
  cat("\n==== B. Timing and magnitude of matched bubbles ====\n")
  print(as.data.frame(pair_summary %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
  cat("\n==== B. Per-bubble magnitude error by overshoot ====\n")
  print(as.data.frame(pair_by_overshoot %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
  cat("\n==== B. Run time (s per incubation) ====\n")
  print(as.data.frame(fl %>% group_by(method) %>% summarise(s = round(mean(seconds), 3))), row.names = FALSE)

  ## --- figures -------------------------------------------------------------------------
  p_flux <- flux_long %>%
    filter(!is.na(est), true > 0, est > 0) %>%
    ggplot(aes(true, est)) +
    geom_abline(linetype = 2, colour = "grey40") +
    geom_point(aes(colour = factor(noise_sd)), alpha = 0.5, size = 0.8) +
    scale_x_log10() + scale_y_log10() +
    facet_grid(component ~ method) +
    labs(x = "True flux (nmol m-2 s-1)", y = "Estimated flux", colour = "Noise (ppb)",
         title = "Flux estimates vs truth (non-positive estimates omitted)") +
    theme_bw(base_size = 8) + theme(legend.position = "bottom")

  p_relerr <- flux_long %>%
    filter(!is.na(rel_err)) %>%
    ggplot(aes(method, rel_err, fill = method)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_boxplot(outlier.size = 0.3) +
    coord_cartesian(ylim = c(-1, 2)) +
    facet_grid(component ~ n_bubbles, labeller = label_both) +
    labs(x = NULL, y = "Relative error (est - true) / true") +
    theme_bw(base_size = 8) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                    legend.position = "none")

  p_det <- det_by_snr %>%
    group_by(method, snr) %>%
    summarise(recall = mean(recall), precision = mean(precision, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_longer(c(recall, precision)) %>%
    ggplot(aes(snr, value, colour = method)) +
    geom_line() + geom_point(size = 1) + scale_x_log10() +
    facet_wrap(~ name) +
    labs(x = "Mean bubble step / noise SD", y = NULL, title = "Bubble detection") +
    theme_bw(base_size = 9)

  p_count <- det_per_inc %>%
    group_by(method, n_bubbles) %>%
    summarise(n_est = mean(n_est), .groups = "drop") %>%
    ggplot(aes(factor(n_bubbles), n_est, fill = method)) +
    geom_col(position = "dodge") +
    geom_point(data = distinct(det_per_inc, n_bubbles), inherit.aes = FALSE,
               aes(factor(n_bubbles), n_bubbles), shape = 4, size = 3, colour = "red") +
    labs(x = "True bubbles per incubation (red x)", y = "Mean detected events",
         title = "Number of events") + theme_bw(base_size = 9)

  p_timing <- pairs %>%
    ggplot(aes(method, timing_err, fill = method)) +
    geom_hline(yintercept = 0, linetype = 2) + geom_boxplot(outlier.size = 0.3) +
    labs(x = NULL, y = "Timing error (s)", title = "Bubble timing (matched events)") +
    theme_bw(base_size = 9) + theme(legend.position = "none")

  p_mag <- pairs %>%
    ggplot(aes(method, mag_rel_err, fill = method)) +
    geom_hline(yintercept = 0, linetype = 2) + geom_boxplot(outlier.size = 0.3) +
    coord_cartesian(ylim = c(-1, 2)) +
    facet_wrap(~ paste("overshoot =", overshoot)) +
    labs(x = NULL, y = "Relative error of bubble magnitude",
         title = "Bubble magnitude (matched events)") +
    theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                    legend.position = "none")

  pdf(file.path(CFG$out_dir, "B_synthetic_figures.pdf"), width = 11, height = 8)
  print(p_flux); print(p_relerr); print(p_det); print(p_count); print(p_timing); print(p_mag)
  ## example incubations: one per bubble count, mid noise, overshoot 1
  mid <- function(x) { u <- sort(unique(x[!is.na(x)])); u[ceiling(length(u) / 2)] }
  ex_ids <- truth_flux %>%
    filter(rep == 1, noise_sd == mid(noise_sd),
           is.na(overshoot) | overshoot == max(overshoot, na.rm = TRUE),
           is.na(bubble_mean) | bubble_mean == mid(bubble_mean)) %>%
    distinct(n_bubbles, .keep_all = TRUE) %>% pull(id)
  inc_by_id <- setNames(syn$incs, vapply(syn$incs, `[[`, "", "id"))
  for (i in ex_ids) {
    tf <- truth_flux[truth_flux$id == i, ]
    print(plot_incubation(inc_by_id[[i]], ev, truth_ev[truth_ev$id == i, ],
                          title = sprintf("%s: noise %s ppb, slope %s ppb/s, %s bubbles (mean %s ppb), overshoot %s",
                                          i, tf$noise_sd, tf$slope, tf$n_bubbles, tf$bubble_mean, tf$overshoot)))
  }
  dev.off()
}


## ============================================================================
## PART C. REAL INCUBATIONS
## ============================================================================

load_real_incubations <- function(rc = CFG$real) {
  e <- new.env(); load(rc$data_file, envir = e)
  mydata_all <- e[[rc$data_object]]; myauxfile <- e[[rc$aux_object]]
  myauxfile$obs.length <- myauxfile[[rc$aux_duration_col]]
  ids <- sort(unique(myauxfile$UniqueID))
  if (!is.null(rc$ids)) ids <- intersect(ids, rc$ids)
  if (!is.null(rc$n_max) && rc$n_max < length(ids)) { set.seed(1); ids <- sort(sample(ids, rc$n_max)) }

  incs <- list(); skipped <- list()
  for (i in ids) {
    aux_i <- myauxfile[myauxfile$UniqueID == i, ]
    dat_i <- mydata_all[mydata_all$UniqueID == i, ]
    if (!nrow(dat_i)) { skipped[[i]] <- data.frame(UniqueID = i, reason = "no measurements"); next }
    td <- as.numeric(difftime(min(dat_i$POSIX.time), aux_i$start.time[1], units = "secs"))
    if (td < -rc$time_tolerance_s) {
      skipped[[i]] <- data.frame(UniqueID = i, reason = sprintf("aux start %.0f s before data", -td)); next
    }
    inc <- tryCatch({
      IDed <- quietly(autoID(inputfile = dat_i, auxfile = aux_i, shoulder = rc$shoulder))$value
      make_incubation(i, IDed, rc$prec_default)
    }, error = function(err) { skipped[[i]] <<- data.frame(UniqueID = i, reason = conditionMessage(err)); NULL })
    if (!is.null(inc)) incs[[i]] <- inc
  }
  skipped <- if (length(skipped)) bind_rows(skipped) else
    data.frame(UniqueID = character(0), reason = character(0))
  list(incs = incs, skipped = skipped)
}

if (CFG$run_real) {

  stopifnot(exists("autoID"))
  f_real <- file.path(CFG$out_dir, "C_real_results.rds")
  if (CFG$recompute || !file.exists(f_real)) {
    real <- load_real_incubations()
    cat("\nReal incubations loaded:", length(real$incs), "| skipped:", nrow(real$skipped), "\n")
    t_start <- Sys.time()
    res_real <- run_all(real$incs, "real")
    cat("Run time:", format(round(Sys.time() - t_start, 1)), "\n")
    saveRDS(list(real = real, res = res_real), f_real)
  } else {
    tmp <- readRDS(f_real); real <- tmp$real; res_real <- tmp$res
  }

  flr <- res_real$flux %>% mutate(method = factor(method, levels = METHOD_LEVELS))
  evr <- res_real$events
  ids_r <- names(real$incs)
  ref_m <- CFG$real$reference_method

  ## --- C1. Availability and events ------------------------------------------------
  avail <- flr %>%
    group_by(method) %>%
    summarise(n = n(), errors = sum(status == "error"),
              total_available = mean(!is.na(total)),
              diffusive_available = mean(!is.na(diffusive)),
              ebullition_available = mean(!is.na(ebullition)),
              share_with_events = mean(n_events > 0, na.rm = TRUE),
              median_events = median(n_events, na.rm = TRUE),
              mean_events = mean(n_events, na.rm = TRUE),
              median_ebullition_fraction = median(ifelse(total > 0, ebullition / total, NA), na.rm = TRUE),
              s_per_incubation = mean(seconds), .groups = "drop")

  ## --- C2. Total flux vs endpoint (model-free mass balance) ----------------------------
  ep <- flr %>% filter(method == "endpoint") %>% select(id, endpoint = total)
  vs_endpoint <- flr %>%
    filter(method != "endpoint") %>%
    left_join(ep, by = "id") %>%
    group_by(method) %>%
    summarise(n = sum(is.finite(total) & is.finite(endpoint)),
              median_ratio = median(total / endpoint, na.rm = TRUE),
              iqr_ratio_low = quantile(total / endpoint, 0.25, na.rm = TRUE),
              iqr_ratio_high = quantile(total / endpoint, 0.75, na.rm = TRUE),
              ccc = ccc(total, endpoint),
              share_above_1.2 = mean(total / endpoint > 1.2, na.rm = TRUE), .groups = "drop")

  ## --- C3. Pairwise agreement with the reference method --------------------------
  pw <- flr %>%
    select(id, method, total, diffusive, ebullition) %>%
    pivot_longer(c(total, diffusive, ebullition), names_to = "component", values_to = "est")
  pw_ref <- pw %>% filter(method == ref_m) %>% select(id, component, ref = est)
  pairwise <- pw %>%
    filter(method != ref_m, method != "endpoint") %>%
    inner_join(pw_ref, by = c("id", "component")) %>%
    group_by(method, component) %>%
    filter(any(!is.na(est))) %>%
    summarise(n = sum(is.finite(est) & is.finite(ref)),
              median_ratio = median(est / ref, na.rm = TRUE),
              spearman = suppressWarnings(cor(est, ref, method = "spearman", use = "complete.obs")),
              ccc = ccc(est, ref), .groups = "drop")

  ## --- C4. Event agreement between methods -------------------------------------------
  ev_methods_r <- setdiff(CFG$methods, "goFlux_noSep")
  ref_ev <- evr %>% filter(method == ref_m) %>% select(id, t_event, magnitude_ppb)
  ev_agree <- evaluate_events(evr, ref_ev, ids_r, setdiff(ev_methods_r, ref_m))
  ev_agree_summary <- ev_agree$per_inc %>%
    group_by(method) %>%
    summarise(events_ref = sum(n_ref), events_method = sum(n_est), matched = sum(TP),
              share_ref_matched = if (sum(n_ref) > 0) sum(TP) / sum(n_ref) else NA_real_,
              share_method_matched = if (sum(n_est) > 0) sum(TP) / sum(n_est) else NA_real_,
              .groups = "drop") %>%
    left_join(ev_agree$pairs %>% group_by(method) %>%
                summarise(med_dt_s = median(t_est - t_ref),
                          med_mag_ratio = median(mag_est / mag_ref), .groups = "drop"),
              by = "method")
  ## incubation-level agreement on "has ebullition"
  has_ev <- flr %>% filter(method %in% ev_methods_r) %>%
    transmute(id, method, has = n_events > 0) %>%
    pivot_wider(names_from = method, values_from = has)

  ## --- save + print -------------------------------------------------------------------------
  write.csv(real$skipped, file.path(CFG$out_dir, "C_skipped_incubations.csv"), row.names = FALSE)
  write.csv(flr, file.path(CFG$out_dir, "C_flux_per_incubation.csv"), row.names = FALSE)
  write.csv(evr, file.path(CFG$out_dir, "C_events_per_incubation.csv"), row.names = FALSE)
  write.csv(avail, file.path(CFG$out_dir, "C_availability_events.csv"), row.names = FALSE)
  write.csv(vs_endpoint, file.path(CFG$out_dir, "C_total_vs_endpoint.csv"), row.names = FALSE)
  write.csv(pairwise, file.path(CFG$out_dir, "C_pairwise_vs_reference.csv"), row.names = FALSE)
  write.csv(ev_agree_summary, file.path(CFG$out_dir, "C_event_agreement.csv"), row.names = FALSE)
  write.csv(has_ev, file.path(CFG$out_dir, "C_has_ebullition_by_method.csv"), row.names = FALSE)

  cat("\n==== C. Availability and events (real) ====\n")
  print(as.data.frame(avail %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
  cat("\n==== C. Total flux vs endpoint estimate ====\n")
  print(as.data.frame(vs_endpoint %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
  cat("\n==== C. Agreement with", ref_m, "====\n")
  print(as.data.frame(pairwise %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
  cat("\n==== C. Event agreement with", ref_m, "====\n")
  print(as.data.frame(ev_agree_summary %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)

  ## --- figures ----------------------------------------------------------------------------
  p_ep <- flr %>% filter(method != "endpoint") %>% left_join(ep, by = "id") %>%
    filter(total > 0, endpoint > 0) %>%
    ggplot(aes(endpoint, total)) +
    geom_abline(linetype = 2, colour = "grey40") +
    geom_point(alpha = 0.4, size = 0.8, colour = "#E29338") +
    scale_x_log10() + scale_y_log10() + facet_wrap(~ method) +
    labs(x = "Endpoint total flux (nmol m-2 s-1)", y = "Method total flux",
         title = "Total CH4 flux vs model-free endpoint estimate") + theme_bw(base_size = 9)

  p_pw <- pw %>%
    filter(method != ref_m, method != "endpoint") %>%
    inner_join(pw_ref, by = c("id", "component")) %>%
    filter(est > 0, ref > 0) %>%
    ggplot(aes(ref, est)) +
    geom_abline(linetype = 2, colour = "grey40") +
    geom_point(alpha = 0.4, size = 0.7) +
    scale_x_log10() + scale_y_log10() +
    facet_grid(component ~ method) +
    labs(x = paste(ref_m, "(nmol m-2 s-1)"), y = "Method", title = paste("Agreement with", ref_m)) +
    theme_bw(base_size = 8)

  p_nev <- flr %>% filter(method %in% ev_methods_r) %>%
    ggplot(aes(pmin(n_events, 20), fill = method)) +
    geom_histogram(binwidth = 1, position = "dodge") +
    labs(x = "Events per incubation (capped at 20)", y = "Incubations",
         title = "Number of detected ebullition events") + theme_bw(base_size = 9)

  p_frac <- flr %>% filter(method %in% ev_methods_r, total > 0) %>%
    mutate(frac = pmin(pmax(ebullition / total, 0), 1)) %>%
    ggplot(aes(method, frac, fill = method)) + geom_boxplot(outlier.size = 0.3) +
    labs(x = NULL, y = "Ebullition / total", title = "Ebullition fraction") +
    theme_bw(base_size = 9) + theme(legend.position = "none")

  pdf(file.path(CFG$out_dir, "C_real_figures.pdf"), width = 11, height = 8)
  print(p_ep); print(p_pw); print(p_nev); print(p_frac)
  set.seed(2)
  for (i in sample(ids_r, min(CFG$real$n_example_plots, length(ids_r))))
    print(plot_incubation(real$incs[[i]], evr, title = i))
  dev.off()
}

cat("\nOutputs in", normalizePath(CFG$out_dir), "\n")
