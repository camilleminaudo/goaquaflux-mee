###############################################################################
## compare_bubble_detection.R   (supplementary; development benchmark)
##
## Benchmark of the bubble-detection STEP ALONE, used while developing
## find.bubbles() to compare detection variants (dispersion metric, window
## size, magnitude model) with the detection of FluxSeparator. It is not
## needed to reproduce the figures and tables of the manuscript, which rely on
## compare_goAquaFlux_FluxSeparator_MSP.R (full flux estimates of all methods).
## Its synthetic design differs from that benchmark (ramp = 1 s; noise 2, 20
## and 200 ppb; slopes 0.05, 0.5 and 3 ppb s-1).
##
## Two ebullition-detection approaches are compared on synthetic floating-
## chamber incubations with a known ground truth:
##
##   * goFlux::find.bubbles()  (rolling dispersion + adaptive threshold +
##                              step / step+re-equilibration regression)
##   * FluxSeparator            (running variance of concentration + absolute
##                              cutoff + last-minus-min magnitude)
##
## Each synthetic incubation = baseline + linear diffusive trend + bubbles
## (short ramp, settled step, optional exponentially decaying overshoot)
## + Gaussian (optionally AR(1)) noise. All concentrations are in ppb.
##
## Metrics (per scenario x detector):
##   recall            fraction of true bubbles covered by a detected event
##   fp_per_inc        detected events that cover no true bubble, per incubation
##   false_alarm       share of bubble-FREE incubations with >= 1 detection
##   merge_ratio       n detected events / n true bubbles (< 1 = merging)
##   relbias_total     (sum est. magnitudes - sum true settled steps) / sum true
##                     -> this is the quantity that becomes the ebullitive flux
##   relbias_event     per-event relative bias, one-to-one matches only
##
## Detectors compared (edit `detectors` in section 4):
##   FS_default   FluxSeparator as shipped (cutoffs tuned for ppm sensors:
##                data are passed in ppm, magnitudes converted back to ppb)
##   FS_oracle    FluxSeparator with cutoffs scaled using the TRUE noise and
##                slope (best case; not available in practice)
##   gF_diff      find.bubbles, method "diff", window 30
##   gF_diff_w15  find.bubbles, method "diff", window 15 (default used in the
##                manuscript)
##   gF_variance  find.bubbles, method "variance", window 30
##   gF_step      find.bubbles, "diff", plain step model, no ramp exclusion
##                (ablation of the re-equilibration model)
##
## FluxSeparator is re-implemented at event level (fs_events()) because
## ebullitive_flux() only returns per-cycle sums. Section 7 (optional) checks
## that this re-implementation reproduces ebullitive_flux() exactly.
###############################################################################


## ---- 0. Settings -----------------------------------------------------------

## find.bubbles() is taken from goFlux, loaded by R/setup.R (installed package,
## or a local clone given by the environment variable GOFLUX_DIR).
source(file.path("R", "setup.R"))

## Optional: path to FluxSeparator's R/ folder (v2.0.0, source code from
## https://github.com/JonasStage/FluxSeparator) to run the fidelity check of
## section 7. Set to NULL to skip.
FLUXSEP_R_DIR <- NULL

N_REP    <- 30          # replicate incubations per scenario
N_CORES  <- 1           # > 1 uses parallel::mclapply (not on Windows)
OUT_DIR  <- file.path("results", "supplementary_bubble_detection_benchmark")
BASE_SEED <- 20260924
MATCH_TOL <- 5          # s; tolerance when matching true bubbles to events

suppressPackageStartupMessages({
  library(zoo)      # used by find.bubbles
  library(TTR)      # runVar, as in FluxSeparator
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

stopifnot(exists("find.bubbles"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


## ---- 1. Scenario grid ------------------------------------------------------

## Noise spans laser analysers (~2 ppb) to low-cost / DIY sensors (~200 ppb).
## Slopes span low to very high diffusive CH4 accumulation.
## `ramp` (s) = time for a bubble to reach the transient peak at the analyser.
## Results for the "diff" metric are sensitive to it (a step spread over r
## increments gives ~r times less squared increment): consider re-running with
## ramp = 1 and ramp = 6 as a sensitivity analysis.
fixed <- list(duration = 600, dt = 1, C0 = 2000, tau = 15, ramp = 1,
              edge = 60, min_spacing = 20, size_sdlog = 0.5, ar1 = 0)

grid_bub <- expand.grid(noise_sd    = c(2, 20, 200),
                        slope       = c(0.05, 0.5, 3),     # ppb s-1
                        bubble_mean = c(50, 250, 1500),    # ppb (settled step)
                        n_bubbles   = c(1, 4, 10),
                        overshoot   = c(0, 1))             # peak excess / step
grid_null <- expand.grid(noise_sd = c(2, 20, 200),
                         slope    = c(0.05, 0.5, 3),
                         bubble_mean = NA, n_bubbles = 0, overshoot = NA)
scenarios <- rbind(grid_bub, grid_null)
scenarios$scenario <- seq_len(nrow(scenarios))
cat(nrow(scenarios), "scenarios x", N_REP, "replicates =",
    nrow(scenarios) * N_REP, "incubations\n")


## ---- 2. Simulator ----------------------------------------------------------

## Shape of one bubble as seen by the analyser, u = time since bubble onset.
## Linear ramp to the transient peak (M + A), then exponential decay of the
## overshoot A towards the settled step M.
bubble_shape <- function(u, M, A, tau, ramp) {
  out <- numeric(length(u))
  if (ramp > 0) {
    r <- u >= 0 & u < ramp
    out[r] <- (M + A) * u[r] / ramp
  }
  p <- u >= ramp
  out[p] <- M + A * exp(-(u[p] - ramp) / tau)
  out
}

## Random bubble onsets in [lower, upper] with a minimum spacing.
place_bubbles <- function(n, lower, upper, min_spacing, max_try = 5000) {
  if (n == 0) return(numeric(0))
  for (k in seq_len(max_try)) {
    tb <- sort(runif(n, lower, upper))
    if (n == 1 || all(diff(tb) >= min_spacing)) return(round(tb))
  }
  stop("Could not place ", n, " bubbles with spacing ", min_spacing)
}

simulate_incubation <- function(sc, fx = fixed) {
  t  <- seq(0, fx$duration - fx$dt, by = fx$dt)
  n  <- sc$n_bubbles
  tb <- place_bubbles(n, fx$edge, fx$duration - fx$edge, fx$min_spacing)
  ## Log-normal sizes with arithmetic mean = bubble_mean
  M  <- if (n > 0) sc$bubble_mean *
    exp(rnorm(n, -fx$size_sdlog^2 / 2, fx$size_sdlog)) else numeric(0)

  signal <- fx$C0 + sc$slope * t
  for (i in seq_len(n)) {
    signal <- signal + bubble_shape(t - tb[i], M[i], sc$overshoot * M[i],
                                    fx$tau, fx$ramp)
  }
  eps <- if (fx$ar1 > 0) {
    as.numeric(arima.sim(list(ar = fx$ar1), n = length(t),
                         sd = sc$noise_sd * sqrt(1 - fx$ar1^2)))
  } else rnorm(length(t), 0, sc$noise_sd)

  list(data  = data.frame(Etime = t, CH4dry_ppb = signal + eps),
       truth = data.frame(t_bubble = tb, magnitude = M))
}


## ---- 3. Detector wrappers --------------------------------------------------

empty_events <- function() data.frame(start = numeric(0), end = numeric(0),
                                      magnitude = numeric(0))

## Event-level re-implementation of FluxSeparator v2.0.0
## (detect_bubbles() + the classification / summary steps of
## ebullitive_flux()), for ONE pump cycle sampled at ~1 Hz. Quirks of the
## original are kept on purpose (marked "as in FS").
fs_events <- function(time, conc, runvar_cutoff = 0.5, IndexSpan = 30,
                      concentration_diffusion_cutoff = 1,
                      top_selection = c("last", "max"),
                      runvar_window = 5, time_gap_seconds = 30,
                      min_obs = 100) {
  top_selection <- match.arg(top_selection)
  ok <- !is.na(conc); time <- time[ok]; conc <- conc[ok]
  if (length(conc) <= min_obs) return(empty_events())

  rv <- TTR::runVar(conc, n = runvar_window)          # right-aligned, n-1
  flagged <- which(rv > runvar_cutoff)
  flagged <- flagged[-1]   # as in FS: drop_na() on the lagged time diff
  if (length(flagged) == 0) return(empty_events())

  ids <- unique(unlist(lapply(flagged, function(x) {
    v <- (x - IndexSpan):(x + IndexSpan); v[v > 0] })))
  ids <- sort(ids[ids <= length(conc)])

  tt <- time[ids]; cc <- conc[ids]
  td <- c(NA, diff(tt))
  keep <- !is.na(td)       # as in FS: first padded row dropped again
  tt <- tt[keep]; cc <- cc[keep]; td <- td[keep]
  if (length(tt) == 0) return(empty_events())
  grp <- 1 + cumsum(td > time_gap_seconds)

  ev <- lapply(split(seq_along(tt), grp), function(i) {
    x <- cc[i]; ti <- tt[i]
    if (!(x[1] < x[length(x)])) return(NULL)         # net increase only
    top <- if (top_selection == "last") x[length(x)] else max(x)
    d   <- top - min(x)
    if (!(d > concentration_diffusion_cutoff &&
          ti[which.min(x)] < ti[which.max(x)])) return(NULL)
    data.frame(start = min(ti), end = max(ti), magnitude = d)
  })
  ev <- do.call(rbind, ev)
  if (is.null(ev)) empty_events() else ev
}

gf_events <- function(time, conc, ...) {
  df  <- data.frame(Etime = time, CH4dry_ppb = conc)
  res <- tryCatch(suppressWarnings(
    find.bubbles(df, bubble_source = "CH4dry_ppb", ...)),
    error = function(e) { message("find.bubbles error: ", conditionMessage(e)); NULL })
  if (is.null(res) || nrow(res) == 0) return(empty_events())
  data.frame(start = res$start, end = res$end, magnitude = res$magnitude)
}


## ---- 4. Detectors to compare -----------------------------------------------
## Each takes the simulated data frame and the scenario row, returns events
## (start, end, magnitude in ppb).

detectors <- list(
  FS_default = function(d, sc) {
    ev <- fs_events(d$Etime, d$CH4dry_ppb / 1000)     # FS defaults assume ppm
    ev$magnitude <- ev$magnitude * 1000
    ev
  },
  FS_oracle = function(d, sc) {
    s_step <- sc$slope * fixed$dt
    fs_events(d$Etime, d$CH4dry_ppb,
              ## E[runVar] of a ramp + noise = 2.5 s^2 + sigma^2 (n = 5);
              ## 6 sigma^2 keeps pure-noise exceedances rare over 600 samples
              runvar_cutoff = 2.5 * s_step^2 + 6 * sc$noise_sd^2,
              ## reject diffusion accumulated over the padded window + noise
              concentration_diffusion_cutoff =
                s_step * 2 * 30 + 5 * sc$noise_sd)
  },
  gF_diff     = function(d, sc) gf_events(d$Etime, d$CH4dry_ppb,
                                          window.size = 30, method = "diff"),
  gF_diff_w15 = function(d, sc) gf_events(d$Etime, d$CH4dry_ppb,
                                          window.size = 15, method = "diff"),
  gF_variance = function(d, sc) gf_events(d$Etime, d$CH4dry_ppb,
                                          window.size = 30, method = "variance"),
  gF_step     = function(d, sc) gf_events(d$Etime, d$CH4dry_ppb,
                                          window.size = 30, method = "diff",
                                          magnitude.model = "step",
                                          exclude.ramp = FALSE)
)


## ---- 5. Scoring ------------------------------------------------------------

score_events <- function(ev, truth, tol = MATCH_TOL) {
  n_true <- nrow(truth); n_det <- nrow(ev)
  cover <- if (n_true > 0 && n_det > 0) {
    outer(truth$t_bubble, ev$start - tol, ">=") &
      outer(truth$t_bubble, ev$end + tol, "<=")      # [true x event]
  } else matrix(FALSE, n_true, n_det)

  inc <- data.frame(
    n_true    = n_true,
    n_det     = n_det,
    n_matched = sum(rowSums(cover) > 0),
    n_fp      = sum(colSums(cover) == 0),
    sum_true  = sum(truth$magnitude),
    sum_est   = sum(ev$magnitude))                    # incl. false positives

  ## One-to-one matches for per-event magnitude bias
  ev_rows <- NULL
  if (n_true > 0 && n_det > 0) {
    one <- which(cover, arr.ind = TRUE)
    one <- one[rowSums(cover)[one[, 1]] == 1 & colSums(cover)[one[, 2]] == 1, ,
               drop = FALSE]
    if (nrow(one) > 0) {
      ev_rows <- data.frame(true_mag = truth$magnitude[one[, 1]],
                            est_mag  = ev$magnitude[one[, 2]])
    }
  }
  list(inc = inc, ev = ev_rows)
}


## ---- 6. Run ----------------------------------------------------------------

jobs <- expand.grid(scenario = scenarios$scenario, rep = seq_len(N_REP))

run_job <- function(j) {
  sc <- scenarios[scenarios$scenario == jobs$scenario[j], ]
  rownames(sc) <- NULL
  set.seed(BASE_SEED + j)                 # same data for every detector
  sim <- simulate_incubation(sc)
  out_inc <- list(); out_ev <- list()
  for (dn in names(detectors)) {
    t0 <- proc.time()[["elapsed"]]
    ev <- detectors[[dn]](sim$data, sc)
    el <- proc.time()[["elapsed"]] - t0
    s  <- score_events(ev, sim$truth)
    out_inc[[dn]] <- cbind(job = j, rep = jobs$rep[j], sc, detector = dn,
                           s$inc, seconds = el)
    if (!is.null(s$ev)) out_ev[[dn]] <- cbind(job = j, sc, detector = dn, s$ev)
  }
  list(inc = do.call(rbind, out_inc), ev = do.call(rbind, out_ev))
}

t_start <- Sys.time()
res <- if (N_CORES > 1 && .Platform$OS.type != "windows") {
  parallel::mclapply(seq_len(nrow(jobs)), run_job, mc.cores = N_CORES)
} else {
  pb <- txtProgressBar(0, nrow(jobs), style = 3)
  r <- lapply(seq_len(nrow(jobs)), function(j) { setTxtProgressBar(pb, j); run_job(j) })
  close(pb); r
}
cat("\nRun time:", format(round(Sys.time() - t_start, 1)), "\n")

inc <- bind_rows(lapply(res, `[[`, "inc"))
evs <- bind_rows(lapply(res, `[[`, "ev"))
inc$detector <- factor(inc$detector, levels = names(detectors))
if (nrow(evs)) evs$detector <- factor(evs$detector, levels = names(detectors))

write.csv(inc, file.path(OUT_DIR, "per_incubation.csv"), row.names = FALSE)
write.csv(evs, file.path(OUT_DIR, "per_event_matches.csv"), row.names = FALSE)


## ---- 6b. Summaries ---------------------------------------------------------

summ_bub <- inc %>%
  filter(n_bubbles > 0) %>%
  group_by(detector, noise_sd, slope, bubble_mean, n_bubbles, overshoot) %>%
  summarise(recall        = sum(n_matched) / sum(n_true),
            fp_per_inc    = mean(n_fp),
            merge_ratio   = sum(n_det) / sum(n_true),
            relbias_total = median((sum_est - sum_true) / sum_true),
            relbias_q25   = quantile((sum_est - sum_true) / sum_true, 0.25),
            relbias_q75   = quantile((sum_est - sum_true) / sum_true, 0.75),
            .groups = "drop")

summ_null <- inc %>%
  filter(n_bubbles == 0) %>%
  group_by(detector, noise_sd, slope) %>%
  summarise(false_alarm      = mean(n_det > 0),
            fp_per_inc       = mean(n_det),
            false_dC_ppb_med = median(sum_est),
            .groups = "drop")

summ_event <- if (nrow(evs)) evs %>%
  group_by(detector, noise_sd, slope, bubble_mean, overshoot) %>%
  summarise(n_matches     = n(),
            relbias_event = median((est_mag - true_mag) / true_mag),
            .groups = "drop") else NULL

write.csv(summ_bub,  file.path(OUT_DIR, "summary_bubbles.csv"),  row.names = FALSE)
write.csv(summ_null, file.path(OUT_DIR, "summary_null.csv"),     row.names = FALSE)
if (!is.null(summ_event))
  write.csv(summ_event, file.path(OUT_DIR, "summary_event_bias.csv"), row.names = FALSE)

## Compact console overview, marginalised over the other factors.
## Printed twice: all bubble scenarios, and only those where the mean bubble
## step is >= 10 x noise SD (i.e. where detection should be feasible).
make_overview <- function(x) {
  x %>%
    group_by(detector) %>%
    summarise(recall          = round(sum(n_matched) / sum(n_true), 3),
              fp_per_inc      = round(mean(n_fp), 2),
              merge_ratio     = round(sum(n_det) / sum(n_true), 2),
              med_relbias     = round(median((sum_est - sum_true) / sum_true), 3),
              med_abs_relbias = round(median(abs(sum_est - sum_true) / sum_true), 3),
              ms_per_inc      = round(1000 * mean(seconds), 1), .groups = "drop") %>%
    left_join(summ_null %>% group_by(detector) %>%
                summarise(false_alarm_null = round(mean(false_alarm), 3)),
              by = "detector")
}
cat("\n==== Overview: all bubble scenarios pooled ====\n")
print(as.data.frame(make_overview(filter(inc, n_bubbles > 0))), row.names = FALSE)
cat("\n==== Overview: scenarios with bubble_mean / noise_sd >= 10 ====\n")
print(as.data.frame(make_overview(filter(inc, n_bubbles > 0,
                                         bubble_mean / noise_sd >= 10))),
      row.names = FALSE)

cat("\n==== Recall by noise x slope ====\n")
print(as.data.frame(inc %>% filter(n_bubbles > 0) %>%
  group_by(detector, noise_sd, slope) %>%
  summarise(recall = round(sum(n_matched) / sum(n_true), 2), .groups = "drop") %>%
  pivot_wider(names_from = detector, values_from = recall)), row.names = FALSE)

cat("\n==== Median relative bias of total bubble dC (step/noise >= 10), by overshoot x n_bubbles ====\n")
print(as.data.frame(inc %>% filter(n_bubbles > 0, bubble_mean / noise_sd >= 10) %>%
  group_by(detector, overshoot, n_bubbles) %>%
  summarise(relbias = round(median((sum_est - sum_true) / sum_true), 3),
            .groups = "drop") %>%
  pivot_wider(names_from = detector, values_from = relbias)), row.names = FALSE)


## ---- 6c. Figures -----------------------------------------------------------

lab_noise <- function(x) paste0("noise ", x, " ppb")
lab_slope <- function(x) paste0("slope ", x, " ppb/s")

p_recall <- summ_bub %>%
  group_by(detector, noise_sd, slope, bubble_mean) %>%
  summarise(recall = mean(recall), .groups = "drop") %>%
  mutate(snr = bubble_mean / noise_sd) %>%
  ggplot(aes(snr, recall, colour = detector)) +
  geom_line() + geom_point() +
  scale_x_log10() +
  facet_grid(slope ~ ., labeller = labeller(slope = lab_slope)) +
  labs(x = "Mean bubble step / noise SD (log)", y = "Recall",
       title = "Detection rate") + theme_bw()

p_null <- summ_null %>%
  ggplot(aes(factor(noise_sd), fp_per_inc, fill = detector)) +
  geom_col(position = "dodge") +
  facet_wrap(~ slope, labeller = labeller(slope = lab_slope)) +
  labs(x = "Noise SD (ppb)", y = "Detections per bubble-free incubation",
       title = "False positives without ebullition") + theme_bw()

p_bias <- inc %>%
  filter(n_bubbles > 0) %>%
  mutate(relbias = (sum_est - sum_true) / sum_true,
         overshoot = paste0("overshoot = ", overshoot)) %>%
  ggplot(aes(detector, relbias, fill = detector)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.size = 0.4) +
  coord_cartesian(ylim = c(-1, 1.5)) +
  facet_grid(overshoot ~ slope, labeller = labeller(slope = lab_slope)) +
  labs(y = "Relative bias of total ebullitive dC", x = NULL,
       title = "Magnitude bias (all incubations with bubbles)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                     legend.position = "none")

p_merge <- summ_bub %>%
  group_by(detector, n_bubbles, noise_sd) %>%
  summarise(merge_ratio = mean(merge_ratio), .groups = "drop") %>%
  ggplot(aes(factor(n_bubbles), merge_ratio, fill = detector)) +
  geom_col(position = "dodge") + geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~ noise_sd, labeller = labeller(noise_sd = lab_noise)) +
  labs(x = "True bubbles per incubation", y = "Detected events / true bubbles",
       title = "Event counting (merging < 1 < splitting)") + theme_bw()

## Example traces with truth and detections (one panel per detector)
plot_example <- function(sc_row, seed = 1) {
  sc <- scenarios[sc_row, ]
  set.seed(seed)
  sim <- simulate_incubation(sc)
  lv <- names(detectors)
  evs_ex <- bind_rows(lapply(lv, function(dn) {
    e <- detectors[[dn]](sim$data, sc)
    if (nrow(e)) cbind(detector = dn, e) else NULL
  }))
  trace <- tidyr::crossing(sim$data, detector = lv)
  trace$detector <- factor(trace$detector, levels = lv)
  p <- ggplot(trace, aes(Etime, CH4dry_ppb))
  if (nrow(evs_ex)) {
    evs_ex$detector <- factor(evs_ex$detector, levels = lv)
    p <- p + geom_rect(data = evs_ex, inherit.aes = FALSE,
                       aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                       fill = "steelblue", alpha = 0.25)
  }
  p + geom_line(linewidth = 0.3) +
    geom_vline(xintercept = sim$truth$t_bubble, colour = "red", linetype = 2) +
    facet_wrap(~ detector, ncol = 1) +
    labs(title = sprintf("noise %s ppb, slope %s ppb/s, step %s ppb, n = %s, overshoot %s",
                         sc$noise_sd, sc$slope, sc$bubble_mean, sc$n_bubbles, sc$overshoot),
         subtitle = "red dashed = true onsets; blue = detected events",
         x = "Etime (s)", y = "CH4 (ppb)") + theme_bw()
}
ex_rows <- with(scenarios, c(
  which(noise_sd == 20 & slope == 0.5 & bubble_mean == 250 & n_bubbles == 4  & overshoot == 1),
  which(noise_sd == 2  & slope == 3   & bubble_mean == 50  & n_bubbles == 4  & overshoot == 0),
  which(noise_sd == 20 & slope == 0.5 & bubble_mean == 250 & n_bubbles == 10 & overshoot == 1)))

pdf(file.path(OUT_DIR, "benchmark_figures.pdf"), width = 10, height = 7)
print(p_recall); print(p_null); print(p_bias); print(p_merge)
for (r in ex_rows) print(plot_example(r))
dev.off()
cat("\nOutputs written to", normalizePath(OUT_DIR), "\n")


## ---- 7. Optional: fidelity check of fs_events() vs FluxSeparator -----------
## Runs FluxSeparator::ebullitive_flux() (sourced from FLUXSEP_R_DIR, so the
## package itself need not be installed) on a sample of incubations and
## compares its per-cycle sum_bubbles_concentration with sum(fs_events()).

if (!is.null(FLUXSEP_R_DIR)) {
  suppressPackageStartupMessages({ library(purrr); library(rlang); library(magrittr) })
  fs_env <- new.env()
  for (f in c("utils-bubble-detection.R", "utils-smoothing.R", "ebullitive_flux.R"))
    sys.source(file.path(FLUXSEP_R_DIR, f), envir = fs_env)

  chk <- lapply(sample(nrow(jobs), min(60, nrow(jobs))), function(j) {
    sc <- scenarios[scenarios$scenario == jobs$scenario[j], ]
    set.seed(BASE_SEED + j)
    d <- simulate_incubation(sc)$data
    ## ppb data with the FS_oracle cutoffs, so that detections do occur
    rvc <- 2.5 * sc$slope^2 + 6 * sc$noise_sd^2
    cdc <- sc$slope * 60 + 5 * sc$noise_sd
    fsdf <- data.frame(datetime = as.POSIXct("2026-01-01", tz = "UTC") + d$Etime,
                       PumpCycle = 1, station = "sim", tempC = 20, sensor = 1,
                       pred_CH4 = d$CH4dry_ppb)
    ref <- suppressWarnings(suppressMessages(
      fs_env$ebullitive_flux(fsdf, show_plots = FALSE, runvar_cutoff = rvc,
                             concentration_diffusion_cutoff = cdc)))
    mine <- sum(fs_events(d$Etime, d$CH4dry_ppb, runvar_cutoff = rvc,
                          concentration_diffusion_cutoff = cdc)$magnitude)
    data.frame(job = j, FluxSeparator = ref$sum_bubbles_concentration[1],
               fs_events = mine)
  })
  chk <- bind_rows(chk)
  cat("\n==== Fidelity check (ppb, FS_oracle cutoffs) ====\n")
  cat("Incubations checked:", nrow(chk),
      "| max abs difference:", signif(max(abs(chk$FluxSeparator - chk$fs_events)), 3),
      "| incubations with detections:", sum(chk$FluxSeparator > 0), "\n")
}
