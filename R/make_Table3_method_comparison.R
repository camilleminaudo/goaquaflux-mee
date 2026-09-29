###############################################################################
## make_Table3_method_comparison.R
##
## Computes Table 3 and every number quoted in Section 5 ("Comparison with
## other existing methods") from the raw results of
## compare_goAquaFlux_FluxSeparator_MSP.R. No method is re-run here.
##
## Synthetic incubations (known truth):
##   * flux errors: relative error e = (estimate - truth) / truth, summarised
##     by its median (bias) and the median of |e|, for the total, diffusive
##     and ebullitive fluxes; share of incubations without a total flux;
##   * conventional flux reported as diffusive: median e with 1 and 10 bubbles;
##   * bubble identification: true and detected bubbles are matched one-to-one
##     when their times differ by <= 15 s (greedy, closest pairs first);
##     recall = matched / true, precision = matched / detected, events per
##     bubble-free incubation, timing error and relative magnitude error of
##     matched bubbles;
##   * detection limit: share of isolated bubbles (1 per incubation) detected,
##     by bubble size relative to the noise SD;
##   * mass-balance check of goAquaFlux (ebullition.check): total-flux error
##     of incubations failing / passing the closure, among incubations with
##     detected bubbles;
##   * sensitivity to the detection window (goAquaFlux_w30).
## Field incubations (no truth):
##   * share of incubations with ebullition and without a flux estimate;
##   * ebullitive share of the CH4 flux; conventional / diffusive flux ratio;
##   * share of FluxSeparator and MSP events coinciding with a goAquaFlux event;
##   * outcome of the ebullition check.
##
## Inputs   results/method_comparison_synthetic/B_synthetic_results.rds
##          results/method_comparison_real_incubations/C_real_results.rds
## Outputs  results/Table3_method_comparison.csv
##          results/Section5_numbers.csv  (quantity, method, value)
##
## Run from the repository root: Rscript R/make_Table3_method_comparison.R
###############################################################################

if (!dir.exists("results")) stop("Please run this script from the repository root.")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

f_syn  <- file.path("results", "method_comparison_synthetic", "B_synthetic_results.rds")
f_real <- file.path("results", "method_comparison_real_incubations", "C_real_results.rds")
match_tol_s <- 15         # as CFG$match_tol_s in compare_goAquaFlux_FluxSeparator_MSP.R
table_methods <- c("goAquaFlux", "goFlux_noSep", "FS_default", "MSP")

`%||%` <- function(a, b) if (is.null(a)) b else a
numbers <- list()        # collected Section 5 numbers
add <- function(quantity, method, value) {
  numbers[[length(numbers) + 1]] <<- data.frame(quantity = quantity, method = method,
                                                value = unname(value))
}
## Formatting for Table 3 (ASCII only; "n/a" = not applicable)
pct_signed <- function(x, digits = 0) {             # signed error, e.g. "+8%", "-11%"
  if (is.na(x)) return("n/a")
  v <- round(100 * x, digits); if (v == 0) v <- 0          # no "-0%"
  paste0(if (v > 0) "+" else "", formatC(v, format = "f", digits = digits), "%")
}
pct_plain <- function(x, digits = 0) {              # share or absolute error, e.g. "14%"
  if (is.na(x)) return("n/a")
  paste0(formatC(round(100 * x, digits), format = "f", digits = digits), "%")
}

## Greedy one-to-one matching of event times (same rule as the benchmark)
match_events <- function(t_ref, t_est, tol = match_tol_s) {
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

## Match the events of `methods` to reference events, incubation by incubation
match_all <- function(events, ref, ids, methods) {
  ev_split  <- split(events, list(events$id, events$method), drop = TRUE)
  ref_split <- split(ref, ref$id)
  empty <- data.frame(t_event = numeric(0), magnitude_ppb = numeric(0))
  per <- pairs <- per_ref <- list()
  for (i in ids) for (m in methods) {
    r <- ref_split[[i]] %||% empty
    e <- ev_split[[paste(i, m, sep = ".")]] %||% empty
    mt <- match_events(r$t_event, e$t_event)
    per[[length(per) + 1]] <- data.frame(id = i, method = m, n_ref = nrow(r),
                                         n_est = nrow(e), TP = nrow(mt))
    if (nrow(mt)) pairs[[length(pairs) + 1]] <- data.frame(
      id = i, method = m, dt = e$t_event[mt[, 2]] - r$t_event[mt[, 1]],
      mag_ref = r$magnitude_ppb[mt[, 1]], mag_est = e$magnitude_ppb[mt[, 2]])
    if (nrow(r)) per_ref[[length(per_ref) + 1]] <- data.frame(
      id = i, method = m, magnitude_ppb = r$magnitude_ppb,
      detected = seq_len(nrow(r)) %in% mt[, 1])
  }
  list(per = bind_rows(per), pairs = bind_rows(pairs), per_ref = bind_rows(per_ref))
}


## ============================================================================
## SYNTHETIC INCUBATIONS
## ============================================================================
B  <- readRDS(f_syn)
fl <- B$res$flux
ev <- B$res$events
tf <- B$syn$truth_flux
te <- B$syn$truth_ev
syn_methods <- intersect(c(table_methods, "goAquaFlux_w30", "endpoint"), unique(fl$method))

## ---- Flux errors ---------------------------------------------------------------
L <- fl %>%
  filter(method %in% syn_methods) %>%
  select(id, method, total, diffusive, ebullition) %>%
  pivot_longer(c(total, diffusive, ebullition), names_to = "component", values_to = "est") %>%
  left_join(tf, by = "id") %>%
  mutate(true = case_when(component == "total"     ~ true_total,
                          component == "diffusive" ~ true_diffusive,
                          TRUE                     ~ true_ebullition),
         e = ifelse(true > 0, (est - true) / true, NA_real_)) %>%
  group_by(method, component) %>% filter(any(!is.na(est))) %>% ungroup()

flux_err <- L %>%
  group_by(method, component) %>%
  summarise(median_e = median(e, na.rm = TRUE),
            median_abs_e = median(abs(e), na.rm = TRUE), .groups = "drop")
for (k in seq_len(nrow(flux_err))) {
  add(paste0("synthetic: median relative error, ", flux_err$component[k], " flux"),
      flux_err$method[k], flux_err$median_e[k])
  add(paste0("synthetic: median absolute relative error, ", flux_err$component[k], " flux"),
      flux_err$method[k], flux_err$median_abs_e[k])
}

no_estimate_syn <- fl %>% filter(method %in% syn_methods) %>%
  group_by(method) %>% summarise(v = mean(is.na(total)), .groups = "drop")
for (k in seq_len(nrow(no_estimate_syn)))
  add("synthetic: share without a total flux", no_estimate_syn$method[k], no_estimate_syn$v[k])

## Diffusive flux of FluxSeparator and MSP missing
for (m in intersect(c("FS_default", "MSP"), syn_methods))
  add("synthetic: share without a diffusive flux", m, mean(is.na(fl$diffusive[fl$method == m])))

## Conventional flux (no separation) reported as a diffusive flux
nosep_as_diff <- fl %>% filter(method == "goFlux_noSep") %>% left_join(tf, by = "id") %>%
  filter(n_bubbles > 0) %>%
  mutate(e = (total - true_diffusive) / true_diffusive) %>%
  group_by(n_bubbles) %>% summarise(median_e = median(e, na.rm = TRUE), .groups = "drop")
for (k in seq_len(nrow(nosep_as_diff)))
  add(paste0("synthetic: conventional flux as diffusive, median relative error, ",
             nosep_as_diff$n_bubbles[k], " bubble(s)"), "goFlux_noSep", nosep_as_diff$median_e[k])

## goAquaFlux ebullitive error by number of bubbles
ebul_by_n <- L %>% filter(method == "goAquaFlux", component == "ebullition", n_bubbles > 0) %>%
  group_by(n_bubbles) %>% summarise(v = median(e, na.rm = TRUE), .groups = "drop")
for (k in seq_len(nrow(ebul_by_n)))
  add(paste0("synthetic: median relative error, ebullitive flux, ", ebul_by_n$n_bubbles[k],
             " bubble(s)"), "goAquaFlux", ebul_by_n$v[k])

## ---- Bubble identification -----------------------------------------------------------
det_methods <- intersect(c("goAquaFlux", "goAquaFlux_w30", "FS_default", "MSP"), syn_methods)
M <- match_all(ev, te, tf$id, det_methods)
per <- M$per %>% left_join(tf, by = "id")

det <- per %>% group_by(method) %>%
  summarise(recall = sum(TP) / sum(n_ref),
            precision = if (sum(n_est) > 0) sum(TP) / sum(n_est) else NA_real_,
            events_bubble_free = mean(n_est[n_bubbles == 0]), .groups = "drop")
for (k in seq_len(nrow(det))) {
  add("synthetic: bubble recall", det$method[k], det$recall[k])
  add("synthetic: bubble precision", det$method[k], det$precision[k])
  add("synthetic: events per bubble-free incubation", det$method[k], det$events_bubble_free[k])
}

recall_by_noise <- per %>% filter(n_bubbles > 0) %>% group_by(method, noise_sd) %>%
  summarise(v = sum(TP) / sum(n_ref), .groups = "drop")
for (k in seq_len(nrow(recall_by_noise)))
  add(paste0("synthetic: recall at noise SD = ", recall_by_noise$noise_sd[k], " ppb"),
      recall_by_noise$method[k], recall_by_noise$v[k])

recall_by_snr <- per %>% filter(n_bubbles > 0) %>%
  mutate(snr = bubble_mean / noise_sd) %>% group_by(method, snr) %>%
  summarise(v = sum(TP) / sum(n_ref), .groups = "drop")
for (k in seq_len(nrow(recall_by_snr)))
  add(paste0("synthetic: recall at mean bubble / noise = ", recall_by_snr$snr[k]),
      recall_by_snr$method[k], recall_by_snr$v[k])

## Timing and magnitude of matched bubbles
mag <- M$pairs %>% group_by(method) %>%
  summarise(median_dt = median(dt),
            median_mag_e = median((mag_est - mag_ref) / mag_ref),
            median_abs_mag_e = median(abs(mag_est - mag_ref) / mag_ref), .groups = "drop")
for (k in seq_len(nrow(mag))) {
  add("synthetic: median timing error of matched bubbles (s)", mag$method[k], mag$median_dt[k])
  add("synthetic: median relative error of bubble magnitude", mag$method[k], mag$median_mag_e[k])
  add("synthetic: median absolute relative error of bubble magnitude", mag$method[k],
      mag$median_abs_mag_e[k])
}

## Detection limit: isolated bubbles, by size relative to the noise
iso <- M$per_ref %>% left_join(tf, by = "id") %>%
  filter(n_bubbles == 1) %>%
  mutate(size_to_noise = cut(magnitude_ppb / noise_sd, c(0, 5, 10, 15, 20, 30, Inf)))
iso_det <- iso %>% group_by(method, size_to_noise) %>%
  summarise(v = mean(detected), n = n(), .groups = "drop")
for (k in seq_len(nrow(iso_det)))
  add(paste0("synthetic: share of isolated bubbles detected, size / noise in ",
             iso_det$size_to_noise[k]), iso_det$method[k], iso_det$v[k])

## ---- Mass-balance check of goAquaFlux ---------------------------------------------------
if ("ebullition_check" %in% names(fl) && any(!is.na(fl$ebullition_check[fl$method == "goAquaFlux"]))) {
  chk <- fl %>% filter(method == "goAquaFlux", n_events > 0, !is.na(total)) %>%
    left_join(tf, by = "id") %>%
    mutate(e_total = (total - true_total) / true_total,
           e_ebul  = ifelse(true_ebullition > 0, (ebullition - true_ebullition) / true_ebullition, NA),
           failed  = startsWith(ebullition_check, "closure"),
           passed  = ebullition_check == "")
  add("synthetic: incubations with detected bubbles failing the closure", "goAquaFlux",
      sum(chk$failed, na.rm = TRUE))
  add("synthetic: share of failed closures with total-flux error > 20%", "goAquaFlux",
      mean(abs(chk$e_total[chk$failed %in% TRUE]) > 0.2))
  add("synthetic: passed closure, median absolute relative error, total flux", "goAquaFlux",
      median(abs(chk$e_total[chk$passed %in% TRUE])))
  add("synthetic: passed closure, median absolute relative error, ebullitive flux", "goAquaFlux",
      median(abs(chk$e_ebul[chk$passed %in% TRUE]), na.rm = TRUE))
} else {
  message("No ebullition check in the synthetic results (older goFlux version): ",
          "the calibration of the mass-balance check is skipped.")
}


## ============================================================================
## FIELD INCUBATIONS
## ============================================================================
C   <- readRDS(f_real)
flr <- C$res$flux
evr <- C$res$events
ids <- names(C$real$incs)
add("field: incubations analysed", "all", length(ids))
add("field: incubations skipped", "all", nrow(C$real$skipped))

W <- flr %>% select(id, method, total, diffusive, ebullition, n_events) %>%
  pivot_wider(names_from = method, values_from = c(total, diffusive, ebullition, n_events))

for (m in intersect(c("goAquaFlux", "FS_default", "MSP"), unique(flr$method)))
  add("field: share of incubations with ebullition", m,
      mean(W[[paste0("n_events_", m)]] > 0, na.rm = TRUE))
for (m in intersect(table_methods, unique(flr$method)))
  add("field: share without a total flux", m, mean(is.na(W[[paste0("total_", m)]])))
if ("FS_default" %in% flr$method)
  add("field: share without a diffusive flux", "FS_default", mean(is.na(W$diffusive_FS_default)))

## goAquaFlux: ebullitive share and conventional / diffusive ratio
bub <- W$n_events_goAquaFlux > 0 & !is.na(W$n_events_goAquaFlux)
add("field: incubations with ebullition", "goAquaFlux", sum(bub))
fr <- W$ebullition_goAquaFlux / W$total_goAquaFlux
add("field: median ebullitive share of the CH4 flux (bubbling incubations)", "goAquaFlux",
    median(fr[bub & W$total_goAquaFlux > 0], na.rm = TRUE))
add("field: ebullitive share of the CH4 flux summed over the dataset", "goAquaFlux",
    sum(W$ebullition_goAquaFlux, na.rm = TRUE) / sum(W$total_goAquaFlux, na.rm = TRUE))
add("field: incubations with ebullition > 10% of the CH4 flux", "goAquaFlux",
    sum(fr > 0.1 & fr <= 1, na.rm = TRUE))
ok <- bub & !is.na(W$diffusive_goAquaFlux) & W$diffusive_goAquaFlux > 0 & W$total_goFlux_noSep > 0
r <- W$total_goFlux_noSep[ok] / W$diffusive_goAquaFlux[ok]
add("field: median ratio conventional / goAquaFlux diffusive flux (bubbling)", "goAquaFlux", median(r))
add("field: share of bubbling incubations with that ratio > 2", "goAquaFlux", mean(r > 2))

## Events coinciding with a goAquaFlux event
ref <- evr %>% filter(method == "goAquaFlux") %>% select(id, t_event, magnitude_ppb)
MR  <- match_all(evr, ref, ids, intersect(c("FS_default", "MSP"), unique(evr$method)))
agree <- MR$per %>% group_by(method) %>% summarise(v = sum(TP) / sum(n_est), .groups = "drop")
for (k in seq_len(nrow(agree)))
  add("field: share of the method's events coinciding with a goAquaFlux event",
      agree$method[k], agree$v[k])

## Outcome of the ebullition check
if ("ebullition_check" %in% names(flr)) {
  ec <- flr$ebullition_check[flr$method == "goAquaFlux"]
  tab <- table(ifelse(is.na(ec), "NA", ifelse(ec == "", "passed", ec)))
  for (k in names(tab)) add(paste0("field: ebullition check = ", k), "goAquaFlux", tab[[k]])
}


## ============================================================================
## OUTPUT
## ============================================================================
numbers <- bind_rows(numbers)
write.csv(numbers, file.path("results", "Section5_numbers.csv"), row.names = FALSE)

get <- function(q, m) { v <- numbers$value[numbers$quantity == q & numbers$method == m]; if (length(v)) v[1] else NA }
fmt_err <- function(comp, m) {                      # "median error (median |error|)"
  e <- get(paste0("synthetic: median relative error, ", comp, " flux"), m)
  a <- get(paste0("synthetic: median absolute relative error, ", comp, " flux"), m)
  if (is.na(e)) "n/a" else paste0(pct_signed(e), " (", pct_plain(a), ")")
}
table3 <- data.frame(row = c(
  "Synthetic: total flux", "Synthetic: diffusive flux", "Synthetic: ebullitive flux",
  "Synthetic: no flux estimate", "Synthetic: bubble recall / precision",
  "Synthetic: events per bubble-free incubation", "Synthetic: bubble magnitude error",
  "Field: incubations with ebullition", "Field: no flux estimate"))
for (m in table_methods) {
  rec <- get("synthetic: bubble recall", m)
  ev0 <- get("synthetic: events per bubble-free incubation", m)
  me  <- get("synthetic: median relative error of bubble magnitude", m)
  table3[[m]] <- c(
    fmt_err("total", m), fmt_err("diffusive", m), fmt_err("ebullition", m),
    pct_plain(get("synthetic: share without a total flux", m), 1),
    if (is.na(rec)) "n/a" else sprintf("%.2f / %.2f", rec, get("synthetic: bubble precision", m)),
    if (is.na(ev0)) "n/a" else sprintf("%.2f", ev0),
    if (is.na(me)) "n/a" else paste0(pct_signed(me), " (", pct_plain(
      get("synthetic: median absolute relative error of bubble magnitude", m)), ")"),
    pct_plain(get("field: share of incubations with ebullition", m)),
    pct_plain(get("field: share without a total flux", m)))
}
## Conventional flux reported as a diffusive flux (footnote a of Table 3)
table3$goFlux_noSep[2] <- paste0(
  pct_signed(get("synthetic: conventional flux as diffusive, median relative error, 1 bubble(s)", "goFlux_noSep")),
  " to ",
  pct_signed(get("synthetic: conventional flux as diffusive, median relative error, 10 bubble(s)", "goFlux_noSep")),
  " (a)")

print(table3, row.names = FALSE)
write.csv(table3, file.path("results", "Table3_method_comparison.csv"), row.names = FALSE)
cat("\nWritten: results/Table3_method_comparison.csv and results/Section5_numbers.csv\n")
