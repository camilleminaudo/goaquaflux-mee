###############################################################################
## make_Fig4_goFlux_vs_goAquaFlux.R
##
## Field-scale comparison of CH4 and CO2 fluxes estimated with and without
## flux separation (manuscript Section 4, Figure 4):
##   * goFlux + best.flux on the whole incubation (conventional calculation,
##     no separation);
##   * goAquaFlux (bubbles detected on CH4; diffusive flux on a bubble-free
##     window, ebullitive flux from the bubble steps).
##
## Steps
##   1. for each incubation of data/many_incubations.RData: identify the
##      measurement window (autoID) and compute both estimates for CH4 and
##      CO2 (skipped when recompute = FALSE: saved results are reloaded);
##   2. summary numbers quoted in Section 4;
##   3. Figure 4: total CH4 flux (a) and diffusive CO2 flux (b), with vs
##      without separation, point size = ebullitive share of the CH4 flux.
##
## Incubations are skipped when they have no measurements, or when the start
## time of the auxiliary file is more than 5 s before the first measurement
## (clock mismatch, which would truncate the incubation). Skipped and failed
## incubations are logged with the reason.
##
## Inputs   data/many_incubations.RData
## Outputs  results/goFlux_vs_goAquaFlux/data_goFlux_vs_goAquaFlux.RData
##          results/goFlux_vs_goAquaFlux/Fig4_summary.csv  (Section 4 numbers)
##          results/figures/Fig4_goFlux_vs_goAquaFlux.jpeg / .svg
##
## Units: CH4 fluxes in nmol m-2 s-1, CO2 fluxes in umol m-2 s-1 (goFlux).
##
## Run from the repository root: Rscript R/make_Fig4_goFlux_vs_goAquaFlux.R
###############################################################################

source(file.path("R", "setup.R"))   # paths, goFlux, BUBBLE_WINDOW

suppressPackageStartupMessages({
  library(ggplot2)
  library(egg)       # theme_article()
  library(dplyr)
})

recompute   <- TRUE     # FALSE: reload the saved fluxes and only redo the figure
time_tol_s  <- 5        # max. lag of aux start time before the first measurement
out_dir     <- file.path(PATHS$results, "goFlux_vs_goAquaFlux")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
res_file    <- file.path(out_dir, "data_goFlux_vs_goAquaFlux.RData")


## ---- 1. Fluxes with and without separation -------------------------------------
if (recompute) {

  load(PATHS$real_data)                         # mydata_all, myauxfile
  myauxfile$obs.length <- myauxfile$duration    # observation length for autoID()

  CH4_goFlux <- CH4_goAquaFlux <- CO2_goFlux <- CO2_goAquaFlux <- NULL
  skipped <- NULL                               # incubations not processed, with reason

  list_ids <- sort(unique(myauxfile$UniqueID))
  for (i in list_ids) {
    message("processing ", i)
    aux_i  <- myauxfile[myauxfile$UniqueID == i, ]
    data_i <- mydata_all[mydata_all$UniqueID == i, ]

    if (nrow(data_i) == 0) {
      skipped <- rbind(skipped, data.frame(UniqueID = i, reason = "no measurements"))
      next
    }

    ## Clock check: skip if the auxiliary start time precedes the data
    time_diff <- as.numeric(difftime(min(data_i$POSIX.time), aux_i$start.time[1],
                                     units = "secs"))
    if (time_diff < -time_tol_s) {
      skipped <- rbind(skipped, data.frame(
        UniqueID = i, reason = sprintf("aux start time %.0f s before data", -time_diff)))
      next
    }

    IDed <- autoID(inputfile = data_i, auxfile = aux_i, shoulder = 0)

    tryCatch({
      ## CO2: no separation, and separation with bubbles detected on CH4
      co2_nosep <- best.flux(goFlux(dataframe = IDed, gastype = "CO2dry_ppm",
                                    H2O_col = "H2O_ppm"))
      co2_sep   <- goAquaFlux(dataframe = IDed, gastype = "CO2dry_ppm",
                              bubble_gas = "CH4dry_ppb", use_bubble_detection = TRUE,
                              H2O_col = "H2O_ppm", bubble.method = "diff",
                              bubble.window.size = BUBBLE_WINDOW)
      ## CH4: no separation, and separation
      ch4_nosep <- best.flux(goFlux(dataframe = IDed, gastype = "CH4dry_ppb",
                                    H2O_col = "H2O_ppm"))
      ch4_sep   <- goAquaFlux(dataframe = IDed, gastype = "CH4dry_ppb",
                              bubble_gas = "CH4dry_ppb", use_bubble_detection = TRUE,
                              H2O_col = "H2O_ppm", bubble.method = "diff",
                              bubble.window.size = BUBBLE_WINDOW)

      ## Results are appended only when all four estimates succeeded, so that
      ## the four tables always contain the same incubations.
      CO2_goFlux     <- rbind(CO2_goFlux, co2_nosep)
      CO2_goAquaFlux <- rbind(CO2_goAquaFlux, co2_sep$flux_summary)
      CH4_goFlux     <- rbind(CH4_goFlux, ch4_nosep)
      CH4_goAquaFlux <- rbind(CH4_goAquaFlux, ch4_sep$flux_summary)
    }, error = function(e) {
      skipped <<- rbind(skipped, data.frame(UniqueID = i, reason = conditionMessage(e)))
      warning("Error processing ", i, ": ", conditionMessage(e))
    })
  }

  message("\nProcessed: ", nrow(CH4_goAquaFlux), " incubations; skipped or failed: ",
          if (is.null(skipped)) 0 else nrow(skipped))
  save(CH4_goFlux, CH4_goAquaFlux, CO2_goFlux, CO2_goAquaFlux, skipped, file = res_file)
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

} else {
  load(res_file)
}


## ---- 2. Merge and summary numbers (Section 4) -------------------------------------
## Estimates are matched by UniqueID.
ch4 <- CH4_goAquaFlux %>%
  mutate(flux_no_separation = CH4_goFlux$best.flux[match(UniqueID, CH4_goFlux$UniqueID)],
         bubbling = !is.na(flux_ebullition) & flux_ebullition > 0,
         ## ebullitive share of the total CH4 flux; values outside [0, 1]
         ## (e.g. negative diffusive flux) are not interpretable and set to NA
         ebullition_contribution = flux_ebullition / flux_total,
         ebullition_contribution = ifelse(ebullition_contribution < 0 |
                                            ebullition_contribution > 1,
                                          NA, ebullition_contribution))

co2 <- CO2_goAquaFlux %>%
  mutate(flux_no_separation = CO2_goFlux$best.flux[match(UniqueID, CO2_goFlux$UniqueID)],
         ebullition_contribution = ch4$ebullition_contribution[match(UniqueID, ch4$UniqueID)],
         bubbling = ch4$bubbling[match(UniqueID, ch4$UniqueID)])

b <- ch4$bubbling
ratio_ch4 <- ch4$flux_no_separation[b] / ch4$flux_diffusive[b]      # conventional / diffusive
ratio_ch4 <- ratio_ch4[is.finite(ratio_ch4) & ratio_ch4 > 0]
ratio_co2 <- abs(co2$flux_no_separation[co2$bubbling %in% TRUE]) /
  abs(co2$flux_diffusive[co2$bubbling %in% TRUE])
ratio_co2 <- ratio_co2[is.finite(ratio_co2)]

summary_s4 <- data.frame(
  quantity = c(
    "incubations processed",
    "incubations skipped or failed",
    "incubations with >= 1 detected bubble",
    "share of incubations with >= 1 detected bubble",
    "incubations with ebullition > 10% of the total CH4 flux",
    "median ebullitive share of the CH4 flux (bubbling incubations)",
    "ebullitive share of the CH4 flux summed over the dataset",
    "median ratio conventional / goAquaFlux diffusive CH4 flux (bubbling)",
    "share of bubbling incubations with that ratio > 2",
    "median ratio |conventional| / |goAquaFlux diffusive| CO2 flux (bubbling)"),
  value = c(
    nrow(ch4),
    if (is.null(skipped)) 0 else nrow(skipped),
    sum(b),
    mean(b),
    sum(ch4$ebullition_contribution > 0.1, na.rm = TRUE),
    median(ch4$ebullition_contribution[b], na.rm = TRUE),
    sum(ch4$flux_ebullition, na.rm = TRUE) / sum(ch4$flux_total, na.rm = TRUE),
    median(ratio_ch4),
    mean(ratio_ch4 > 2),
    median(ratio_co2)))
print(summary_s4, row.names = FALSE)
write.csv(summary_s4, file.path(out_dir, "Fig4_summary.csv"), row.names = FALSE)


## ---- 3. Figure 4 ------------------------------------------------------------------------
lab_s <- expression(atop("Ebullition ratio",
                         F[CH[4]~ebull]^goAquaFlux*" / "*F[CH[4]~tot]^goAquaFlux))

theme_fig4 <- theme_article() +
  theme(legend.position      = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background    = element_rect(fill = alpha("white", 0.75), colour = NA),
        legend.title         = element_text(size = 9, lineheight = 0.9),
        legend.text          = element_text(size = 8),
        legend.key.height    = unit(0.9, "lines"),
        axis.title           = element_text(lineheight = 1.0))

## (a) CH4: total flux with vs without separation. Log axes: non-positive
## fluxes are not shown.
plt_ch4 <- ggplot(ch4, aes(flux_no_separation, flux_total, size = ebullition_contribution)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_point(alpha = 0.75, shape = 21, fill = "#F2B876", colour = "#E29338", stroke = 0.5) +
  scale_x_log10() + scale_y_log10() +
  scale_size_continuous(name = lab_s, range = c(1.5, 4),
                        labels = scales::percent_format(accuracy = 1)) +
  coord_fixed() +
  labs(x = expression(atop(F[CH[4]~tot]^goFlux*" [nmol "*m^-2*" "*s^-1*"]",
                           "no flux separation")),
       y = expression(atop(F[CH[4]~tot]^goAquaFlux*" [nmol "*m^-2*" "*s^-1*"]",
                           "with flux separation"))) +
  theme_fig4

## (b) CO2: diffusive flux (goAquaFlux) vs conventional flux. CO2 fluxes can
## be negative (uptake); incubations with a total CO2 flux <= 1e-5 are not
## shown, and absolute values are plotted on log axes.
plt_co2 <- ggplot(co2[!is.na(co2$flux_total) & co2$flux_total > 1e-5, ],
                  aes(abs(flux_no_separation), abs(flux_diffusive),
                      size = ebullition_contribution)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_point(alpha = 0.75, shape = 21, fill = "#1C7293", colour = "#0D4861", stroke = 0.5) +
  scale_x_log10() + scale_y_log10() +
  scale_size_continuous(name = lab_s, range = c(1.5, 4),
                        labels = scales::percent_format(accuracy = 1)) +
  coord_fixed() +
  labs(x = expression(atop(F[CO[2]]^goFlux*" ["*mu*"mol "*m^-2*" "*s^-1*"]",
                           "no flux separation")),
       y = expression(atop(F[CO[2]]^goAquaFlux*" ["*mu*"mol "*m^-2*" "*s^-1*"]",
                           "with flux separation"))) +
  theme_fig4

fig4 <- ggpubr::ggarrange(plt_ch4, plt_co2, nrow = 1, labels = c("a", "b"),
                          align = "h", common.legend = TRUE, legend = "bottom")

ggsave(plot = fig4, filename = "Fig4_goFlux_vs_goAquaFlux.jpeg", path = PATHS$figures,
       width = 6, height = 3.2, units = "in", scale = 1.2, dpi = 300)
ggsave(plot = fig4, filename = "Fig4_goFlux_vs_goAquaFlux.svg", path = PATHS$figures,
       width = 6, height = 3.2, units = "in", scale = 1.2)
cat("Figure 4 written to", PATHS$figures, "\n")
