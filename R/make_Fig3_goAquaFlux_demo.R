###############################################################################
## make_Fig3_goAquaFlux_demo.R
##
## Worked example of the manuscript (Section 4, Figure 3): one field
## incubation with two ebullition events, processed with goAquaFlux for CH4
## (the bubble gas) and CO2 (bubbles detected on CH4), and compared with
## conventional goFlux flux calculation without separation.
##
## Steps
##   1. load the field data and select the example incubation;
##   2. identify the measurement window (goFlux::autoID);
##   3. goAquaFlux for CH4 and CO2 (flux separation);
##   4. the same without separation (goFlux LM/HM + best.flux on the whole
##      incubation), for comparison;
##   5. Figure 3: flux.plot.aqua() panels for CH4 and CO2, combined with one
##      shared legend (combine_aqua_plots()).
##
## Inputs   data/many_incubations.RData
## Outputs  results/figures/Fig3_worked_example.svg / .jpeg
##          results/Fig3_worked_example_fluxes.csv  (fluxes quoted in Section 4)
##
## Run from the repository root: Rscript R/make_Fig3_goAquaFlux_demo.R
###############################################################################

source(file.path("R", "setup.R"))                # paths, goFlux, BUBBLE_WINDOW
source(file.path("R", "combine_aqua_plots.R"))   # one legend for several panels

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})


## ---- 1. Data and example incubation ----------------------------------------------
load(PATHS$real_data)                         # mydata_all, myauxfile
myauxfile$obs.length <- myauxfile$duration    # observation length used by autoID()

example_id <- "s1-cu-a2-1-o-d-06:36"          # incubation shown in Figure 3
myaux <- myauxfile[myauxfile$UniqueID %in% example_id, ]


## ---- 2. Measurement window ----------------------------------------------------------
IDed <- autoID(inputfile = mydata_all, auxfile = myaux,
               shoulder = 0, deadband = 0, crop.end = 0)

## One precision value per gas for the whole incubation: the largest value
## reported by the analyser over the incubation (conservative minimum
## detectable flux).
IDed$CH4_prec <- max(IDed$CH4_prec, na.rm = TRUE)
IDed$CO2_prec <- max(IDed$CO2_prec, na.rm = TRUE)


## ---- 3. Flux separation with goAquaFlux ---------------------------------------------
## Bubbles are detected on CH4 for both gases.
flux_sep_ch4 <- goAquaFlux(dataframe = IDed, gastype = "CH4dry_ppb",
                           use_bubble_detection = TRUE, bubble_gas = "CH4dry_ppb",
                           bubble.method = "diff", bubble.window.size = BUBBLE_WINDOW)

flux_sep_co2 <- goAquaFlux(dataframe = IDed, gastype = "CO2dry_ppm",
                           use_bubble_detection = TRUE, bubble_gas = "CH4dry_ppb",
                           bubble.method = "diff", bubble.window.size = BUBBLE_WINDOW)


## ---- 4. Conventional flux calculation (no separation) ---------------------------------
## goAquaFlux with bubble detection turned off fits goFlux (LM and HM) and
## best.flux on the whole incubation, i.e. the conventional calculation.
flux_nosep_ch4 <- goAquaFlux(dataframe = IDed, gastype = "CH4dry_ppb",
                             use_bubble_detection = FALSE, bubble_gas = "CH4dry_ppb")
flux_nosep_co2 <- goAquaFlux(dataframe = IDed, gastype = "CO2dry_ppm",
                             use_bubble_detection = FALSE, bubble_gas = "CH4dry_ppb")

summarise_fluxes <- function(res, gas, method) {
  s <- res$flux_summary
  data.frame(gas = gas, method = method,
             model = if ("model" %in% names(s)) s$model[1] else NA,
             flux_total = s$flux_total[1], SE_total = s$SE_total[1],
             flux_diffusive = s$flux_diffusive[1], SE_diffusive = s$SE_diffusive[1],
             flux_ebullition = s$flux_ebullition[1], SE_ebullition = s$SE_ebullition[1])
}

fluxes <- rbind(summarise_fluxes(flux_sep_ch4,   "CH4dry_ppb", "goAquaFlux (separation)"),
                summarise_fluxes(flux_nosep_ch4, "CH4dry_ppb", "no separation"),
                summarise_fluxes(flux_sep_co2,   "CO2dry_ppm", "goAquaFlux (separation)"),
                summarise_fluxes(flux_nosep_co2, "CO2dry_ppm", "no separation"))
## Units: CH4 in nmol m-2 s-1, CO2 in umol m-2 s-1
print(fluxes, row.names = FALSE)
write.csv(fluxes, file.path(PATHS$results, "Fig3_worked_example_fluxes.csv"),
          row.names = FALSE)


## ---- 5. Figure 3 ------------------------------------------------------------------------
p_ch4 <- flux.plot.aqua(flux.results.ls = flux_sep_ch4, dataframe = IDed,
                        gastype = "CH4dry_ppb")
p_co2 <- flux.plot.aqua(flux.results.ls = flux_sep_co2, dataframe = IDed,
                        gastype = "CO2dry_ppm")

## The incubation identifier is given in the caption of the manuscript, so the
## panel titles are removed.
no_title <- function(p) p + labs(title = NULL)

## One legend key for the diffusive fit (each panel caption names its model)
fig <- combine_aqua_plots(list(no_title(p_ch4[[1]]), no_title(p_co2[[1]])),
                          ncol = 2, merge_fit_labels = TRUE)

## The legend is the third element of the layout: it gets an empty tag.
fig3 <- fig + plot_annotation(tag_levels = list(c("a", "b", "")))

## Full page width (7 in); `scale` enlarges the canvas so that panel headers
## and the legend fit without truncation.
ggsave(plot = fig3, filename = "Fig3_worked_example.svg", path = PATHS$figures,
       width = 7, height = 3.2, units = "in", scale = 1.45)
ggsave(plot = fig3, filename = "Fig3_worked_example.jpeg", path = PATHS$figures,
       width = 7, height = 3.2, units = "in", scale = 1.45, dpi = 300)

cat("Figure 3 written to", PATHS$figures, "\n")
