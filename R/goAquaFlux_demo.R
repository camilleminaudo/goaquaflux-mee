
# ============================================================================
# SCRIPT:
# ============================================================================
#
# DESCRIPTION:
# This script ...
#
# ============================================================================

# Clear workspace and console
rm(list = ls())
cat("\014")

# ---- Loading libraries and functions ----

library(ggplot2)     # For plotting (ggplot, geom_point, scale_*_log10, etc.)
library(egg)         # For theme_article()
library(goFlux)      # For autoID, goFlux, goAquaFlux, best.flux

repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))

devtools::load_all("C:/Projects/myGit/goFlux")   # path to the package root


# ---- Setting paths to directories ----
data_path <- paste0(repo_root,"/data")
results_path <- paste0(repo_root,"/results")


# ---- Loading data ----
setwd(data_path)
load("many_incubations.RData")

myauxfile$obs.length <- myauxfile$duration

# ---- selected incubations ----

selct <- c(
  # "s1-da-p1-8-o-d-10:05",
  # "s1-cu-a2-16-o-d-11:58",
  # "s3-ca-r1-1-o-d-07:46",
  # "s2-ri-a2-15-o-d-11:56",
  "s1-cu-a2-1-o-d-06:36"   # chosen example for publication
)

myaux <- myauxfile[which(myauxfile$UniqueID %in% selct),]

# Automatic selection of data
IDed <- autoID(inputfile = mydata_all, auxfile = myaux, shoulder = 0, deadband = 0, crop.end = 0)
print(IDed)



IDed$CH4_prec <- max(IDed$CH4_prec)
IDed$CO2_prec <- max(IDed$CO2_prec)

flux_sep_ch4 <- goAquaFlux(dataframe = IDed,
                           gastype = "CH4dry_ppb",
                           use_bubble_detection = T,
                           bubble.method = "diff",
                           bubble_gas = "CH4dry_ppb")

flux_sep_co2 <- goAquaFlux(dataframe = IDed, gastype = "CO2dry_ppm",
                           use_bubble_detection = T, bubble.method = "diff", bubble_gas = "CH4dry_ppb")

# Create plots
p_sep_ch4 <- flux.plot.aqua(
  flux.results = flux_sep_ch4,
  dataframe = IDed,
  gastype = "CH4dry_ppb")


p_sep_co2 <- flux.plot.aqua(
  flux.results = flux_sep_co2,
  dataframe = IDed,
  gastype = "CO2dry_ppm")





flux_NO_sep_ch4 <- goAquaFlux(dataframe = IDed, gastype = "CH4dry_ppb",
                              use_bubble_detection = F, bubble.method = "diff", bubble_gas = "CH4dry_ppb")

flux_NO_sep_co2 <- goAquaFlux(dataframe = IDed, gastype = "CO2dry_ppm",
                              use_bubble_detection = F, bubble.method = "diff", bubble_gas = "CH4dry_ppb")

# Create plots
p_NO_sep_ch4 <- flux.plot.aqua(
  flux.results = flux_NO_sep_ch4,
  dataframe = IDed,
  gastype = "CH4dry_ppb",
  plot.display = NULL)

p_NO_sep_co2 <- flux.plot.aqua(
  flux.results = flux_NO_sep_co2,
  dataframe = IDed,
  gastype = "CO2dry_ppm",
  plot.display = NULL)



p_sep_ch4_notitle <- lapply(p_sep_ch4, function(p) p + labs(title = NULL))
p_NO_sep_ch4_notitle <- lapply(p_NO_sep_ch4, function(p) p + labs(title = NULL))
p_sep_co2_notitle <- lapply(p_sep_co2, function(p) p + labs(title = NULL))
p_NO_sep_co2_notitle <- lapply(p_NO_sep_co2, function(p) p + labs(title = NULL))






ggsave(plot = p_sep_ch4_notitle[1], filename = "Fig2a_chosen_incubation.svg", path = results_path,
       width = 5, height = 3.5, dpi = 300, units = 'in', scale = 1)


ggsave(plot = p_sep_co2_notitle[1], filename = "Fig2b_chosen_incubation.svg", path = results_path,
       width = 5, height = 3.5, dpi = 300, units = 'in', scale = 1)



library(patchwork)
setwd(repo_root)
source("R/combine_aqua_plots.R")

fig <- combine_aqua_plots(list(p_sep_ch4_notitle[[1]], p_sep_co2_notitle[[1]]),
                          ncol = 2)



# the legend counts as a third element, so give it an empty tag
fig2 <- fig + plot_annotation(tag_levels = list(c("a", "b", "")))




ggsave(plot = fig2, filename = "Fig2_chosen_incubation.svg", path = results_path,
       width = 8, height = 3.5, dpi = 300, units = 'in', scale = 1)

ggsave(plot = fig2, filename = "Fig2_chosen_incubation.jpeg", path = results_path,
       width = 8, height = 3.5, dpi = 300, units = 'in', scale = 1)




fig <- wrap_plots(c(p_NO_sep_ch4_notitle[1],p_sep_ch4_notitle[1],
                    p_NO_sep_co2_notitle[1],p_sep_co2_notitle[1]), ncol = 2) +
  plot_layout(guides = "collect", tag_level = "new") &
  theme(legend.position = "bottom")

fig2 <- fig + plot_annotation(tag_levels = "a")

