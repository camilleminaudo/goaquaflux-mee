
# ============================================================================
# SCRIPT: Comparison of goFlux vs goAquaFlux for CH4 and CO2 flux estimates
# ============================================================================
#
# DESCRIPTION:
# This script processes a batch of aquatic chamber incubations to quantify
# greenhouse gas (CH4 and CO2) fluxes using two methodological approaches:
#   - goFlux: traditional flux calculation (no flux separation)
#   - goAquaFlux: flux separation method accounting for ebullitive pathways
#
# The analysis compares flux estimates between methods, calculates the 
# contribution of ebullition to total CH4 flux, and generates publication-ready
# comparison plots. Includes robust error handling for missing or problematic
# data files, and tracks failed incubations for quality control.
#
# OUTPUT:
#   - CH4_goFlux, CH4_goAquaFlux: CH4 flux results (both methods)
#   - CO2_goFlux, CO2_goAquaFlux: CO2 flux results (both methods)
#   - failed_incubations: log of incubations that failed processing
#   - df.no_measurements: log of incubations lacking required data
#   - Fig3_goFlux_vs_goAquaFlux_R4Cs.jpeg: comparison figure
#
# NOTES:
#   - Requires goFlux package with autoID, goFlux, goAquaFlux, best.flux functions
#   - Tolerance for timing mismatches between auxiliary file and measurements: ±5 sec
#   - Ebullition ratios <0 or >1 are treated as invalid and set to NA
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


compute_it = F


# ---- Loading data ----
setwd(data_path)
load("many_incubations.RData")

myauxfile$obs.length <- myauxfile$duration

# ---- Running for all data available ----

if (compute_it){
  
  CH4_goFlux <- CH4_goAquaFlux <- CO2_goFlux <- CO2_goAquaFlux <- df.no_measurements <- NULL
  failed_incubations <- NULL  # Track incubations with timing issues
  
  list_ids <- sort(unique(myauxfile$UniqueID))
  for (k in seq_along(list_ids)){
    i = list_ids[k]
    message(paste0("processing ", i))
    
    myauxfile_i <- myauxfile[which(myauxfile$UniqueID==i),]
    
    if(dim(myauxfile_i)[1]==0){
      message(paste0("Could not find corresponding myauxfile for ", i))
      df.no_measurements <- rbind(df.no_measurements,
                                  data.frame(UniqueID=i,
                                             message = "no myauxfile"))
    } else {
      mydata <- mydata_all[which(mydata_all$UniqueID==i),]
      
      if(dim(mydata)[1]==0){
        df.no_measurements <- rbind(df.no_measurements,
                                    data.frame(UniqueID=i,
                                               message = "no measurements"))
      } else {
        
        # ---- Check time tolerance ----
        data_start_time <- min(mydata$POSIX.time)
        aux_start_time <- myauxfile_i$start.time[1]
        time_diff <- as.numeric(difftime(data_start_time, aux_start_time, units = "secs"))
        
        if(time_diff < -5){  # If aux start time is more than 5 secs before data
          message(paste0("WARNING: Auxiliary start time is ", abs(time_diff), 
                         " seconds before first measurement. Skipping incubation ", i))
          failed_incubations <- rbind(failed_incubations,
                                      data.frame(UniqueID = i,
                                                 aux_start_time = aux_start_time,
                                                 data_start_time = data_start_time,
                                                 time_diff_secs = time_diff,
                                                 message = "start time before data range"))
          next  # Skip to next iteration
        }
        
        # ---- Process normally ----
        IDed <- autoID(inputfile = mydata, auxfile = myauxfile_i, shoulder = 0)
        
        tryCatch({
          CO2_goFlux_i <- goFlux(dataframe = IDed, gastype = "CO2dry_ppm", H2O_col = "H2O_ppm")
          CO2_goFlux_i <- best.flux(CO2_goFlux_i)
          
          CO2_goAquaFlux_i <- goAquaFlux(dataframe = IDed, gastype = "CO2dry_ppm", bubble_gas = "CH4dry_ppb",
                                         use_bubble_detection = TRUE, H2O_col = "H2O_ppm", 
                                         bubble.method = "diff")
          
          CO2_goFlux <- rbind(CO2_goFlux, CO2_goFlux_i)
          CO2_goAquaFlux <- rbind(CO2_goAquaFlux, CO2_goAquaFlux_i$flux_summary)
          
          CH4_goFlux_i <- goFlux(dataframe = IDed, gastype = "CH4dry_ppb", H2O_col = "H2O_ppm")
          CH4_goFlux_i <- best.flux(CH4_goFlux_i)
          
          CH4_goAquaFlux_i <- goAquaFlux(dataframe = IDed, gastype = "CH4dry_ppb", 
                                         use_bubble_detection = TRUE, H2O_col = "H2O_ppm", 
                                         bubble.method = "diff")
          
          CH4_goFlux <- rbind(CH4_goFlux, CH4_goFlux_i)
          CH4_goAquaFlux <- rbind(CH4_goAquaFlux, CH4_goAquaFlux_i$flux_summary)
          
        }, error = function(e){
          failed_incubations <<- rbind(failed_incubations,
                                       data.frame(UniqueID = i,
                                                  aux_start_time = NA,
                                                  data_start_time = NA,
                                                  time_diff_secs = NA,
                                                  message = conditionMessage(e)))
          warning("Error processing ", i, ": ", conditionMessage(e))
        })
      }
    }
  }
  
  # ---- Summary ----
  message("\n--- Processing Summary ---")
  if(!is.null(failed_incubations)){
    message("Failed incubations: ", nrow(failed_incubations))
    print(failed_incubations)
  }
  
  
  setwd(dir = results_path)
  save(list = c("df.no_measurements", "failed_incubations",
                "CH4_goFlux", "CH4_goAquaFlux", "CO2_goFlux", "CO2_goAquaFlux"),
       file = "data_goFlux_vs_goAquaFlux.RData")
  
} else {
  setwd(dir = results_path)
  load("data_goFlux_vs_goAquaFlux.RData")
}



CH4_goAquaFlux$best.flux_no_separation = CH4_goFlux$best.flux
CH4_goAquaFlux$LM.flux_no_separation = CH4_goFlux$LM.flux

# CH4_goAquaFlux$UniqueID[which(CH4_goAquaFlux$flux_total<0)]


CH4_goAquaFlux_sel <- CH4_goAquaFlux#[CH4_goAquaFlux$flux_total>0 & CH4_goAquaFlux$flux_diffusive>0,]
# CH4_goAquaFlux_sel <- CH4_goAquaFlux_sel[!is.na(CH4_goAquaFlux_sel$flux_total),]

CH4_goAquaFlux_sel$ebullition_contribution <- CH4_goAquaFlux_sel$flux_ebullition/CH4_goAquaFlux_sel$flux_total
CH4_goAquaFlux_sel$ebullition_contribution[CH4_goAquaFlux_sel$ebullition_contribution<0] <- NA
CH4_goAquaFlux_sel$ebullition_contribution[CH4_goAquaFlux_sel$ebullition_contribution>1] <- NA


summary(CH4_goAquaFlux_sel$ebullition_contribution )
summary(CH4_goAquaFlux_sel$flux_diffusive )



lab_x <- expression(atop(F[CH[4]~tot]^goFlux*" [nmol "*m^-2*" "*s^-1*"]",
                         "no flux separation"))
lab_y <- expression(atop(F[CH[4]~tot]^goAquaFlux*" [nmol "*m^-2*" "*s^-1*"]",
                         "with flux separation"))
lab_s <- expression(atop("Ebullition ratio",
                         F[CH[4]~ebull]^goAquaFlux*" / "*F[CH[4]~tot]^goAquaFlux))

plt_ch4 <- ggplot(CH4_goAquaFlux_sel,
                  aes(best.flux_no_separation, flux_total, size = ebullition_contribution)) +
  geom_abline(slope = 1, intercept = 0,
              # linetype = "dashed",
              colour = "grey40", linewidth = 0.4) +
  geom_point(alpha = 0.75, shape = 21, fill = "#F2B876", colour = "#E29338",
             stroke = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  scale_size_continuous(name = lab_s, range = c(1.5, 4),
                        labels = scales::percent_format(accuracy = 1)) +
  coord_fixed() +
  labs(x = lab_x, y = lab_y) +
  theme_article() +
  theme(
    legend.position      = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background    = element_rect(fill = alpha("white", 0.75), colour = NA),
    legend.title         = element_text(size = 9, lineheight = 0.9),
    legend.text          = element_text(size = 8),
    legend.key.height    = unit(0.9, "lines"),
    axis.title           = element_text(lineheight = 1.0)
  )


# how many incubations show ebullition ratio > 10%?
length(which(CH4_goAquaFlux_sel$ebullition_contribution>0.1))
length(CH4_goAquaFlux_sel$UniqueID)






# ------- CO2 flux



CO2_goAquaFlux$best.flux_no_separation = CO2_goFlux$best.flux
CO2_goAquaFlux$LM.flux_no_separation = CO2_goFlux$LM.flux

# CO2_goAquaFlux$UniqueID[which(CO2_goAquaFlux$flux_total<0)]


CO2_goAquaFlux_sel <- CO2_goAquaFlux#[abs(CO2_goAquaFlux$flux_diffusive)>1e-5,]
# CO2_goAquaFlux_sel <- CO2_goAquaFlux_sel[!is.na(CO2_goAquaFlux_sel$flux_total),]

CO2_goAquaFlux_sel$ebullition_contribution <- CH4_goAquaFlux_sel$ebullition_contribution[match(CO2_goAquaFlux_sel$UniqueID, 
                                                                                               CH4_goAquaFlux_sel$UniqueID)]

lab_x <- expression(atop(F[CO[2]]^goFlux*" [mmol "*m^-2*" "*s^-1*"]",
                         "no flux separation"))
lab_y <- expression(atop(F[CO[2]]^goAquaFlux*" [mmol "*m^-2*" "*s^-1*"]",
                         "with flux separation"))
lab_s <- expression(atop("Ebullition ratio",
                         F[CH[4]~ebull]^goAquaFlux*" / "*F[CH[4]~tot]^goAquaFlux))

plt_co2 <- ggplot(CO2_goAquaFlux_sel[CO2_goAquaFlux_sel$flux_total>1e-5,],
                  aes(abs(best.flux_no_separation), abs(flux_diffusive), size = ebullition_contribution)) +
  geom_abline(slope = 1, intercept = 0,
              # linetype = "dashed",
              colour = "grey40", linewidth = 0.4) +
  geom_point(alpha = 0.75, shape = 21, fill = "#1C7293", colour = "#0D4861",
             stroke = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  scale_size_continuous(name = lab_s, range = c(1.5, 4),
                        labels = scales::percent_format(accuracy = 1)) +
  coord_fixed() +
  labs(x = lab_x, y = lab_y) +
  theme_article() +
  theme(
    legend.position      = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background    = element_rect(fill = alpha("white", 0.75), colour = NA),
    legend.title         = element_text(size = 9, lineheight = 0.9),
    legend.text          = element_text(size = 8),
    legend.key.height    = unit(0.9, "lines"),
    axis.title           = element_text(lineheight = 1.0)
  )


# how many incubations show ebullition ratio > 10%?
length(which(CO2_goAquaFlux_sel$ebullition_contribution>0.1))
length(CO2_goAquaFlux_sel$UniqueID)



plt <- ggpubr::ggarrange( plt_ch4, plt_co2,
                         nrow = 1, labels = c("a","b"), 
                         align = "h", common.legend = T, legend = "bottom")

plt

ggsave(plot = plt, filename = "Fig3_goFlux_vs_goAquaFlux_R4Cs.jpeg", path = results_path,
       width = 6, height = 3.2, dpi = 300, units = 'in', scale = 1.2)





# -----------------


library(dplyr)
library(ggplot2)
library(egg)

d_ch4 <- CH4_goAquaFlux_sel %>%
  # filter(best.flux_no_separation > 0, flux_total > 0) %>%
  mutate(
    ratio = flux_total / best.flux_no_separation,
    # keep the exact zeros as their own group; collapse the sparse upper tail
    eb_bin = cut(ebullition_contribution,
                 breaks = c(-Inf, 0, 1),
                 labels = c("0\n(none)", ">0"),
                 # breaks = c(-Inf, 0, 0.25, 0.5, 0.75, 1),
                 # labels = c("0\n(none)", "0-0.25", "0.25-0.5", "0.5-0.75",
                 #            "0.75-1"),
                 include.lowest = TRUE),
  )

n_lab <- d_ch4 %>% count(eb_bin) %>% mutate(lab = paste0("n = ", n))

plt_ratios_ch4 <- ggplot(d_ch4, aes(eb_bin, ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_boxplot(outlier.shape = NA, width = 0.55,
               fill = "#F2B876", colour = "grey30", linewidth = 0.35) +
  geom_jitter(width = 0.14, height = 0,
              alpha = 0.6, size = 1.6) +
  geom_text(data = n_lab, aes(x = eb_bin, y = Inf, label = lab),
            vjust = 1.4, size = 2.8, colour = "grey35", inherit.aes = FALSE) +
  scale_y_log10() +
  labs(
    x = expression(atop("Ebullition ratio",
                        F[CH[4]~ebull]^goAquaFlux*" / "*F[CH[4]~tot]^goAquaFlux)),
    y = expression(atop(F[CH[4]~tot]^goAquaFlux*" / "*F[CH[4]~tot]^goFlux,
                        "(flux ratio, with / without separation)"))
  ) +
  theme_article() +
  theme(axis.title = element_text(lineheight = 1.0),
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.title = element_text(size = 8, lineheight = 0.9),
        legend.text = element_text(size = 8))




d_co2 <- CO2_goAquaFlux_sel %>%
  mutate(
    ratio = abs(flux_diffusive) / abs(best.flux_no_separation),
    # keep the exact zeros as their own group; collapse the sparse upper tail
    eb_bin = cut(ebullition_contribution,
                 breaks = c(-Inf, 0, 1),
                 labels = c("0\n(none)", ">0"),
                 include.lowest = TRUE)
  )

n_lab_co2 <- d_co2 %>% count(eb_bin) %>% mutate(lab = paste0("n = ", n))

plt_ratios_co2 <- ggplot(d_co2, aes(eb_bin, ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_boxplot(outlier.shape = NA, width = 0.55,
               fill = "#E8F0F4", colour = "grey30", linewidth = 0.35) +
  geom_jitter(width = 0.14, height = 0,
              alpha = 0.6, size = 1.6) +
  geom_text(data = n_lab_co2, aes(x = eb_bin, y = Inf, label = lab),
            vjust = 1.4, size = 2.8, colour = "grey35", inherit.aes = FALSE) +
  scale_y_log10() +
  labs(
    x = expression(atop("Ebullition ratio",
                        F[CH[4]~ebull]^goAquaFlux*" / "*F[CH[4]~tot]^goAquaFlux)),
    y = expression(atop(F[CO2[2]]^goAquaFlux*" / "*F[CO2[2]]^goFlux,
                        "(flux ratio, with / without separation)"))
  ) +
  theme_article() +
  theme(axis.title = element_text(lineheight = 1.0),
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.title = element_text(size = 8, lineheight = 0.9),
        legend.text = element_text(size = 8))



plt <- ggpubr::ggarrange(plt_ratios_co2, plt_ratios_ch4,
                         nrow = 1, labels = c("a","b"), 
                         align = "h", common.legend = T, legend = "bottom")

plt
