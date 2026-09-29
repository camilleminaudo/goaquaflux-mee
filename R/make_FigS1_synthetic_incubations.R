###############################################################################
## make_FigS1_synthetic_incubations.R
##
## Figure S1 (Appendix S1): construction of the synthetic incubations used in
## the benchmark of compare_goAquaFlux_FluxSeparator_MSP.R.
##   (a) response of the headspace CH4 concentration to one bubble, without
##       and with a transient overshoot (Eq. S4);
##   (b) noise-free incubation: diffusive accumulation + four bubbles;
##   (c) the same incubation with the three noise levels of the design.
##
## Inputs   R/synthetic_incubations.R (bubble_shape())
## Output   results/figures/FigS1_synthetic_incubations.png
##
## Run from the repository root: Rscript R/make_FigS1_synthetic_incubations.R
###############################################################################

if (!file.exists(file.path("R", "synthetic_incubations.R"))) {
  stop("Please run this script from the repository root.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})
source(file.path("R", "synthetic_incubations.R"))   # bubble_shape()

th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 9))

## Constants shared with the benchmark (CFG$syn in compare_goAquaFlux_FluxSeparator_MSP.R)
s <- list(duration = 600, dt = 1, C0 = 2000, ramp = 3, tau = 15)
## ---- (a) Response to one bubble, without and with overshoot -------------------
u <- seq(-10, 80, by = 0.25)
da <- rbind(data.frame(u = u, b = bubble_shape(u, 300, 0, s$tau, s$ramp), case = "overshoot = 0"),
            data.frame(u = u, b = bubble_shape(u, 300, 300, s$tau, s$ramp), case = "overshoot = 1"))
pa <- ggplot(da, aes(u, b, colour = case)) +
  geom_hline(yintercept = 300, linetype = 3, colour = "grey40") +
  geom_line(linewidth = 0.7) +
  annotate("text", x = 50, y = 270, label = "settled step M", size = 2.8, colour = "grey30") +
  annotate("segment", x = 0, xend = 3, y = 630, yend = 630, colour = "grey30",
           arrow = arrow(ends = "both", length = unit(1.2, "mm"))) +
  annotate("text", x = 12, y = 630, label = "ramp r = 3 s", size = 2.8, hjust = 0, colour = "grey30") +
  annotate("text", x = 22, y = 470, label = "decay of A,  tau = 15 s", size = 2.8, hjust = 0, colour = "grey30") +
  scale_colour_manual(values = c("#0D4861", "#E29338"), name = NULL) +
  labs(title = "(a) Response to one bubble (M = 300 ppb)", x = "Time since bubble onset (s)",
       y = "CH4 increase (ppb)") + th + theme(legend.position = c(0.78, 0.25),
                                            legend.background = element_blank())
## ---- (b) Noise-free incubation with four bubbles ----------------------------------
## Bubble times and sizes are fixed for illustration; in the benchmark they are
## random (R/synthetic_incubations.R).
set.seed(12)
t <- seq(0, s$duration - s$dt, by = s$dt)
tb <- c(118, 243, 371, 468); M <- c(210, 420, 160, 330); ov <- 1; slope <- 1
base <- s$C0 + slope * t
y0 <- base
for (i in seq_along(tb)) y0 <- y0 + bubble_shape(t - tb[i], M[i], ov * M[i], s$tau, s$ramp)
db <- rbind(data.frame(t = t, y = base, what = "diffusive accumulation (C0 + s t)"),
            data.frame(t = t, y = y0, what = "+ bubbles (noise-free)"))
pb <- ggplot(db, aes(t, y, linetype = what)) +
  geom_vline(xintercept = tb, colour = "red", linetype = 2, linewidth = 0.3) +
  geom_line(linewidth = 0.5) +
  annotate("text", x = tb + 4, y = c(2060, 2150, 2060, 2150), label = paste0("M", 1:4, " = ", M),
           hjust = 0, size = 2.5, colour = "red") +
  scale_linetype_manual(values = c(1, 2), name = NULL) +
  labs(title = "(b) Noise-free: s = 1 ppb/s, 4 bubbles (ppb), overshoot = 1",
       x = "Etime (s)", y = "CH4 (ppb)") + th +
  theme(legend.position = c(0.25, 0.85), legend.background = element_blank())
## ---- (c) Same incubation with the three noise levels of the design ----------------
dc <- do.call(rbind, lapply(c(2, 10, 50), function(sd) {
  set.seed(3); data.frame(t = t, y = y0 + rnorm(length(t), 0, sd), noise = paste0("noise SD = ", sd, " ppb")) }))
dc$noise <- factor(dc$noise, levels = unique(dc$noise))
pc <- ggplot(dc, aes(t, y)) +
  geom_vline(xintercept = tb, colour = "red", linetype = 2, linewidth = 0.3) +
  geom_line(linewidth = 0.25) + facet_wrap(~ noise, nrow = 1) +
  labs(title = "(c) Same incubation with white measurement noise (the three noise levels of the grid)",
       x = "Etime (s)", y = "CH4 (ppb)") + th
## ---- Assemble and save ----------------------------------------------------------------
out_file <- file.path("results", "figures", "FigS1_synthetic_incubations.png")
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
png(out_file, width = 7, height = 7.2, units = "in", res = 300)
grid.arrange(arrangeGrob(pa, pb, ncol = 2, widths = c(1, 1.35)), pc, nrow = 2, heights = c(1, 1))
invisible(dev.off())
cat("Figure S1 written to", out_file, "\n")
