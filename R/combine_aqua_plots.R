## Combine several flux.plot.aqua() figures with a single, complete legend.
##
## patchwork's guides = "collect" only merges legends that are drawn
## identically. Figures of different gases or incubations contain different
## layers and categories (e.g. bubble fits on CH4 only, "outside diffusive
## window" on CO2 only), so their legends differ and are all kept. Here the
## scales of all figures are first given the union of their categories, then
## one legend is drawn once below the panels.
combine_aqua_plots <- function(plots, ncol = length(plots), merge_fit_labels = FALSE,
                               legend_height = 0.12) {
  built <- lapply(plots, ggplot2::ggplot_build)

  ## Union of the categories, labels and values used by one scale
  collect <- function(aes) {
    lims <- character(0); labs <- character(0); vals <- NULL
    for (b in built) {
      sc <- b$plot$scales$get_scales(aes)
      if (is.null(sc)) next
      l <- sc$get_limits(); lab <- sc$get_labels(l); pal <- sc$palette(length(l))
      for (i in seq_along(l)) {
        j <- match(l[i], lims)
        if (is.na(j)) {
          lims <- c(lims, l[i]); labs <- c(labs, lab[i])
          vals <- c(vals, stats::setNames(pal[[l[i]]], l[i]))
        } else if (nchar(lab[i]) > nchar(labs[j])) {
          labs[j] <- lab[i]   # keep the more explicit label, e.g. "ebullition event (CH4)"
        }
      }
    }
    list(limits = lims, labels = stats::setNames(labs, lims), values = vals)
  }

  obs_order <- c("retained", "outside diffusive window", "discarded")
  shp <- collect("shape"); alp <- collect("alpha")
  obs_lims <- obs_order[obs_order %in% union(shp$limits, alp$limits)]
  key_col  <- c("retained" = "grey15", "outside diffusive window" = "grey75",
                "discarded" = "grey15")[obs_lims]

  col <- collect("colour"); fil <- collect("fill")
  col_breaks <- col$limits; col_labels <- col$labels
  if (merge_fit_labels) {
    ## one key for the diffusive fit, whichever model was selected
    fit_lv <- grep("^diffusive fit \\(", col_breaks, value = TRUE)
    if (length(fit_lv) > 1) {
      col_breaks <- setdiff(col_breaks, fit_lv[-1])
      col_labels[fit_lv[1]] <- "diffusive fit (best model)"
    }
  }

  common <- list(
    ggplot2::scale_shape_manual(NULL, values = shp$values[obs_lims],
                                limits = obs_lims, drop = FALSE),
    ggplot2::scale_alpha_manual(NULL, values = alp$values[obs_lims],
                                limits = obs_lims, drop = FALSE),
    ggplot2::scale_colour_manual(NULL, values = col$values, limits = col$limits,
                                 breaks = col_breaks, labels = col_labels[col_breaks],
                                 drop = FALSE),
    ggplot2::scale_fill_manual(NULL, values = fil$values, limits = fil$limits,
                               labels = fil$labels, drop = FALSE),
    ggplot2::guides(
      shape  = ggplot2::guide_legend(order = 1, override.aes = list(
        colour = unname(key_col), size = 2)),
      alpha  = ggplot2::guide_legend(order = 1),
      colour = ggplot2::guide_legend(order = 2),
      fill   = ggplot2::guide_legend(order = 3, override.aes = list(alpha = 0.35))))
  plots <- lapply(plots, function(p) suppressMessages(p + common))

  ## The legend, taken from the first figure: all its scales now hold every
  ## category, so it is complete.
  g   <- ggplot2::ggplotGrob(plots[[1]] + ggplot2::theme(legend.position = "bottom"))
  idx <- grep("^guide-box", g$layout$name)
  idx <- idx[!vapply(g$grobs[idx], inherits, logical(1), "zeroGrob")][1]
  legend <- g$grobs[[idx]]

  panels <- lapply(plots, function(p) p + ggplot2::theme(legend.position = "none"))
  patchwork::wrap_plots(panels, ncol = ncol) / patchwork::wrap_elements(full = legend) +
    patchwork::plot_layout(heights = c(1, legend_height))
}
