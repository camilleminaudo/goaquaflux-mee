###############################################################################
## setup.R
##
## Shared set-up sourced at the top of every analysis script of this
## repository. It
##   1. checks that the working directory is the repository root,
##   2. defines the paths to the data and results folders,
##   3. loads goFlux (which contains goAquaFlux) and checks its version,
##   4. defines the settings shared by all scripts (bubble detection window).
##
## All scripts must be run from the repository root, e.g.
##   * in RStudio: open goaquaflux-mee.Rproj, then source the script;
##   * in a terminal: Rscript R/<script>.R  (from the repository root).
##
## goFlux is loaded from the installed package (see README.md for the version
## to install). To use a local clone of the goFlux source instead, set the
## environment variable GOFLUX_DIR to its root folder before sourcing, e.g.
##   Sys.setenv(GOFLUX_DIR = "path/to/goFlux")
## in which case devtools::load_all() is used.
###############################################################################

## ---- 1. Repository root -----------------------------------------------------
if (!file.exists(file.path("R", "setup.R")) || !dir.exists("data")) {
  stop("Please run the scripts from the repository root (the folder that ",
       "contains README.md, R/, data/ and results/). In RStudio, open ",
       "goaquaflux-mee.Rproj first.", call. = FALSE)
}

## ---- 2. Paths -----------------------------------------------------------------
PATHS <- list(
  data      = "data",
  results   = "results",
  figures   = file.path("results", "figures"),
  real_data = file.path("data", "many_incubations.RData")
)
dir.create(PATHS$figures, showWarnings = FALSE, recursive = TRUE)

## ---- 3. goFlux / goAquaFlux ---------------------------------------------------
if (!exists("goAquaFlux", mode = "function")) {
  goflux_dir <- Sys.getenv("GOFLUX_DIR", unset = "")
  if (nzchar(goflux_dir)) {
    devtools::load_all(goflux_dir, quiet = TRUE)
  } else {
    suppressPackageStartupMessages(library(goFlux))
  }
}
if (!exists("goAquaFlux", mode = "function")) {
  stop("goAquaFlux() was not found: install the goFlux version given in ",
       "README.md, or set GOFLUX_DIR to a local clone of goFlux.", call. = FALSE)
}

## Internal goFlux functions used by the benchmark scripts. They are not
## exported by the package, so they are not visible after library(goFlux)
## (devtools::load_all() exposes them); they are taken from its namespace.
##   flux.term()    : chamber flux term (compare_goAquaFlux_FluxSeparator_MSP.R)
##   find.bubbles() : bubble detection   (compare_bubble_detection.R)
for (fn in c("flux.term", "find.bubbles")) {
  if (!exists(fn, mode = "function")) assign(fn, utils::getFromNamespace(fn, "goFlux"))
}

## ---- 4. Shared settings ---------------------------------------------------------
## Width (observations on the 1-s grid) of the rolling window used by
## find.bubbles() to detect ebullition. It is passed explicitly in every call
## so that the results do not depend on the default of the installed goFlux
## version. 15 is the default used throughout the manuscript.
BUBBLE_WINDOW <- 15