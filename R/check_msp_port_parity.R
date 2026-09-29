###############################################################################
## check_msp_port_parity.R
##
## Checks that the R port of MethaneSignalProcessor (R/MethaneSignalProcessor.R)
## reproduces the outputs of the original Python tool (Cardona et al., 2026;
## https://github.com/ACBonet/MethaneSignalProcessor) on the same input files.
## This validates the MSP results used in compare_goAquaFlux_FluxSeparator_MSP.R
## (Appendix S1, Section S3.6).
##
## Procedure
##   1. Put the same MSP input files (tab-separated .txt, e.g. the example data
##      distributed with MSP) in
##        results/msp_parity/MSP_python/Raw data/
##        results/msp_parity/MSP_R/Raw data/
##   2. Run the Python tool in results/msp_parity/MSP_python/ (python main.py),
##      which writes "Processed data/".
##   3. Run the R port from the repository root:
##        Rscript R/MethaneSignalProcessor.R --dir results/msp_parity/MSP_R
##   4. Run this script from the repository root:
##        Rscript R/check_msp_port_parity.R
##
## For every input file, the script compares
##   * the text summary (<file>_results.txt): identical or not, with the first
##     differing lines printed when they differ;
##   * the processed CH4 series (column "CH4_final (ppm)" of
##     <file>_processed.csv): maximum absolute difference, against `tol`.
##
## Output  results/msp_parity/msp_port_parity.csv
###############################################################################

## ---- Configuration -------------------------------------------------------------
py_dir <- file.path("results", "msp_parity", "MSP_python", "Processed data")  # Python outputs
r_dir  <- file.path("results", "msp_parity", "MSP_R", "Processed data")       # R outputs
tol    <- 1e-9                            # accepted difference on CH4_final (ppm)
out_csv <- file.path("results", "msp_parity", "msp_port_parity.csv")

## ---- Files processed by the R port ---------------------------------------------
files <- sub("_results\\.txt$", "",
             list.files(file.path(r_dir, "results"), "_results\\.txt$"))
if (!length(files)) stop("No *_results.txt found in ", file.path(r_dir, "results"))

## ---- Comparison, file by file ----------------------------------------------------
cmp <- do.call(rbind, lapply(files, function(f) {
  rp <- file.path(py_dir, "results", paste0(f, "_results.txt"))
  cp <- file.path(py_dir, "data", paste0(f, "_processed.csv"))
  if (!file.exists(rp)) {
    return(data.frame(file = f, results_txt = "missing in Python",
                      max_abs_diff_CH4_final = NA))
  }

  ## Text summaries (readLines also handles Windows line endings)
  txt_py <- readLines(rp, encoding = "UTF-8", warn = FALSE)
  txt_r  <- readLines(file.path(r_dir, "results", paste0(f, "_results.txt")),
                      encoding = "UTF-8", warn = FALSE)
  same <- identical(txt_py, txt_r)
  if (!same) {
    n <- seq_len(min(length(txt_py), length(txt_r)))
    k <- which(txt_py[n] != txt_r[n])
    message("\n", f, ": first differing lines",
            "\n  PY: ", paste(head(txt_py[k], 3), collapse = "\n  PY: "),
            "\n  R : ", paste(head(txt_r[k], 3), collapse = "\n  R : "))
  }

  ## Processed CH4 series
  dpy <- read.csv(cp, check.names = FALSE)
  dr  <- read.csv(file.path(r_dir, "data", paste0(f, "_processed.csv")),
                  check.names = FALSE)
  d <- if (nrow(dpy) == nrow(dr)) {
    max(abs(dpy[["CH4_final (ppm)"]] - dr[["CH4_final (ppm)"]]), na.rm = TRUE)
  } else NA

  data.frame(file = f, results_txt = if (same) "identical" else "DIFFERENT",
             max_abs_diff_CH4_final = d)
}))

cmp$CH4_final_ok <- !is.na(cmp$max_abs_diff_CH4_final) & cmp$max_abs_diff_CH4_final < tol
print(cmp, row.names = FALSE)
write.csv(cmp, out_csv, row.names = FALSE)
cat("\n", sum(cmp$results_txt == "identical"), "of", nrow(cmp),
    "summaries identical;", sum(cmp$CH4_final_ok), "of", nrow(cmp),
    "processed series within", tol, "ppm.\nWritten to", out_csv, "\n")
