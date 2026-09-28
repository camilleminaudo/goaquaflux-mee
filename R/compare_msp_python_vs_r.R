## Compare MethaneSignalProcessor outputs: Python vs R port
## ---- CONFIG ------------------------------------------------------------------
py_dir <- "MSP_python/Processed data"   # output folder of the Python run
r_dir  <- "MSP_R/Processed data"        # output folder of the R run
tol    <- 1e-9                          # tolerance on CH4_final (ppm)
## ---------------------------------------------------------------------------------

files <- sub("_results\\.txt$", "", list.files(file.path(r_dir, "results"), "_results\\.txt$"))
if (!length(files)) stop("No *_results.txt found in ", file.path(r_dir, "results"))

cmp <- do.call(rbind, lapply(files, function(f) {
  rp <- file.path(py_dir, "results", paste0(f, "_results.txt"))
  cp <- file.path(py_dir, "data", paste0(f, "_processed.csv"))
  if (!file.exists(rp)) return(data.frame(file = f, results_txt = "missing in Python", max_abs_diff_CH4_final = NA))
  txt_py <- readLines(rp, encoding = "UTF-8", warn = FALSE)   # handles CRLF from Windows
  txt_r  <- readLines(file.path(r_dir, "results", paste0(f, "_results.txt")), encoding = "UTF-8", warn = FALSE)
  same <- identical(txt_py, txt_r)
  if (!same) {
    k <- which(txt_py[seq_len(min(length(txt_py), length(txt_r)))] != txt_r[seq_len(min(length(txt_py), length(txt_r)))])
    message("\n", f, ": first differing lines\n  PY: ", paste(head(txt_py[k], 3), collapse = "\n  PY: "),
            "\n  R : ", paste(head(txt_r[k], 3), collapse = "\n  R : "))
  }
  dpy <- read.csv(cp, check.names = FALSE)
  dr  <- read.csv(file.path(r_dir, "data", paste0(f, "_processed.csv")), check.names = FALSE)
  d <- if (nrow(dpy) == nrow(dr)) max(abs(dpy[["CH4_final (ppm)"]] - dr[["CH4_final (ppm)"]]), na.rm = TRUE) else NA
  data.frame(file = f, results_txt = if (same) "identical" else "DIFFERENT", max_abs_diff_CH4_final = d)
}))

cmp$CH4_final_ok <- !is.na(cmp$max_abs_diff_CH4_final) & cmp$max_abs_diff_CH4_final < tol
print(cmp, row.names = FALSE)
