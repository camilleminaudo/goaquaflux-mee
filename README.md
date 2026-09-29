# goAquaFlux: separating bubbles from diffusion in aquatic chamber flux measurements

Code, data and results accompanying the manuscript *"goAquaFlux: separating bubbles from diffusion in aquatic chamber flux measurements"*, submitted to *Methods in Ecology and Evolution*.

goAquaFlux is part of the R package **goFlux** (Rheault et al., 2024). It partitions the greenhouse-gas fluxes measured with floating chambers into a diffusive and an ebullitive (bubble) component. This repository contains the scripts that produce every figure, table and number of the manuscript and its Supporting Information, the field data they use, and the results they produce.

---

## Repository structure

```
.
├── README.md                     this file
├── goaquaflux-mee.Rproj          RStudio project (sets the working directory)
├── LICENSE
├── R/
│   ├── setup.R                                  shared set-up, sourced by the scripts
│   ├── synthetic_incubations.R                  generator of synthetic incubations
│   ├── combine_aqua_plots.R                     helper: panels with one shared legend
│   ├── MethaneSignalProcessor.R                 R port of MethaneSignalProcessor (MSP)
│   ├── make_Fig3_goAquaFlux_demo.R              Figure 3 and Section 4 worked example
│   ├── make_Fig4_goFlux_vs_goAquaFlux.R         Figure 4 and Section 4 field results
│   ├── compare_goAquaFlux_FluxSeparator_MSP.R   benchmark of Section 5 / Appendix S1
│   ├── make_Table3_method_comparison.R          Table 3 and all Section 5 numbers
│   ├── make_FigS1_synthetic_incubations.R       Figure S1
│   ├── check_msp_port_parity.R                  check of the MSP R port against Python
│   └── compare_bubble_detection.R               supplementary development benchmark
├── data/
│   ├── README.md                                description of the field data
│   └── many_incubations.RData                   field incubations (Cabrera-Brufau et al., 2026)
└── results/                                     outputs of the scripts (see below)
```

## Where each result of the manuscript comes from

| Manuscript item | Script | Output (in `results/`) |
|---|---|---|
| Figure 1 (workflow schematic) | drawn outside R | – |
| Figure 2 (bubble detection illustration) | [to be added] | – |
| Figure 3, Section 4 worked example | `make_Fig3_goAquaFlux_demo.R` | `figures/Fig3_worked_example.*`, `Fig3_worked_example_fluxes.csv` |
| Figure 4, Section 4 field results | `make_Fig4_goFlux_vs_goAquaFlux.R` | `figures/Fig4_goFlux_vs_goAquaFlux.*`, `goFlux_vs_goAquaFlux/Fig4_summary.csv` |
| Section 5 benchmark (raw results) | `compare_goAquaFlux_FluxSeparator_MSP.R` | `method_comparison_synthetic/`, `method_comparison_real_incubations/`, `A_capabilities.csv` |
| Table 3 and numbers of Section 5 | `make_Table3_method_comparison.R` | `Table3_method_comparison.csv`, `Section5_numbers.csv` |
| Figure S1 | `make_FigS1_synthetic_incubations.R` | `figures/FigS1_synthetic_incubations.png` |
| Appendix S1, check of the MSP port | `check_msp_port_parity.R` | `msp_parity/msp_port_parity.csv` |

`compare_bubble_detection.R` is a benchmark of the bubble-detection step alone, used while developing `find.bubbles()`. It is not needed to reproduce the manuscript.

## Requirements

- **R** (≥ 4.3; the R and package versions used are recorded in the `sessionInfo.txt` files written next to the results).
- **goFlux** with goAquaFlux, version [version / commit]:
  ```r
  remotes::install_github("Qepanna/goFlux@[commit]")
  ```
  To use a local copy of the goFlux source instead, set `Sys.setenv(GOFLUX_DIR = "path/to/goFlux")` before running a script (`devtools::load_all()` is then used).
- **FluxSeparator v2.0.0** (Sø et al., 2024), for the benchmark only:
  ```r
  remotes::install_github("JonasStage/FluxSeparator@v2.0.0")
  ```
- **MethaneSignalProcessor** (Cardona et al., 2026) is included as an R port (`R/MethaneSignalProcessor.R`); Python is needed only for `check_msp_port_parity.R`.
- CRAN packages: `dplyr`, `tidyr`, `purrr`, `ggplot2`, `zoo`, `TTR`, `broom`, `lubridate`, `patchwork`, `egg`, `ggpubr`, `gridExtra`, `scales`, `svglite`, plus `remotes` and, optionally, `devtools`:
  ```r
  install.packages(c("dplyr", "tidyr", "purrr", "ggplot2", "zoo", "TTR", "broom",
                     "lubridate", "patchwork", "egg", "ggpubr", "gridExtra",
                     "scales", "svglite", "remotes"))
  ```

## How to reproduce the results

All scripts are run from the repository root: open `goaquaflux-mee.Rproj` in RStudio and source the scripts, or run them from a terminal, e.g. `Rscript R/make_Fig3_goAquaFlux_demo.R`.

| Step | Script | Approximate run time |
|---|---|---|
| 1 | `make_Fig3_goAquaFlux_demo.R` | < 1 min |
| 2 | `make_Fig4_goFlux_vs_goAquaFlux.R` | ~10 min |
| 3 | `compare_goAquaFlux_FluxSeparator_MSP.R` | ~1 h (synthetic) + ~15 min (field) |
| 4 | `make_Table3_method_comparison.R` | < 1 min (needs step 3) |
| 5 | `make_FigS1_synthetic_incubations.R` | < 1 min |
| 6 | `check_msp_port_parity.R` | < 1 min (see the procedure in the script header) |

The long computations are cached: in `make_Fig4_goFlux_vs_goAquaFlux.R` and `compare_goAquaFlux_FluxSeparator_MSP.R`, setting `recompute = FALSE` reloads the saved results in `results/` and only redoes the tables and figures.

**Reproducibility notes**

- The synthetic incubations are fully determined by the settings in `CFG$syn` of `compare_goAquaFlux_FluxSeparator_MSP.R`: incubation *j* is generated after `set.seed(20260928 + j)`.
- The detection window of goAquaFlux (15 observations, `BUBBLE_WINDOW` in `R/setup.R`) is passed explicitly in every call, so results do not depend on the default of the installed goFlux version.
- All methods are applied to the same pre-processed series, and all fluxes are converted with the same goFlux flux term. CH₄ fluxes are in nmol m⁻² s⁻¹, CO₂ fluxes in µmol m⁻² s⁻¹.
- Error columns in the result tables (e.g. `med_rel_err`) are unitless fractions: 0.011 means +1.1%.

## Results

`results/` contains the outputs of the scripts, as used in the manuscript:

- `figures/`: Figures 3, 4 and S1;
- `goFlux_vs_goAquaFlux/`: fluxes of the field incubations with and without separation, and the Section 4 summary;
- `method_comparison_synthetic/`: benchmark on synthetic incubations (truth, estimates of every method, summary tables, figures);
- `method_comparison_real_incubations/`: benchmark on the field incubations;
- `Table3_method_comparison.csv`, `Section5_numbers.csv`: Table 3 and the numbers quoted in Section 5;
- `A_capabilities.csv`: capabilities of goAquaFlux, FluxSeparator and MSP;
- `msp_parity/`: check of the MSP R port against the Python implementation.

Appendix S1 (Supporting Information) describes the benchmark design and lists every output file.

## Data

`data/many_incubations.RData` contains the floating-chamber incubations of Cabrera-Brufau et al. (2026), in goFlux format. See `data/README.md`.

## Licence

The code is released under the MIT licence (see `LICENSE`). `R/MethaneSignalProcessor.R` is a port of MethaneSignalProcessor (Cardona et al., 2026) and remains subject to the licence of the original software. The field data are from Cabrera-Brufau et al. (2026); please cite that publication when reusing them.

## References

- Cabrera-Brufau et al. (2026). Assessing the effects of restoration and conservation on gaseous carbon fluxes and climate mitigation potential across six European coastal wetlands. *Ecological Engineering*, 232, 108080. https://doi.org/10.1016/j.ecoleng.2026.108080
- Cardona, A., Butturini, A., & Fonollosa, J. (2026). MethaneSignalProcessor (MSP): Automated discrimination of diffusive and ebullitive methane fluxes at the water–air interface from time-series data. *Ecological Informatics*, 95, 103781. https://doi.org/10.1016/j.ecoinf.2026.103781
- Rheault, K., Christiansen, J. R., & Larsen, K. S. (2024). goFlux: A user-friendly way to calculate GHG fluxes yourself, regardless of user experience. *Journal of Open Source Software*, 9(96), 6393. https://doi.org/10.21105/joss.06393
- Sø, J. S., Sand-Jensen, K., & Kragh, T. (2024). Self-made equipment for automatic methane diffusion and ebullition measurements from aquatic environments. *Journal of Geophysical Research: Biogeosciences*. https://doi.org/10.1029/2024JG008035
