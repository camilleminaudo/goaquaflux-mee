# Field data

`many_incubations.RData` contains the floating-chamber incubations used in Sections 4 and 5 of the manuscript. They come from:

> Cabrera-Brufau et al. (2026). Assessing the effects of restoration and conservation on gaseous carbon fluxes and climate mitigation potential across six European coastal wetlands. *Ecological Engineering*, 232, 108080. https://doi.org/10.1016/j.ecoleng.2026.108080

Please cite that publication when reusing the data.

## Content

The file contains two data frames in goFlux format.

**`mydata_all`**: gas concentration time series (1 Hz), one row per measurement, with the columns used here:

| Column | Description | Unit |
|---|---|---|
| `UniqueID` | incubation identifier | – |
| `POSIX.time` | time stamp | POSIXct |
| `CO2dry_ppm` | CO₂ dry mole fraction | ppm |
| `CH4dry_ppb` | CH₄ dry mole fraction | ppb |
| `H2O_ppm` | water vapour mole fraction | ppm |
| `CO2_prec`, `CH4_prec`, `H2O_prec` | analyser precision | ppm, ppb, ppm |

**`myauxfile`**: one row per incubation (auxiliary information), with the columns used here:

| Column | Description | Unit |
|---|---|---|
| `UniqueID` | incubation identifier | – |
| `start.time` | start of the incubation | POSIXct |
| `duration` | length of the incubation | s |
| `Vtot` | total chamber volume | L |
| `Area` | chamber basal area | cm² |
| `Tcham` | chamber temperature | °C |
| `Pcham` | chamber pressure | kPa |
| `gas_analiser` | gas analyser (LI-COR, Los Gatos, Picarro) | – |

Other columns describe the sampling design (site, habitat, chamber type, light conditions); see Cabrera-Brufau et al. (2026).

The scripts copy `duration` to `obs.length` and use `goFlux::autoID()` to extract each incubation from `mydata_all`.
