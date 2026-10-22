# SimLT: Simulation of Survival Times from a Life Table (R package)

[![R package](https://img.shields.io/badge/language-R-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

`SimLT` is an R package for simulating survival times from a population life
table. It is primarily intended for use in **relative (net) survival** simulation
studies, where background mortality must be incorporated by drawing event times
consistent with nationally or regionally stratified life tables.

Life tables are typically stratified by characteristics such as age, calendar
year, sex, and deprivation level. `SimLT` constructs a **piecewise-constant
hazard function** from the life table rates for each combination of
characteristics, and uses this to simulate individual times to event. An England
life table (`England_LT_2010_2015_dep_gor.txt`), stratified by age, year,
deprivation, and Government Office Region, is bundled directly in the package as
an illustrative dataset.

## Installation

```r
# install.packages("pak")
pak::pak("FJRubio67/SimLT")
library(SimLT)
```

## Main functions

| Function | Description |
|---|---|
| `sim_pophaz` | Simulate background survival times from a life table |
| `simDesMatrix` | Construct a design matrix of life table characteristics for simulation |
| `hpwexp` | Piecewise-exponential hazard function |
| `chpwexp` | Cumulative piecewise-exponential hazard function |
| `rsim.pwexp` | Random sampling from a piecewise-exponential distribution |

For full documentation: `?sim_pophaz`, `?simDesMatrix`, `?hpwexp`, `?chpwexp`,
`?rsim.pwexp`

A PDF manual is also available in the repository: [`simLT.pdf`](simLT.pdf).

## Tutorial

- [Simulating times to event from a life table](https://fjrubio-simlt.share.connect.posit.cloud) —
  illustrative tutorial showing how to construct the piecewise-constant hazard and simulate from a stratified life table.
- [SurvLT: Calculating survival probabilities from a Life Table](https://fjrubio-survlt.share.connect.posit.cloud/) - a tutorial for calculating survival probabilities (at specific time points) for a population of interest, based on the corresponding life tables.

## Related resources

- [HazReg](https://github.com/FJRubio67/HazReg) — parametric hazard-based
  regression models for overall and relative survival (R package)
- [LBANS](https://github.com/FJRubio67/LBANS) — link-based additive net survival
  models; uses life-table background rates of the same form
- [IFNS](https://github.com/FJRubio67/IFNS) — individual frailty net survival
  models (R package)

## License

This package is licensed under the [MIT License](LICENSE).
