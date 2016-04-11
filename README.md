# isocyclr

[![Build Status](https://travis-ci.org/sebkopf/isocyclr.svg)](https://travis-ci.org/sebkopf/isocyclr)
[![codecov.io](https://codecov.io/github/sebkopf/isocyclr/coverage.svg?branch=master)](https://codecov.io/github/sebkopf/isocyclr?branch=master)


This R package is intended as an educational tool to facilitate the modeling of isotopic effects through biogeochemical pathways and cycles including non-steady state and steady-state solutions. Core functionality for the construction of reaction networks and the corresponding system of differential equations is fully implemented and can be explored for multiple isotope systems in paralle. However, this is still very much a work in progress and only generic custom and standard reaction types are available at this point as well as rudimentary computation of steady-state solutions. All implemented functionality is automatically tested with pretty good coverage to avoid updates accidentally breaking existing code.

## Installation

```
install.packages("devtools")
devtools::install_github("sebkopf/isocyclr")
```

## Example
