# Multivariate censored extremes

R code accompanying *Statistics of multivariate extremes under random censoring* by Martin Bladt.

The code implements tail inference for a multivariate distribution under componentwise random right-censoring. Fixing a direction reduces the problem exactly to a univariate right-censored one, so an ordinary Kaplan–Meier product limit estimates the joint tail probability in any dimension, with no multivariate survival function and no smoothing. It includes the two marginal standardizations, Greenwood confidence intervals, the stability rule for threshold selection, competing estimators, simulation studies, and a flood-insurance application.

## Files

- `functions.R`: estimation, inference, simulation, and plotting functions.
- `simulations.R`: Monte Carlo studies of the paper and the Supplementary Material.
- `application.R`: application to Hurricane Ian flood claims from the OpenFEMA NFIP file.

The full simulation uses 500 replications in each scenario and may take substantial time; set `MC_REPS` lower for a quick pass. Figures and numbers are written to `figures/` and `results/`.

`application.R` downloads the claims on first use and caches them. It reads the `v2/FimaNfipClaims` endpoint, which FEMA will update in favour of `v3/NfipClaims`.
