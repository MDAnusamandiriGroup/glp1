# Expected numerical smoke output

The earlier uncalibrated smoke values are superseded. An independent numerical
implementation of the current equations gives the following deterministic
base-case checksum. Small solver-level differences are expected.

| Quantity | Expected value |
|---|---:|
| Calibrated placebo first-MACE hazard | 0.025291895/year |
| Calibrated latent treated-state HR | 0.789379725 |
| Annualized SELECT discontinuation probability | 0.089399611 |
| Placebo risk at 39.8 months | 0.079650040 |
| Semaglutide risk at 39.8 months | 0.064637056 |
| No-treatment first MACE at 5 years | 11,700.035/100,000 |
| Full-persistence first MACE at 5 years | 9,356.595/100,000 |
| SELECT-like first MACE at 5 years | 9,634.483/100,000 |
| Real-world first MACE at 5 years | 10,106.924/100,000 |
| Real-world MACE prevented | 1,593.111/100,000 |
| Real-world benefit retained vs full persistence | 0.679817 |
| Real-world benefit retained vs SELECT-like implementation | 0.771276 |

The current pipeline also writes its own checks to:

- `results/shared/select_calibration_validation.csv`;
- `results/shared/calibrated_parameters.csv`;
- `results/shared/run_metadata.csv`.

Required pass conditions are:

1. modeled placebo and semaglutide first-MACE risks match the SELECT calibration
   targets to numerical tolerance;
2. the calibrated treated-state HR is below the treatment-policy HR of 0.80;
3. total population is conserved to numerical tolerance;
4. increasing persistence never reduces benefit retained;
5. EJPC and DOM outputs differ only through their declared benchmark.

Run `source("tests/testthat.R")` after any model change. The automated checksum
allows a tolerance of 0.05 events per 100,000. Final manuscript values must
come from a run with `N_PSA=10000`.
