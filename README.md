# How much cardiovascular benefit survives real-world GLP-1 discontinuation?

**A dynamic compartmental model of semaglutide persistence, re-initiation, and cardiometabolic rebound**

*Manuscript in preparation (2026).*

Reproducible R proof-of-concept estimating how discontinuation, post-cessation
rebound, and re-initiation attenuate semaglutide's population-level prevention
of first major adverse cardiovascular events (MACE).

## Primary question

Among adults matching the core SELECT population (age >=45 years, BMI >=27
kg/m2, established cardiovascular disease, and no diabetes), what fraction of
the cardiovascular benefit expected under trial-like or full treatment
persistence survives under real-world discontinuation and re-initiation?

The repository has one model engine and two reporting profiles:

- **EJPC:** real-world implementation versus a SELECT-like implementation;
- **DOM:** real-world persistence versus full treatment persistence.

Both profiles use the same calibrated parameters and Monte Carlo draws. The
different denominators are explicit outputs, not two competing model versions.

The first paper should report four primary outputs:

1. first MACE prevented per 100,000 initiators;
2. fatal-or-undetermined first MACE prevented per 100,000 initiators;
3. fraction of the SELECT-like and full-persistence benefit retained;
4. minimum one-year persistence required to retain 50%, 75%, and 90% of the
   journal-profile benchmark benefit.

This is a semi-Markov approximation implemented as an ODE model. Eight rebound
stages retain time-since-discontinuation information that a single `R` state
would lose.

## Model states

- `U`: eligible and untreated
- `T`: currently treated
- `R1`-`R8`: discontinued, ordered by time since cessation
- `M_nonfatal`: absorbed after a first nonfatal MACE
- `D_fatal_mace`: absorbed after fatal/undetermined first MACE
- `D_ncv`: absorbed after non-cardiovascular death before a first MACE

See [docs/model-specification.md](docs/model-specification.md) for the equations,
estimands, structural assumptions, and parameter provenance.
See [docs/publication-profiles.md](docs/publication-profiles.md) for the EJPC-to-
DOM submission switch and the exact output mapping.
An independent numerical checksum is in
[docs/expected-smoke-output.md](docs/expected-smoke-output.md).

## Released results and reproducibility

This public snapshot contains aggregate tables reported as originating from
the prior 10,000-draw analysis in `results/ejpc/` and the derived NHANES denominators in
`data/derived/`. These aggregate files can be inspected without
rerunning the model. Independent end-to-end verification has not been performed
on this public snapshot. Raw NHANES transport files and the full 10,000-row PSA RDS
are intentionally not tracked.

To reproduce the full EJPC analysis from public source data, open the project
root in RStudio and run:

```r
Sys.setenv(N_PSA = 10000)
source("run_ejpc.R")
```

For a fast installation check, use `N_PSA = 100`. Smoke-test output is not
manuscript output.

An optional supplementary post-processing script is available at
`analysis/07_ejpc_editor_strengthening.R`. It consumes the released
aggregate tables and does **not** rerun the PSA.

To generate the DOM reporting profile using the same model engine and draws:

```r
source("run_dom.R")
```

The pipelines install missing CRAN packages, download public NHANES files, build
the SELECT-like population denominator, calibrate the model to both SELECT
trial arms, run deterministic scenarios and probabilistic sensitivity analysis,
and generate journal-specific figures and tables locally. The public ZIP
includes currently available aggregate EJPC tables and one supplementary EJPC
figure; it does not contain the complete main figure set or precomputed DOM outputs.

## Public data downloaded automatically

The script `data-raw/download_nhanes.R` downloads these five official NHANES
August 2021-August 2023 transport files:

| File | Purpose | Direct data | Documentation |
|---|---|---|---|
| `DEMO_L.XPT` | age, survey design, phlebotomy weight | [XPT](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/DEMO_L.XPT) | [codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/DEMO_L.htm) |
| `BMX_L.XPT` | measured BMI | [XPT](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/BMX_L.XPT) | [codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/BMX_L.htm) |
| `DIQ_L.XPT` | self-reported diabetes | [XPT](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/DIQ_L.XPT) | [codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/DIQ_L.htm) |
| `GHB_L.XPT` | HbA1c, to exclude HbA1c >=6.5% | [XPT](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/GHB_L.XPT) | [codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/GHB_L.htm) |
| `MCQ_L.XPT` | prior MI and stroke | [XPT](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/MCQ_L.XPT) | [codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/MCQ_L.htm) |

No CDC WONDER download is required in version 1. Baseline MACE and competing
death rates are anchored directly to SELECT placebo person-time, which avoids
mixing general-population mortality with a secondary-prevention cohort.

NHANES cannot identify the full SELECT definition of symptomatic peripheral
arterial disease, recent events within 60 days, or every exclusion criterion.
The NHANES estimate is therefore a transparent, conservative operational
denominator based on prior MI or stroke; it is not an exact replication of
trial eligibility.

## Evidence anchors

- SELECT MACE effect and event counts: [NEJM, DOI 10.1056/NEJMoa2307563](https://doi.org/10.1056/NEJMoa2307563)
- SELECT endpoint person-time table: [NCBI Bookshelf Table 2](https://www.ncbi.nlm.nih.gov/books/NBK619926/table/tr8269868382084182_ch01_t02/)
- Real-world discontinuation and re-initiation: [JAMA Network Open, DOI 10.1001/jamanetworkopen.2024.57349](https://doi.org/10.1001/jamanetworkopen.2024.57349)
- Nonlinear post-cessation regain: [eClinicalMedicine, DOI 10.1016/j.eclinm.2026.103796](https://doi.org/10.1016/j.eclinm.2026.103796)
- SELECT long-term weight change: [Nature Medicine, DOI 10.1038/s41591-024-02996-7](https://doi.org/10.1038/s41591-024-02996-7)
- SELECT eligibility details: [Diabetes Care / PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11282385/)

## Repository structure

```text
glp1/
├── R/                         model, simulation, PSA, and plotting functions
├── analysis/                  numbered end-to-end analyses
│   └── 06_ejpc_robustness.R   robustness + figure refresh; reuses saved PSA
├── config/parameters.csv      all fixed and uncertain inputs
├── data-raw/                  download script only
├── data/raw/                  downloaded XPT files; git-ignored
├── data/derived/              released aggregate NHANES estimates
├── docs/                      model specification and data dictionary
├── figures/                   generated figures
│   ├── ejpc/                  SELECT-like benchmark figures
│   └── dom/                   full-persistence benchmark figures (generated locally)
├── results/                   released aggregate tables; large RDS ignored
│   ├── shared/                calibration, shared simulations, PSA draws
│   ├── ejpc/                  EJPC-specific retained-benefit outputs
│   └── dom/                   DOM-specific retained-benefit outputs (generated locally)
├── tests/testthat/            mass-balance and monotonicity tests
├── run_all.R                  generate both publication profiles
├── run_ejpc.R                 generate EJPC outputs only
├── run_dom.R                  generate DOM outputs only
└── glp1.Rproj
```

## Interpretation guardrails

- `HR=0.80` is the randomized SELECT treatment-policy effect and is not used
  directly as a full-adherence effect. The latent treated-state effect is
  calibrated so that a SELECT-like discontinuation trajectory reproduces the
  reported trial arms. This avoids applying real-world discontinuation on top
  of an uncorrected intention-to-treat effect.
- The mapping from remaining weight-loss benefit to remaining cardiovascular
  benefit after cessation is an explicit structural assumption. It is varied
  in sensitivity analysis through `effect_weight_exponent`.
- `D_fatal_mace` is fatal or undetermined **first MACE**, not all-cause deaths.
- The 64.8% discontinuation and 36.3% re-initiation estimates come from mixed
  real-world GLP-1 RA use among patients without diabetes, not a semaglutide-only
  SELECT-eligible cohort. This transportability limitation belongs prominently
  in the paper.

## GitHub setup

Create an empty public GitHub repository with no generated README or licence,
then run from this project root:

```bash
git init -b main
git add .
git commit -m "Initial reproducible GLP-1 persistence model"
git remote add origin https://github.com/MDAnusamandiriGroup/glp1.git
git push -u origin main
```

Before pushing, run `Rscript tests/testthat.R` and inspect `git status`. Raw
NHANES files, full Monte Carlo objects, manuscript drafts, and local R state are
excluded by `.gitignore`. For browser uploads: create the public `glp1` repository, use **Add file → Upload files**,
select all contents of this folder (not the outer folder), commit, and check
the repository URL in a private browser window. No manuscript drafts are included.
