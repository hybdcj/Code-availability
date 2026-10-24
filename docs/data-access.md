# Data Sources and Local Layout

The study uses CHARLS, HRS (including RAND harmonized data), ELSA, and SHARE/easySHARE. Individual-level data must be obtained from the relevant providers and used under their conditions. This repository does not distribute raw or derived participant data, model objects containing participant rows, or temporary data caches.

| Cohort | Analysis entry point | Input layout |
|---|---|---|
| CHARLS | `analysis/core/01_build_2020.R`, `03_waves_build.R` | Wave-specific modules under `<project_root>/data/raw/`; see scripts for exact names and capitalization |
| HRS / RAND HRS | `analysis/core/23_external_cohort_harmonize.R` | `<external_staging_root>/hrs/H_HRS_d.dta` and `randhrs1992_2022v1.dta` |
| ELSA | Same as above | `<external_staging_root>/elsa/gh_elsa_h.dta` |
| SHARE / easySHARE | Same as above | `<external_staging_root>/share/easySHARE_rel9-0-0.dta` |

These filenames identify the inputs expected by the study scripts. They do not imply unrestricted download access or interchangeability with other releases. Consult the provider documentation for the actual versions, variable dictionaries, and access conditions. Derived external-cohort data for the main analysis are written to `<project_root>/external_data/`.

`docs/input-references.csv` is a static index of input expressions, not participant data. See `reproduction.md` for third-round temporary aliases and manual preparation steps.

`.gitignore` excludes common data, model, and figure formats, but ignore rules do not replace inspection before committing. Any format, including CSV, can contain participant information. The repository CSV files contain source manifests, input-expression indexes, compatibility mappings, or dependency version records only.
