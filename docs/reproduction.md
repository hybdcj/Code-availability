# Reproduction Setup and Execution

## 1. Environment

The checked environment uses Windows and R 4.6.1. `environment/r-packages.csv` records an installed-package snapshot, not a `renv.lock` covering all transitive dependencies. `install_dependencies.R` installs missing direct dependencies without forcing downgrades. Build a separate isolated environment from the snapshot if exact versions are required.

Python helpers require Python 3.10 or later and use the standard library. The measurement-summary script additionally requires openpyxl; see `requirements.txt`. UTF-8 support is required. Original R scripts retain Windows Chinese locale settings for compatibility with historical data. TIFF rendering requires ragg; some historical figures use SimSun. Adjust locales/fonts on other operating systems and verify the results before assuming equivalent cross-platform behavior.

## 2. Configure paths

Copy `config.example.json` to `config.local.json` and supply absolute paths on your machine. Use a new workspace for `project_root`, intermediate data storage for `scratch_root`, extracted external raw files for `external_staging_root`, and your actual R package library for `r_library`.

```powershell
python tools/prepare_workspace.py --config config.local.json
```

The tool copies only text scripts and variable-mapping documentation, replaces known local path prefixes, and writes a preparation manifest. It does not copy, download, or upload data or execute analyses. It stops if a destination source file already has different contents. `source-manifest.csv` verifies the distributed source bytes; prepared scripts have different hashes after path replacement.

## 3. Data and derived inputs

Prepare authorized data as described in `data-access.md`. `input-references.csv` lists read expressions and line numbers; examine script context for dynamically constructed paths, dynamic reads, and report caches. Scripts use fixed dated filenames. Do not rename inputs casually or reuse aliases for RDS files from different periods or definitions.

The main development route is 01 → 03 → identifier verification/33 → 04 → 19 → 20/21. Run other descriptive and longitudinal modules according to their inputs. Script 33 reads `waves_cognition_orig.rds` from scratch storage and produces `waves_cognition_fixed.rds`. After checking its log and identifier consistency, place the corrected result at the five-wave input location before running dependent modules. Repository tools do not perform this replacement automatically.

The external route is 23 → 24 → 25/26, followed by freezing and sensitivity modules using existing development results. This is module navigation, not a substitute for checking each script's inputs. Auditing frozen outputs and refitting from raw data are different operations. Keep version identifiers for old and new results to avoid overwriting manuscript results.

Third-round scripts also use short filenames under `scratch_root/tmp`:

| Alias | Source |
|---|---|
| `waves_cognition.rds` | Verified `output/02_five_wave_data_and_quality_control/waves_cognition_2026-09-13.rds` |
| `cohort_change.rds` | `output/04_longitudinal_cohort/cohort_change_2026-09-13.rds` |
| `a07_memory_v3.rds` | `output/10_screening_index/screening_index_memory_v3_2026-09-26.rds` |
| `ext_hrs.rds` / `ext_elsa.rds` / `ext_share.rds` | Corresponding external-cohort main-analysis RDS, produced by 23 |
| `ext_*_sens.rds` | Corresponding external-cohort sensitivity risk set, produced by 23 |
| `F16_base_*.rds` | Preparation outputs from F16_build |
| `F22_d.rds` | Preparation output from F22a; its historical comparison is not valid original-model validation |
| `F26_bdsi_audit.rds` | Audit output from F26 |

Historical temporary aliases also include `w1_hs.dta`, `w2_hs.dta`, `w3_hs.dta`, `w1_biomarker.dta`, `w3_biomarker.dta`, `w4_cognition.dta`, `w4_hsf.dta`, and `w5_hsf.dta`. Their complete automated preparation chain is not included. Identify the correct modules using variable inventories and wave documentation before staging them locally; similar filenames are insufficient. F29/F30 and related update scripts require existing historical outputs and should not be run as initial reconstruction steps.

## 4. Run and check modules in stages

Example single-module commands after path configuration; these do not imply that the required inputs are ready:

```powershell
Rscript "E:/memory-reproduction/project/code/19_screening_index_memory_v3.R"
Rscript "E:/memory-reproduction/project/16.manuscript_revision_round_three/_audit/scripts/F24_recompute_screening.R" "E:/memory-reproduction/project/16.manuscript_revision_round_three"
```

Before execution from the project root, verify inputs, missing-data rules, output locations, and the intended analysis stage. Check sample/event counts, differences between model and score samples, thresholds, and ties. Scripts write and update outputs. Run formal analyses in a separate directory and compare them with frozen results. Do not commit participant-level outputs or model objects to GitHub.

## 5. F39 subgroup and missingness supplements

The source code for the S37 subgroup statistics and S38 missingness statistics is included in `analysis/revisions/F39_TRIPOD.R`. Run it after the verified five-wave data, development model, and external S0 data are available, then run `F39_check_calibration.R`. See [the F39 reproduction guide](f39-reproduction.md) for exact paths, execution order, output files, and historical consistency checks. These scripts generate statistical tables and reports; they do not insert formatted tables into manuscript Word files.

## 6. Validation boundaries

Release checks cover source integrity, R/Python syntax, workspace preparation, and repository contents. They do not rerun the original statistics, modify the manuscript, or demonstrate clean-environment end-to-end reproduction. See `validation.md` and `known-limitations.md`.
