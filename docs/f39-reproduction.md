# F39 subgroup and missing-data analyses

The F39 source scripts generate the subgroup statistics underlying Supplementary Table S37 and the variable-level missingness statistics underlying S38. They read existing derived data and apply the frozen model; they do not refit the prediction model. Calibration regressions are evaluation steps. The scripts generate CSV, optional XLSX, and a Markdown report, not the formatted manuscript Word tables.

## Source and adaptation

The scripts were recovered from the local F39 analysis dated 2026-10-17. Only local junction paths, introductory path comments, and text encoding/line endings were adapted for this repository. Statistical expressions and historical output names are preserved. The duplicate Chinese-named main script has identical source bytes and is not distributed a second time.

Original local SHA-256 checksums, before path adaptation:

| Source | SHA-256 |
|---|---|
| `F39_TRIPOD.R` | `bb40c3ef854077f045029ff12aa483199528d3deebc1fa69ed5fc54cfeed8d94` |
| `F39_check_calibration.R` | `851c8d095cbb32c38c4c6b157fffb0540c9f15f2bf2fc35c4cd4a3ab2cc18241` |

The adapted files are covered by `source-manifest.csv`. Neither participant data nor model objects, local outputs, or run logs are included in this addition.

## Prepare inputs

Run `tools/prepare_workspace.py` as described in [reproduction.md](reproduction.md). It installs both scripts under `<project_root>/16.manuscript_revision_round_three/_audit/scripts/` and replaces the canonical project and R-library prefixes. The old `D:/xu_*` junctions are not required.

Required inputs under the configured project root:

| Input | Producer or prerequisite |
|---|---|
| `output/02_five_wave_data_and_quality_control/waves_cognition_2026-09-13.rds` | Script 03; verified identifiers and required `female`, `rural`, cognition, predictor and survey fields |
| `output/10_screening_index/screening_index_memory_v3_2026-09-26.rds` | Script 19; frozen `model_u` and development data |
| `external_data/HRS/analysis/HRS_analysis_2026-09-26.rds` | Script 23, S0 risk set |
| `external_data/ELSA/analysis/ELSA_analysis_2026-09-26.rds` | Script 23, S0 risk set |
| `external_data/SHARE/analysis/SHARE_analysis_2026-09-26.rds` | Script 23, S0 risk set |

The main script checks frozen coefficients against the saved model, development participant identifiers, scores, and historical sample counts (4,441 participants, 658 events; broad complete-case sample 8,296). These checks intentionally target the manuscript version, not arbitrary new datasets.

## Execute in order

Example after preparing the workspace; replace the example project path with your configuration:

```powershell
Rscript "E:/memory-reproduction/project/16.manuscript_revision_round_three/_audit/scripts/F39_TRIPOD.R"
Rscript "E:/memory-reproduction/project/16.manuscript_revision_round_three/_audit/scripts/F39_check_calibration.R"
```

Outputs are written to `<project_root>/16.manuscript_revision_round_three/F39_TRIPOD_supplement_2026-10-17/`. Re-execution overwrites files in that output directory. Run in a dedicated reproduction workspace.

| Manuscript target | Generated files |
|---|---|
| S37, development subgroups | `F39_23a_CHARLS_dev_subgroups.csv` |
| S37, external S0 subgroups | `F39_23a_external_S0_subgroups.csv` |
| Combined subgroup results | `F39_23a_all_subgroups.csv` |
| S38, CHARLS stages | `F39_20b_CHARLS_missing_by_var.csv` |
| S38, external risk sets | `F39_20b_external_missing_by_var.csv` |
| Calibration cross-check | `F39_check_calibration_O_vs_E.csv` |
| Supporting summaries | `F39_TRIPOD_report_zh.md`, `F39_TRIPOD_tables.xlsx` (when openxlsx is installed), `F39_run_log.txt` |

The calibration cross-check runs after the main script and appends `cal_note` to the development and combined subgroup CSVs. It explains the near-zero apparent education-stratum offset intercepts rather than changing the fitted model.

## Interpretation and coverage

CHARLS includes sex, education, age bands, age below/above 65, and urban/rural groups. The saved external analysis inputs lack sex and rural status, so external subgroup results cover available education and age groups only. S37 uses the continuous-model complete-case external rows (`*_S0_modelComplete`); the CSVs additionally contain score-full-risk-set rows (`*_S0_scoreFull`) and integer-score AUCs. Calibration pertains to continuous-model probabilities, not to the integer score. In score-full rows, probability-based metrics exclude missing predictions whereas `n` describes the full score sample; do not use that full `n` as the denominator for continuous-model metrics.

The code includes more rows and metrics than the manuscript table. Formatting and selecting rows for S37/S38 must preserve the intended populations. Source recovery and software checks do not establish a fresh end-to-end statistical rerun; data authorization, derived-input preparation, and historical version compatibility remain prerequisites.
