# Wave 6 reproduction

This extension reproduces the 2018-to-2021/2022/2023 temporal analyses and the intercept-only internal validation. It requires authorized CHARLS data and existing frozen development artifacts; it is not a raw-data-to-manuscript reconstruction of every historical analysis.

## Prepare the workspace

Use Python 3.14, R, and the R packages listed under `environment/` (including pROC). The extension was tested with Python 3.14 and R 4.6.1 on Windows; other Python versions must satisfy the pinned packages' own requirements and were not tested. Install the additional Python packages with:

```text
python -m pip install -r environment/requirements-wave6.txt
python tools/verify_sources.py
python tools/prepare_workspace.py --config /absolute/path/to/config-local.json
```

Configure an absolute `project_root` outside this repository using the existing configuration example and [preparation instructions](reproduction.md). The preparation tool copies code and protocol summaries only. Put your authorized inputs in that workspace:

| Required file | Location relative to project_root |
| --- | --- |
| 2018 cognition and sample files | data/raw/wave4_2018/Cognition.dta; data/raw/wave4_2018/Sample_Infor.dta |
| Annual Wave 6 cognition and demographics | data/raw/wave6_2021_2023/{2021,2022,2023}/Cognition_Depression.dta and Demographics.dta |
| Wave 6 sample information | data/raw/wave6_2021_2023/2023/Sample_Infor.dta |
| Derived historical cognitive records | output/02_five_wave_data_and_quality_control/waves_cognition_2026-09-13.rds |
| Derived development cohort | output/04_longitudinal_cohort/cohort_change_2026-09-13.rds |
| Frozen fitted model object | output/12_final_model/A14_final_model_2026-09-26.rds |

The R exporters require one matching artifact per historical filename under output. The frozen model object must contain `$model` and its development data (4,441 participants, 658 events). The event exporter derives full-precision coefficients directly from that model. These restricted inputs and generated participant-level exports are not distributed here.

## Execution order

Change the current directory to the configured project_root. Run the following commands in order, stopping if any fails:

```text
Rscript code/export_wave6_qc_inputs_20261113.R
python code/wave6_qc_20261113.py
python code/verify_wave6_qc_20261113.py
Rscript code/export_wave6_event_inputs_20261114.R
python code/wave6_event_audit_20261114.py
python code/validate_wave6_20261114.py
Rscript code/verify_wave6_validation_20261114.R
python code/validate_wave6_intercept_update_20261116.py
Rscript code/_checks/verify_wave6_intercept_update_20261116.R
```

Outputs appear under `output/wave6_qc_20261113`, `output/wave6_events_20261114`, `output/wave6_validation_20261114`, and `output/wave6_intercept_validation_20261116`. Keep these outside the repository: some outputs contain participant identifiers, covariates, predictions, or bootstrap membership.

## Interpretation and expected checks

The annual nonoverlapping evaluation sets contain 4,041, 1,188, and 672 participants, with 785, 269, and 160 events. The four-year analysis uses 1,188 participants in 212 communities. Expected four-year frozen AUC is approximately 0.685815, with intercept update 0.499309 and calibration slope 1.033425.

The intercept-only analysis resamples 212 communities with replacement 10,000 times (seed 20261116), holds every original predictor coefficient fixed, re-estimates only the intercept, and compares bootstrap training loss with loss on the original evaluation set. Add the mean test-minus-training Brier difference to apparent updated Brier loss.

| Quantity | Expected value |
| --- | ---: |
| Frozen Brier score | 0.166192323801 |
| Updated apparent Brier score | 0.160236888712 |
| Mean loss optimism correction | 0.000297193034 |
| Updated optimism-corrected Brier score | 0.160534081746 |
| Corrected minus frozen Brier score | -0.005658242056 |
| Monte Carlo SE of correction | 0.000070513239 |

The Monte Carlo SE describes finite bootstrap simulation error, not sampling uncertainty. No sampling confidence interval is calculated for the corrected estimate or its difference. This is internal validation of an intercept update in the temporal evaluation sample, not validation of the updated model in an independent sample. The frozen model's ranking is unchanged by an intercept shift; measurement differences across waves and the evaluation sample's limitations remain.

## Provenance

The nine source scripts and three English protocol summaries are checksum-tracked in `source-manifest.csv`. [Source provenance](wave6-source-provenance.csv) records the original script hashes and publication adaptations: portable paths, best-effort R locale selection, unambiguous artifact discovery, and exporting model coefficients from the frozen R object. The English protocol summaries reflect dated local analysis specifications; their publication does not constitute prospective registration. Historical dates in filenames identify analysis versions.

See [software validation](validation.md) for what was actually executed for this release and [known limitations](known-limitations.md) for the broader repository.

