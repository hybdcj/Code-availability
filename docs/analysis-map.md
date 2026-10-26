# Analysis Guide and Version Precedence

Filename prefixes are historical project identifiers. They are not a strict execution order. Do not traverse the directory and run every R script as a batch.

| Module | Main scripts, relative to repository root | Dependencies and notes |
|---|---|---|
| 2020 data preparation | [01_build_2020.R](../analysis/core/01_build_2020.R) | Raw CHARLS modules; produces derived 2020 baseline data |
| Five-wave integration | [03_waves_build.R](../analysis/core/03_waves_build.R) | Raw 2011–2018 data and derived 2020 data |
| Participant identifier repair | [33_fix_waves_cognition_idn.R](../analysis/core/33_fix_waves_cognition_idn.R) | After integration; repairs a temporary copy that must be verified before replacement |
| Development source cohort | [04_cohort_change.R](../analysis/core/04_cohort_change.R) | Five-wave data with verified identifiers |
| Measurement and longitudinal analyses | `analysis/core/05*`–`16*` | Risk factors, trajectories, measurement equivalence, incidence, weighting, and attrition; identical numeric prefixes may identify different modules |
| Memory-outcome score development | [19_screening_index_memory_v3.R](../analysis/core/19_screening_index_memory_v3.R) | Main memory-outcome entry point; 18 is a historical script for another cognitive outcome |
| Bootstrap internal validation | [20_internal_validation_bootstrap.R](../analysis/core/20_internal_validation_bootstrap.R), [21_internal_validation_figures.R](../analysis/core/21_internal_validation_figures.R) | Existing score/model objects |
| External data harmonization | [23_external_cohort_harmonize.R](../analysis/core/23_external_cohort_harmonize.R) | Independently authorized HRS/RAND, ELSA, and SHARE data |
| External validation | [24_external_validation.R](../analysis/core/24_external_validation.R) | Analysis data from 23 and frozen model from 19 |
| External figures/tables and heterogeneity | [25_a12_figures_report.R](../analysis/core/25_a12_figures_report.R), [26_a13_cross_cohort_heterogeneity.R](../analysis/core/26_a13_cross_cohort_heterogeneity.R) | External-validation outputs |
| Model freezing and checks | [27_finalize_model_A14.R](../analysis/core/27_finalize_model_A14.R), [28_A14_final_freeze_checklist.R](../analysis/core/28_A14_final_freeze_checklist.R) | Development and validation results |
| Sensitivity analyses | [29_A15_sensitivity_analyses.R](../analysis/core/29_A15_sensitivity_analyses.R) | Distributed version includes corrected row alignment and imputation summaries; F29/F30 update and audit historical products |
| Historical figure/table exports | [31_A16_tables.R](../analysis/core/31_A16_tables.R), [32_A16_figures.R](../analysis/core/32_A16_figures.R) | Apply subsequent corrections; these are not automatically final manuscript figures/tables |
| Third-round period and broader-population preparation | `analysis/revisions/F02*`–`F22*` | Preparation expressions used by later scripts; includes superseded result branches |
| Screening-metric correction | [F24_recompute_screening.R](../analysis/revisions/F24_recompute_screening.R) | References F03/F14/F15/F16/F19/F22; takes the third-round directory as an argument; later BDSI/Tu restrictions still apply |
| BDSI correction | [F26_bdsi_audit.R](../analysis/revisions/F26_bdsi_audit.R) | F22 derived data and relevant raw variables; supersedes historical BDSI comparisons |
| Validation excluding development participants | [F27_nonoverlap_validation.R](../analysis/revisions/F27_nonoverlap_validation.R) | F26 results, development objects, and cohort data |
| Outcome-definition registry | [F28_outcome_registry.R](../analysis/revisions/F28_outcome_registry.R) | Frozen development and external data |
| Internal-performance row-alignment audit | [F29_fix_internal_performance_order.R](../analysis/revisions/F29_fix_internal_performance_order.R) | Existing A15 outputs; not a prerequisite for reconstruction from scratch |
| Imputation-interval audit/update | [F30_audit_imputation_intervals.R](../analysis/revisions/F30_audit_imputation_intervals.R), [F30_update_imputation_intervals.R](../analysis/revisions/F30_update_imputation_intervals.R) | Audit per-imputation performance first, then update outputs; 2.5%–97.5% limits are empirical quantile intervals |
| Period and measurement audit | [F31_audit_period_and_measurement.R](../analysis/revisions/F31_audit_period_and_measurement.R) | F04 preparation expressions and development model |
| Broader-population denominator correction | [F32_correct_population_eligibility_and_denominators.R](../analysis/revisions/F32_correct_population_eligibility_and_denominators.R) | F15 preparation expressions; integer-score eligibility no longer additionally requires complete ADL |
| Revised figures/tables | [F34_update_main_figures.R](../analysis/revisions/F34_update_main_figures.R), [F34_build_supplementary_tables.R](../analysis/revisions/F34_build_supplementary_tables.R), [F36_longitudinal_figures.R](../analysis/revisions/F36_longitudinal_figures.R) | TIFF figures and native-format tables; existing results and some layout caches are required |
| TRIPOD subgroup and missingness supplements (S37/S38) | [F39_TRIPOD.R](../analysis/revisions/F39_TRIPOD.R), then [F39_check_calibration.R](../analysis/revisions/F39_check_calibration.R) | Verified five-wave data, frozen development model, and script 23 external S0 risk sets; see [F39 reproduction](f39-reproduction.md) for inputs, outputs, and population definitions |

Workspace preparation maps `analysis/core/*` to `<project_root>/code/*` and `analysis/revisions/*` to `<project_root>/16.manuscript_revision_round_three/_audit/scripts/*`, keeping inter-script references consistent.
