# Software Validation

## Historical core and F39 release

The following checks cover source integrity and preparation utilities for the historical core and F39 release:

- All 69 R files pass syntax parsing without executing analysis expressions.
- All three Python source files pass AST parsing.
- SHA-256 checks verify all 70 research sources against `source-manifest.csv`.
- The input-reference index matches the recorded source lines, and local Markdown links resolve.
- Workspace preparation creates all 70 source files at the expected locations and replaces the configured absolute path prefixes.
- Repeating preparation with the same configuration is idempotent. The tool refuses to overwrite destination files with different contents.
- Synthetic checks cover screening confusion-matrix calculations, script extraction markers, and interval replacement.
- The F39 addition was checked with isolated subgroup and variable-missingness functions on synthetic data, including absent-variable reporting. Preparation was repeated successfully, migrated F39 paths were checked, and original local source hashes were verified unchanged. The recovered historical analysis was not rerun on participant data.

The checks use Windows and R 4.6.1. The R syntax checker sets `Chinese_China.utf8` explicitly because the default startup locale is unavailable in the checked environment. Python helpers require Python 3.10 or later.

Those historical checks did not download participant data, refit models, run bootstrap or imputation analyses, or regenerate the complete set of manuscript figures. They do not establish end-to-end reproducibility or independently validate the statistical results. See [known-limitations.md](known-limitations.md) for analysis and data dependencies.

## Wave 6 release: 2026-11-16

The nine published Wave 6 analysis and verification scripts were executed in order in an isolated prepared workspace on Windows, Python 3.14, and R 4.6.1. The workspace contained locally authorized raw Wave 6/2018 files and the existing derived historical cohort and frozen model artifacts listed in the [Wave 6 guide](wave6-reproduction.md).

- All 82 source/protocol manifest entries passed SHA-256 verification and workspace preparation.
- All Python sources passed AST parsing; every new R source executed successfully.
- All nine analysis/verification stages exited successfully, including QC, event audit, frozen temporal validation, 10,000 intercept-update bootstrap draws, and independent R verification of the full fit and first 50 draws.
- Forty generated CSV files, including local participant-membership/prediction files and aggregate performance tables, matched the original analysis with relative tolerance 1e-10 and absolute tolerance 1e-12. Identifier columns were compared as strings. Those outputs remain private and are not included in the release.
- Local Markdown file links resolved, and the release delta contained only analysis sources, protocol summaries, dependency specifications, source manifests, and documentation.

R emitted startup locale warnings in this Windows environment; all stages completed and numerical checks passed. Best-effort locale selection is retained for portability.

This validates execution of the new extension from its stated prerequisites. It does not reconstruct the historical derived cohorts or refit the original development model from raw data, establish a fully locked environment, or regenerate all manuscript artifacts. Independent R calculations are computational cross-checks, not independent-sample validation of the updated model.
