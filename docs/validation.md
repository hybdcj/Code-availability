# Software Validation

The following checks cover source integrity and preparation utilities:

- All 69 R files pass syntax parsing without executing analysis expressions.
- All three Python source files pass AST parsing.
- SHA-256 checks verify all 70 research sources against `source-manifest.csv`.
- The input-reference index matches the recorded source lines, and local Markdown links resolve.
- Workspace preparation creates all 70 source files at the expected locations and replaces the configured absolute path prefixes.
- Repeating preparation with the same configuration is idempotent. The tool refuses to overwrite destination files with different contents.
- Synthetic checks cover screening confusion-matrix calculations, script extraction markers, and interval replacement.
- The F39 addition was checked with isolated subgroup and variable-missingness functions on synthetic data, including absent-variable reporting. Preparation was repeated successfully, migrated F39 paths were checked, and original local source hashes were verified unchanged. The recovered historical analysis was not rerun on participant data.

The checks use Windows and R 4.6.1. The R syntax checker sets `Chinese_China.utf8` explicitly because the default startup locale is unavailable in the checked environment. Python helpers require Python 3.10 or later.

These checks do not download participant data, refit models, run bootstrap or imputation analyses, or regenerate the complete set of manuscript figures. They do not establish end-to-end reproducibility or independently validate the statistical results. See [known-limitations.md](known-limitations.md) for analysis and data dependencies.
