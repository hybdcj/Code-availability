# English Release and Historical Compatibility

The 2026-10-24 release translates the repository documentation, source comments, reporting text, Chinese identifiers, and eight script filenames into English. Numerical constants, statistical operators, and analysis control flow are preserved. The original sources remain in Git history at commit `6af6d074369f9a3c2c4b40320528b5edb3776e79`.

## Name mapping

[legacy-name-map.json](legacy-name-map.json) records the original-to-English mappings. Original names are stored as JSON Unicode escapes; a standard JSON reader restores them automatically.

| Key | Purpose |
|---|---|
| `repository_files` | Renamed repository scripts |
| `path_components` | Historical directory and file components |
| `identifiers` | Data-frame columns and other R identifiers |
| `literals` | Changed labels, reporting strings, and path strings |

For example, the revised row-alignment audit is `analysis/revisions/F29_fix_internal_performance_order.R`. The workspace preparation tool places revision scripts in `<project_root>/16.manuscript_revision_round_three/_audit/scripts/` and creates the corresponding `output/` directory.

## Preparing inputs

Use a new workspace and the translated scripts when regenerating derived inputs. Raw cohort variable codes, dataset names in Latin characters, model coefficients, and numerical thresholds have not been translated. Data access requirements are unchanged.

Historical derived CSV/RDS files may contain the former column names, factor levels, or reporting labels. They are not automatically interchangeable with newly generated English outputs. Before reusing a historical cache, consult the mapping and translate a separate copy of the required names and labels, including nested data frames and factor levels where applicable. Do not alter numerical observations or silently overwrite original research data. Model objects and formulas require separate inspection; this release does not claim that renaming their fields is sufficient for migration.

Some Latin-character identifiers and filenames retain historical suffixes such as `_cn` or `_zh` so existing script references remain recognizable. Their repository text has been translated into English. Windows Chinese locale settings and historical font names are retained for compatibility with the original workflow.

## Verification scope

See [validation.md](validation.md). Static comparison checks preserve R token structure, numerical constants, operators, and mapped identifiers; they do not prove statistical equivalence after execution. No participant data, fitted models, or generated manuscript products were needed for this release.
