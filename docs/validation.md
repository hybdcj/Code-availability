# Validation Record

## Original packaging checks: 2026-10-17

- SHA-256 checks verified that 68 research source files (65 R scripts, one Python script, and two variable-mapping documents) matched the original project byte for byte at the original packaging commit.
- All 67 R files, including two repository tools, passed parsing without executing analysis expressions. The Windows startup default C.UTF-8 locale was unavailable; the checker explicitly used Chinese_China.utf8.
- Python sources passed AST syntax checks. Preparation and verification helpers require Python 3.10 or later.
- Workspace preparation passed initial execution, repeated execution with identical configuration, and refusal to overwrite a modified target, in an isolated temporary directory. Original research files were not touched.
- Submitted files were checked for allowed text formats, common credential patterns, and local documentation links. CSV files contained only source manifests, input-expression indexes, and package versions.
- R versions came from installed local packages; the openpyxl version came from the packaging environment. These snapshots were not presented as historical lockfiles.

## English release checks: 2026-10-24

- All 67 R files pass syntax parsing. Static comparison with the original commit preserves token structure, numerical constants, operators, and all 324 mapped identifiers.
- All three Python source files pass AST parsing.
- The manifest verifies all 68 distributed research sources after translation. The input-reference index matches the current source lines.
- All eight renamed script paths and local Markdown links resolve. Repository filenames and text contain no remaining Chinese characters; original names in the compatibility JSON use Unicode escapes.
- Workspace preparation creates all 68 files in the expected English directory layout. Repeated preparation is idempotent, and modified destination files are correctly rejected. Known absolute project prefixes are replaced.

The English release changes source bytes, so the current manifest records the translated distribution. Original source bytes and checksums remain available in Git history. See [english-release.md](english-release.md) for historical input compatibility.

No data download, model refitting, bootstrap, imputation, or complete figure regeneration was performed. No clean-environment end-to-end reproduction is claimed. Manuscript files were not modified, and participant data and derived model objects were not uploaded. Syntax and integrity checks do not independently validate the statistical results.
