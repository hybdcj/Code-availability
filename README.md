# Memory Low-Score Risk Model in Community-Dwelling Middle-Aged and Older Adults: Analysis Code

Analysis code for a memory low-score risk model developed in CHARLS and evaluated in HRS, ELSA, and SHARE.

This repository contains manuscript-related analysis source code, third-round revision scripts, dependency version records, and tools for preparing a reproduction workspace. It does not contain participant data, model objects, manuscript files, or generated figures and tables.

**Start with the [analysis guide](docs/analysis-map.md) and [reproduction instructions](docs/reproduction.md).** Documentation, comments, and reporting text have been translated into English. Internal script numbers are historical tracking identifiers, not study groups. This repository is not a one-command pipeline for reproducing every manuscript result: data access, path configuration, derived inputs, and some historical intermediate products must be prepared separately.

See the [English release and compatibility guide](docs/english-release.md) when reusing historical scripts or derived data with the former names.

## Repository contents

| Path | Contents |
|---|---|
| `analysis/core/` | Data preparation, longitudinal analyses, score development, internal and external validation, and sensitivity analyses |
| `analysis/revisions/` | Third-round supplementary analyses, corrections to screening metrics and analysis definitions, and figure/table source code |
| [docs/analysis-map.md](docs/analysis-map.md) | Main scripts, purposes, dependencies, and version precedence |
| [docs/data-access.md](docs/data-access.md) | Data sources and local input layout |
| [docs/reproduction.md](docs/reproduction.md) | Environment setup, path migration, and staged execution |
| [docs/known-limitations.md](docs/known-limitations.md) | Historical limitations, corrections, and remaining reproduction requirements |
| [docs/input-references.csv](docs/input-references.csv) | Locations of data/script reads, obtained by static scanning; not a complete dependency graph |
| `environment/` | Direct R dependencies, version snapshots, Python requirements, and installation script |
| `tools/` | Source verification, workspace preparation, and R syntax checking |
| [source-manifest.csv](source-manifest.csv) | Workspace paths, byte counts, and SHA-256 checksums of distributed source files |

## Quick start

```powershell
git clone https://github.com/hybdcj/Code-availability.git
cd Code-availability
python tools/verify_sources.py
Copy-Item config.example.json config.local.json
# Edit config.local.json: use a dedicated workspace and the actual local R library path.
python tools/prepare_workspace.py --config config.local.json
Rscript environment/install_dependencies.R
python -m pip install -r environment/requirements.txt
Rscript tools/check_syntax.R .
```

These commands verify and prepare code; they do not fit models. `prepare_workspace.py` copies scripts into the configured workspace and replaces known absolute paths. It refuses to overwrite files with different contents. Then prepare authorized data and derived inputs and execute the relevant modules as described in the [reproduction instructions](docs/reproduction.md).

## Analysis definitions

The primary development analysis uses CHARLS data from 2011–2015. Report the continuous model and integer score separately. The continuous model includes age, education, depression, and ADL; the nonzero items in the integer score are age, education, and depression. External validation applies the frozen model and score and distinguishes outcome threshold definitions. Consult the relevant scripts and manuscript methods for eligibility, missing-data handling, and sensitivity analyses.

Third-round corrections take precedence over historical exports: F24 addresses screening metrics; F26, the BDSI implementation; F27, validation excluding development participants; F29, internal-performance row alignment; F30, interpretation of imputation intervals; and F32, integer-score denominators in the broader population. Historical Tu 2026 comparisons remain exploratory and do not constitute complete external validation of the original model.

## Data and code availability

Individual-level data are not distributed in this repository. Access to CHARLS, HRS/RAND HRS, ELSA, and SHARE/easySHARE data must be obtained from the respective data providers under their applicable conditions. Analysis source code and reproduction documentation are available at https://github.com/hybdcj/Code-availability.

See [validation.md](docs/validation.md) for completed checks. The complete analysis has not been rerun in a clean environment, and dependency snapshots do not constitute a fully locked historical environment. No open-source license has been selected; a license and formal manuscript citation remain to be supplied by the author. Historical source versions remain available in Git history.
