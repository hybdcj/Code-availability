# Memory Low-Score Risk Model in Community-Dwelling Middle-Aged and Older Adults: Analysis Code

Analysis code for a memory low-score risk model developed in CHARLS and evaluated in HRS, ELSA, and SHARE.

This repository contains analysis source code, supplementary analyses, dependency information, and tools for preparing a reproduction workspace. It does not contain participant data, model objects, manuscript files, or generated figures and tables.

**Start with the [analysis guide](docs/analysis-map.md) and [reproduction instructions](docs/reproduction.md).** Script numbers identify analysis modules, not study groups. Data access, path configuration, derived inputs, and intermediate products must be prepared before running the relevant modules.

## Repository contents

| Path | Contents |
|---|---|
| `analysis/core/` | Data preparation, longitudinal analyses, score development, internal and external validation, and sensitivity analyses |
| `analysis/revisions/` | Supplementary analyses, screening-metric audits, and figure/table source code |
| [docs/analysis-map.md](docs/analysis-map.md) | Main scripts, purposes, dependencies, and version precedence |
| [docs/data-access.md](docs/data-access.md) | Data sources and local input layout |
| [docs/reproduction.md](docs/reproduction.md) | Environment setup, path migration, and staged execution |
| [docs/known-limitations.md](docs/known-limitations.md) | Analysis limitations and reproduction requirements |
| [docs/input-references.csv](docs/input-references.csv) | Locations of data/script reads, obtained by static scanning; not a complete dependency graph |
| [docs/f39-reproduction.md](docs/f39-reproduction.md) | F39 source provenance and S37/S38 subgroup and missingness generation workflow |
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

Use F24 for screening metrics, F26 for the BDSI implementation, F27 for validation excluding development participants, F29 for internal-performance row alignment, F30 for imputation intervals, and F32 for integer-score denominators in the broader population. Tu 2026 comparisons are exploratory and do not constitute complete external validation of the original model.

## Data and code availability

Individual-level data are not distributed in this repository. Access to CHARLS, HRS/RAND HRS, ELSA, and SHARE/easySHARE data must be obtained from the respective data providers under their applicable conditions. Analysis source code and reproduction documentation are available at https://github.com/hybdcj/Code-availability.

See [validation.md](docs/validation.md) for software checks and their scope. A complete analysis run in a clean environment has not been verified. Dependency snapshots do not constitute a fully locked environment. Licensing and formal manuscript citation information are not yet provided.
