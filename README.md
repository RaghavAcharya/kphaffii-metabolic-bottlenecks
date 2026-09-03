# Metabolic Upregulation (Paper 1) — Code & Data

Computational code and data accompanying the manuscript
**"Genome-Scale Metabolic Models Identify Metabolic Bottlenecks in Production of
Recombinant Proteins in *Komagataella phaffii*"** (Acharya et al.).

This repository reproduces the modeling and figure pipeline of the paper: adapting
the *K. phaffii* GEM iMT1026v3 to produce four recombinant proteins (VHH, RBD, HSA,
IgG1), identifying overexpression/knockout targets by FSEOF, nominating transcription-
factor targets, and computing per-gene nutrient cost.

> Extracted from the larger `Kpha_GEM_Expansion` research repository; only the files
> required for this paper are included.

## Requirements

| Tool | Version | Used for |
|------|---------|----------|
| MATLAB | R2024a | model build, FBA/FSEOF, cost pipeline |
| Gurobi | 9.5+ | LP/MILP solver for MATLAB |
| COBRA Toolbox | v3.0 | constraint-based modeling |
| RAVEN Toolbox | 2.0 | model I/O (`importModel`) and utilities |
| R | 4.5.0 | figure generation and statistics |

R packages: `tidyverse`, `readxl`, `janitor`, `ggplot2`, `fgsea`/`gmtPathways`, `patchwork`.

## Setup

1. Install MATLAB, Gurobi, COBRA v3.0, and RAVEN 2.0.
2. Open **`setup_solver.m`** and edit it to point to your Gurobi install (and, if not
   already on your MATLAB path, COBRA and RAVEN). Each MATLAB script calls this
   automatically and adds `matlab/functions/` to the path — you do not need to set
   paths per script.
3. Run any script from `matlab/scripts/**` directly (the scripts locate the repo root
   themselves via `mfilename`, so the working directory does not matter).

## Layout

```
Paper1_Code_and_Data/
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── setup_solver.m                  ← edit this once (Gurobi/COBRA/RAVEN paths)
├── matlab/
│   ├── scripts/
│   │   ├── 01_HSA_target_identification/       HSA GEM build + validation + FSEOF (Fig 1, Fig 2a,b)
│   │   ├── 02_common_target_identification/    4-protein GEMs + FSEOF + common targets (Fig 1, Fig 2c,d)
│   │   ├── 03_transcription_factor_targets/    TF–gene network (Fig 3)
│   │   └── 04_nutrient_cost/                   per-gene nutrient cost pipeline (Fig 4, SI §5)
│   └── functions/modelBuilding/    17 custom reaction builders required by the scripts
│                                    (translation, N-/O-glyc, disulfide, COPII transport, PTMs)
├── models/                         base GEM iMT1026v3 (see Provenance)
├── data/                           model-modification + gene-property inputs
├── R/                              figure & statistics notebooks
│   ├── Figs.Rmd                    master figure notebook → publication SVGs
│   ├── geneTargetAnalysis.Rmd      gene-target / GO analysis (Fig 2 support)
│   ├── 03_TF_Graphs.Rmd            TF graphs (Fig 3 support)
│   └── figure_data/                inputs the R notebooks read
└── verification/                   transfer-integrity checks (see below)
```

## Run order and figure map

MATLAB (produces the modeling outputs and intermediate CSVs):

| Script | Produces | Manuscript |
|--------|----------|-----------|
| `matlab/scripts/01_HSA_target_identification/HSA_model_Generation_Verification.m` | HSA GEM, carbon-source validation, FSEOF slopes, gene-target table | Fig 1 (HSA), Fig 2a,b |
| `matlab/scripts/02_common_target_identification/common_target_identification_validation.m` | four-protein GEMs, cross-protein FSEOF (`FSEOF_results_*.csv`), `commonTargets.csv` | Fig 1, Fig 2c,d |
| `matlab/scripts/03_transcription_factor_targets/tf_target_analysis.m` | TF→gene directed-graph figures | Fig 3 |
| `matlab/scripts/04_nutrient_cost/metCostDetermination.m` | per-gene nutrient/protein cost (`synthCost_*.csv`) | Fig 4, SI §5 |
| `matlab/scripts/04_nutrient_cost/metCostKOBenchmark.m` | cost-based KO benchmark | Fig 4 support |

R (renders the final publication figures from the MATLAB/experimental outputs):

| Notebook | Produces |
|----------|----------|
| `R/Figs.Rmd` | Fig 1c/d, Fig 2abc/d, Fig 3b/c/d, Fig 4/4b, Fig 5, and SI cost figures (SVG/PDF, written into `R/`) |
| `R/geneTargetAnalysis.Rmd` | gene-target enrichment/analysis supporting Fig 2 |
| `R/03_TF_Graphs.Rmd` | transcription-factor graphs supporting Fig 3 |

FSEOF itself is implemented inline in scripts 01 and 02 (local function `FSEOF_custom`),
following Choi et al. 2010.

## Data & model provenance

- `models/MODEL1612130000.xml`, `models/mbt212871sup0003AppendixS3.xml` — the base
  *K. phaffii* GEM **iMT1026v3** (Tomàs-Gamisans et al.), two source encodings from
  the EBI BioModels entry <https://www.ebi.ac.uk/biomodels/MODEL1612130000>. Script 01
  loads the first, script 02 the second; they are the same model. Redistribution here
  is under the source's terms — verify licensing before publishing.
- `models/iMT1026v3.mat` — the same model as a MATLAB struct (used by the cost pipeline).
- `data/baselineModifications.xlsx` — metabolites and reactions added to the base model
  (two sheets: *New Metabolites*, *New Reactions*).
- `data/ppPSIMv3_RP.csv` — per-gene protein-property table (sequence, localization,
  PTM counts, signal-peptide/GPI/TM flags) consolidated per SI §5.1.
- `data/geneConversion.csv` — GQ67 ↔ PAS locus-tag map with functional annotation.
- `data/essentiality.csv` — gene-essentiality reference used to annotate targets.

## Verification

After copying/cloning, confirm the transfer is complete (see `verification/`):

- **`verification/verify_files.py`** — cross-platform file/reference check. Run:
  `python verification/verify_files.py` from the repo root. Confirms every expected
  script, function, model, data, and R input is present, and that each path referenced
  by the scripts resolves inside the repo.
- **`verification/verify_matlab_dependencies.m`** — authoritative dependency check.
  Run in MATLAB from the repo root. Uses `matlab.codetools.requiredFilesAndProducts`
  on each entry script to compute the true required-file set, flags anything required
  but missing, reports the external toolboxes each script needs (COBRA/RAVEN/Gurobi),
  and lists any bundled function not actually required (safe-to-prune candidates).

## Citation

If you use this code or data, please cite:

> Acharya, R.; Cha, E.; Atufa, A.; Wang, E.; Narayanan, H.; Sunday, B. C.; Love, J. C. *Genome-Scale Metabolic Models Identify Metabolic Bottlenecks in Production of Recombinant Proteins in Komagataella phaffii.* (2026).

Machine-readable metadata is in `CITATION.cff`. Released under the terms in `LICENSE`.

## Notes

- `matlab/functions/` contains exactly the 17 functions required by the four entry
  scripts, as confirmed by `matlab.codetools.requiredFilesAndProducts` (see
  `verification/verify_matlab_dependencies.m`). Enzyme/proteome-constraint helpers that
  appeared only in a comment in the synthesis-reaction builder were confirmed unused and
  removed during pruning.
- The `*ProductionFluxData_*.csv` and `synthCost_*.csv` files in `04_nutrient_cost/` are
  regenerated by the scripts; shipped copies are provided only for convenience.
- COBRA, RAVEN, and Gurobi are **not** bundled — install them separately and configure
  `setup_solver.m`.
- `.gitignore` excludes regenerated pipeline outputs (`synthCost_*.csv`,
  `*ProductionFluxData_*.csv`, `fluxData.mat`) and figures rendered by `R/Figs.Rmd`
  (`R/*.svg`, `R/*.pdf`, `R/*.png`); copies currently in the folder are provided for
  convenience and are not tracked by git.
