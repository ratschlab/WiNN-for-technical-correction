# WiNN Technical Correction Companion Repository

This repository contains the shareable companion code for the manuscript **"WiNN enables selective correction of run-order drift and batch effects while preserving biological signal in metabolomics data"** by Tanmay Tanna and colleagues. It is organized as a publication-style analysis repository for the parts of the study that can be redistributed: the synthetic simulation benchmark and the public MetaboLights benchmark dataset.

The repository does **not** include the internal clinical FIA-MS cohort used in the manuscript. The paper's data-availability statement is explicit that the internal clinical dataset is not publicly available because of **patient privacy and consent restrictions**. For that reason, this repository only ships code and reproducible notebooks for the simulated data and the public dataset workflow.

## Repository Scope

The materials here cover the two shareable workflows described in the manuscript:

1. **Simulation benchmark**
   The notebook and script in this repository generate the synthetic LC-MS-like benchmark described in the manuscript's *Simulation benchmark* section and *Supplementary Methods S6*. The output includes:
   - a clean ground-truth intensity matrix
   - a matched matrix with injected batch and run-order artefacts
   - sample metadata
   - a feature-level design audit table

2. **Public MTBLS79 benchmark preprocessing**
   The notebook and script in this repository implement the missing-value handling described in *Supplementary Methods S7* for the public MetaboLights study **MTBLS79**. The workflow:
   - reads `Dataset07__SFPM.xlsx`
   - removes the trailing summary row present in the workbook
   - treats technical zeros as missing values
   - retains metabolites observed in at least 80% of samples
   - imputes remaining missing values with `missForest` using `set.seed(42)`
   - writes the processed feature matrix and metadata table used for downstream analysis

## Clinical Data Exclusion

The manuscript evaluates WiNN on three datasets: simulated data, a public repeated multi-batch benchmark, and an internal clinical serum metabolomics cohort. Only the first two are shareable here.

The clinical cohort is excluded from this repository for scientific and governance reasons:

- The study participants are human subjects.
- The manuscript states that the internal clinical dataset is not publicly available because of patient privacy and consent restrictions.
- Releasing raw or processed clinical measurements outside the approved access pathway would conflict with those restrictions.

This repository therefore mirrors the public part of the paper faithfully: it provides the simulation generator, the public-data preprocessing workflow, rendered notebooks, and documentation, while intentionally omitting non-public clinical data files.

## Repository Layout

```text
.
├── README.md
├── LICENSE
├── manuscript/
│   └── winn_manuscript.pdf
├── data/
│   ├── README.md
│   ├── simulated/
│   └── public/
│       ├── raw/
│       └── processed/
├── scripts/
│   ├── simulate_metabolomics_benchmark.R
│   ├── preprocess_mtbls79_public_data.R
│   └── render_notebooks.R
└── notebooks/
    ├── 01_simulation_benchmark_data.Rmd
    ├── 02_public_mtbls79_preprocessing.Rmd
    └── rendered/
```

## Requirements

The notebooks are written in R Markdown and were prepared against a recent R installation with the following package families:

- CRAN: `dplyr`, `ggplot2`, `knitr`, `missForest`, `openxlsx`, `readr`, `scales`, `tibble`, `tidyr`
- Rendering: `rmarkdown`

The notebooks will attempt to install missing CRAN packages automatically. If you prefer to manage the environment yourself, install the required packages before rendering.

## How To Use This Repository

### 1. Simulated benchmark

Generate the simulation benchmark files directly:

```r
source("scripts/simulate_metabolomics_benchmark.R")
simulate_metabolomics_benchmark()
```

This writes the benchmark files to `data/simulated/`.

### 2. Public MTBLS79 preprocessing

Download `Dataset07__SFPM.xlsx` from the public MTBLS79 record:

- MetaboLights record: [MTBLS79](https://www.ebi.ac.uk/metabolights/MTBLS79)

Place the file at:

```text
data/public/raw/Dataset07__SFPM.xlsx
```

Then run:

```r
source("scripts/preprocess_mtbls79_public_data.R")
preprocess_mtbls79_public_data()
```

This writes:

- `data/public/processed/MTBLS79_imputed_data.csv`
- `data/public/processed/MTBLS79_metadata.csv`

## Rendered Notebooks

The repository includes two rendered notebooks in `notebooks/rendered/`:

- `01_simulation_benchmark_data.nb.html`
- `02_public_mtbls79_preprocessing.nb.html`

To re-render them locally:

```bash
Rscript scripts/render_notebooks.R
```

If you want the render script to fetch the public spreadsheet automatically before rendering the public notebook, use:

```bash
Rscript scripts/render_notebooks.R --download-public-data
```

## Notes On Reproducibility

- The simulation benchmark uses fixed seeds for each stochastic component so that the ground truth, batch assignment, and drift generation remain reproducible.
- The public-data preprocessing follows the manuscript directly: the 80% feature-presence rule is applied before `missForest` imputation, and imputation is performed with `set.seed(42)`.
- The public raw spreadsheet is intentionally **not** bundled in the repository. Users should retrieve it directly from MetaboLights to preserve a clean provenance trail to the source record.

## Correspondence

**Author:** Tanmay Tanna  
**Email:** [ttanna@ethz.ch](mailto:ttanna@ethz.ch)
