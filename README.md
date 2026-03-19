# WiNN Technical Correction Companion Repository

Companion code repository for the manuscript *WiNN enables selective correction of run-order drift and batch effects while preserving biological signal in metabolomics data*.

Author: Tanmay Tanna (`ttanna@ethz.ch`)

## Contents

- `manuscript/winn_manuscript.pdf`: manuscript PDF.
- `data/simulated/`: simulated benchmark data.
- `data/public/`: public MTBLS79 input/output locations.
- `notebooks/simulation_comparison.Rmd`: simulated benchmark analysis.
- `notebooks/public_data_comparison.Rmd`: public MTBLS79 benchmark analysis.
- `notebooks/rendered/`: rendered notebook outputs.
- `scripts/`: data-generation and public-data preprocessing utilities.

## Data Preparation

Populate the data directories before rendering the notebooks.

### Simulated benchmark

```bash
Rscript scripts/simulate_metabolomics_benchmark.R data/simulated
```

### Public MTBLS79 benchmark

1. Download `Dataset07__SFPM.xlsx` from the [MetaboLights MTBLS79 record](https://www.ebi.ac.uk/metabolights/MTBLS79).
2. Save it at `data/public/raw/Dataset07__SFPM.xlsx`.
3. Preprocess it with:

```bash
Rscript scripts/preprocess_mtbls79_public_data.R \
  data/public/raw/Dataset07__SFPM.xlsx \
  data/public/processed
```

## Analysis

Run the analyses through the notebooks after the data are available. The notebooks handle dependency installation for the comparison methods, perform subset-based tuning where applicable, and then execute full-data correction and evaluation.

```bash
Rscript -e "rmarkdown::render('notebooks/simulation_comparison.Rmd', output_dir = 'notebooks/rendered')"
Rscript -e "rmarkdown::render('notebooks/public_data_comparison.Rmd', output_dir = 'notebooks/rendered')"
```

`TIGER` and `SERRF` are the main runtime bottlenecks.
