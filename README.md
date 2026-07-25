# WiNN technical-correction benchmarks

This repository accompanies the WiNN manuscript and contains the shareable analysis code for the simulation, MTBLS79, BatchCorrMetabolomics Set 1, Sacurine, and WaveICA plasma benchmarks. It also contains the fixed-parameter ablation and synthetic robustness workflows used for the revision.

The full study includes a sixth, internal clinical FIA-MS cohort. Its row-level data, identifiers, and cohort-specific adapter are deliberately absent. The repository can therefore reproduce the five public or synthetic datasets, while the accompanying results release reports only aggregate clinical endpoints.

## Analysis design

Every empirical benchmark attempts the same nine methods in this order:

1. Raw
2. ComBat
3. QC-RLSC
4. QC-RFSC
5. TIGER
6. SERRF
7. WiNN auto (QC), with supplied batches
8. WiNN auto-batch (QC)
9. WiNN fixed default without QCs

A fixed set of interior QCs is hidden before any method is trained or tuned. The saved holdout files in `config/holdouts/` are shared by the method comparison and the WiNN ablation. All plotting code uses the common colour definitions in `scripts/plot_theme.R`.

## Set up the software

The current runs use R 4.5 and WiNN 0.1.3. In this version, the outlier stage computes each feature's median and MAD once from the original eligible values, identifies both tails from those fixed thresholds, and shrinks upper and lower extremes independently in one non-iterative pass. For an exact rerun, clone the package into the directory expected by the provenance checks, then install it:

```bash
git clone https://github.com/ratschlab/winn.git winn
R CMD INSTALL winn
```

Install the benchmark dependencies once:

```bash
Rscript -e "install.packages(c('BiocManager','digest','dplyr','ggplot2','jsonlite','lmtest','mgcv','missForest','openxlsx','patchwork','pmartR','qcrlscR','remotes','rmarkdown','scales','tibble','tidyr'))"
Rscript -e "BiocManager::install(c('limma','statTarget','sva'), ask=FALSE, update=FALSE)"
Rscript -e "remotes::install_github(c('pmartR/malbacR','HAN-Siyu/TIGER'), upgrade='never')"
```

## Download and preprocess the public datasets

Only small processed matrices and supporting metadata are downloaded. None of these commands downloads raw spectra.

```bash
# MTBLS79
python3 scripts/download_mtbls79.py
Rscript scripts/preprocess_mtbls79_public_data.R

# BatchCorrMetabolomics Set 1
python3 scripts/download_batchcorr_set1.py
Rscript scripts/preprocess_batchcorr_set1.R

# Sacurine / W4M00001 / MTBLS404
python3 scripts/download_sacurine.py
Rscript scripts/preprocess_sacurine.R

# WaveICA 2.0 adenocarcinoma plasma data
python3 scripts/download_waveica_adenocarcinoma.py
Rscript scripts/preprocess_waveica_adenocarcinoma.R
```

Each downloader is checksum-aware and writes a manifest. The preprocessing scripts validate matrix orientation, sample order, batch and run-order fields, QC identities, missingness, and feature variance. Details and dataset-specific limitations are in [docs/DATASETS.md](docs/DATASETS.md).

## Run the four public-data benchmarks

```bash
# MTBLS79: this notebook executes the benchmark and renders its report
Rscript -e "rmarkdown::render('notebooks/public_data_comparison.Rmd', output_dir='notebooks/rendered', clean=TRUE)"

# BatchCorr Set 1
Rscript scripts/run_batchcorr_set1_benchmark.R
Rscript -e "rmarkdown::render('notebooks/batchcorr_set1_comparison.Rmd', output_dir='notebooks/rendered', clean=TRUE)"

# Sacurine and WaveICA
Rscript scripts/run_human_public_benchmark.R --dataset=sacurine
Rscript scripts/run_human_public_benchmark.R --dataset=waveica_adenocarcinoma
Rscript -e "rmarkdown::render('notebooks/sacurine_comparison.Rmd', output_dir='notebooks/rendered', clean=TRUE)"
Rscript -e "rmarkdown::render('notebooks/waveica_adenocarcinoma_comparison.Rmd', output_dir='notebooks/rendered', clean=TRUE)"
```

Use `--force` on a command-line runner only when a deliberate full rerun is intended. Existing results are otherwise reused.

After a WiNN package update, `--winn-only` refreshes the three WiNN variants and regenerates downstream tables and figures while requiring and reusing the existing competitor caches. Reused competitor caches are read-only and are not reserialized, so their file hashes remain stable. The option cannot be combined with `--force`:

```bash
Rscript scripts/run_batchcorr_set1_benchmark.R --winn-only
Rscript scripts/run_human_public_benchmark.R --dataset=sacurine --winn-only
Rscript scripts/run_human_public_benchmark.R --dataset=waveica_adenocarcinoma --winn-only
```

The final 0.1.3 refresh archived the preceding outputs before execution and verified frozen competitor artifacts by SHA-256 before and after the WiNN-only runs.

## Run the WiNN ablation

The cumulative ablation evaluates Raw, outlier handling, selective drift correction, selective batch correction, and PQN shrinkage. A 2 × 2 gate experiment additionally compares selective versus forced-all drift and batch correction.

```bash
for dataset in simulation mtbls79 batchcorr_set1 sacurine waveica; do
  Rscript scripts/run_public_winn_ablations.R --dataset="$dataset"
done
Rscript scripts/summarize_public_winn_ablations.R
Rscript scripts/generate_public_ablation_step_impact_tables.R
```

## Synthetic benchmark and robustness analyses

The revision uses 30 independent simulation realizations and a reviewer-scoped partial-confounding grid. The partial-confounding analysis contains ten seeds, sixteen scenarios per seed, and only Raw plus fixed QC-free WiNN.

```bash
# Generate and execute the 30 canonical simulation realizations
Rscript scripts/robustness/generate_simulation_bundles.R --all
for seed in $(seq -w 1 30); do
  Rscript scripts/robustness/run_canonical_simulation_seed.R --seed-id="SIM${seed}"
done
Rscript scripts/robustness/aggregate_canonical_simulation_seeds.R

# Execute the 160 active partial-confounding scenarios
awk -F, 'NR>1 && $2 ~ /^CONF(0[1-9]|10)_/ {gsub(/\"/,"",$2); print $2}' \
  results/robustness/06_partial_confounding/full_grid/config/scenario_order.csv |
while read -r scenario; do
  Rscript scripts/robustness/run_partial_confounding_scenario.R \
    --scenario-id="$scenario" --analysis-mode=winn_only
done
Rscript scripts/robustness/aggregate_partial_confounding_winn_only_grid.R
```

These loops are intentionally simple local examples. On a scheduler, submit one seed or scenario per array task. The methods within a task are fixed by the script; the partial-confounding task does not execute comparator methods.

For a package-only refresh of an existing complete canonical simulation run, use `--winn-only`. This reruns the three WiNN modes and the WiNN ablation while requiring the six frozen Raw/competitor caches:

```bash
Rscript scripts/robustness/run_canonical_simulation_seed.R --seed-id=SIM01 --winn-only
```

## Repository map

- `scripts/`: download, preprocessing, benchmark, evaluation, and plotting code.
- `scripts/robustness/`: canonical simulation and partial-confounding workflows.
- `notebooks/`: manuscript-facing reports.
- `config/holdouts/`: frozen public hidden-QC assignments.
- `reports/`: suitability audits for the public datasets.
- `data/simulated/`: the original display simulation.
- `archive/`: superseded material retained for provenance, not used by current commands.

Generated downloads, processed matrices, caches, results, and rendered notebooks are ignored by Git. See [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) for validation and provenance conventions.
