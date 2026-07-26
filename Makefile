R ?= Rscript
PROJECT_LIB := $(CURDIR)/Rlib
export R_LIBS_USER := $(PROJECT_LIB)
export WINN_RELEASE_ROOT := $(CURDIR)
export WINN_CANONICAL_SOURCE_ROOT := $(CURDIR)

.PHONY: restore install-winn download preprocess simulations smoke aggregate validate-results report config-manifest

restore:
	$(R) -e 'renv::restore(prompt = FALSE)'

install-winn:
	mkdir -p "$(PROJECT_LIB)"
	R CMD INSTALL --library="$(PROJECT_LIB)" package/winn_0.1.4.tar.gz

download:
	python3 scripts/download_mtbls79.py
	python3 scripts/download_batchcorr_set1.py
	python3 scripts/download_sacurine.py
	python3 scripts/download_waveica_adenocarcinoma.py

preprocess:
	$(R) scripts/preprocess_mtbls79_public_data.R
	$(R) scripts/preprocess_batchcorr_set1.R
	$(R) scripts/preprocess_sacurine.R
	$(R) scripts/preprocess_waveica_adenocarcinoma.R

simulations:
	$(R) scripts/robustness/generate_simulation_bundles.R --all

smoke:
	$(R) analysis/smoke_public_workflow.R

aggregate:
	$(R) analysis/finalize_reuse_audit.R
	$(R) analysis/aggregate_release.R
	$(R) analysis/build_execution_ledger.R
	$(R) analysis/audit_reference_separation.R
	$(R) analysis/validate_release.R

validate-results:
	$(R) analysis/write_config_manifest.R --check
	$(R) analysis/validate_committed_results.R

config-manifest:
	$(R) analysis/write_config_manifest.R

report:
	mkdir -p notebooks/rendered
	$(R) -e 'rmarkdown::render("notebooks/benchmark_summary.Rmd", output_dir="notebooks/rendered", clean=TRUE, envir=new.env(parent=globalenv()))'
