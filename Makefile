R ?= Rscript
PROJECT_LIB := $(CURDIR)/Rlib
export R_LIBS_USER := $(PROJECT_LIB)
export WINN_RELEASE_ROOT := $(CURDIR)
export WINN_CANONICAL_SOURCE_ROOT := $(CURDIR)

.PHONY: restore install-winn download preprocess simulations smoke aggregate-historical validate-results report tables config-manifest

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

aggregate-historical:
	$(R) analysis/finalize_reuse_audit.R
	$(R) analysis/aggregate_release.R
	$(R) analysis/build_execution_ledger.R
	$(R) analysis/audit_reference_separation.R
	$(R) analysis/validate_release.R

validate-results:
	python3 analysis/validate_publication_bundle.py

config-manifest:
	$(R) analysis/write_config_manifest.R

report:
	$(R) manuscript/scripts/generate_revision_main_panels.R
	$(R) manuscript/scripts/generate_revision_supplementary_uncertainty_panels.R
	$(R) manuscript/scripts/generate_revision_s6_panels.R
	python3 manuscript/scripts/assemble_main_figures.py
	python3 manuscript/scripts/assemble_supplementary_figures.py

tables:
	$(R) manuscript/scripts/generate_supplementary_tables.R --final-root=results/final/revision_inputs/supplement_20260727
	$(R) manuscript/scripts/nearest_baseline/generate_validation_tables.R
	python3 manuscript/scripts/generate_review2_support_tables.py
