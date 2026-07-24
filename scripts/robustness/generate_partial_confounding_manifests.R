#!/usr/bin/env Rscript

# Generate the complete paired partial-confounding design without running any
# correction method.  The 320 persisted bundles contain exact phenotype
# allocations, feature truth, reconstruction seeds, diagnostics, and hashes;
# large matrices are reconstructed rather than duplicated across scenarios.

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "scripts/robustness/generate_partial_confounding_manifests.R",
    mustWork = TRUE
  )
}
script_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
source(file.path(script_dir, "simulation_core.R"), local = FALSE)
source(file.path(script_dir, "partial_confounding_core.R"), local = FALSE)

required_packages <- c("digest", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE)
}

sha256_object <- function(object) {
  digest::digest(object, algo = "sha256", serialize = TRUE)
}

relative_to <- function(path, root) {
  path <- normalizePath(path, mustWork = FALSE)
  root <- normalizePath(root, mustWork = TRUE)
  sub(paste0("^", root, "/?"), "", path)
}

write_csv_stable <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, quote = TRUE, na = "")
}

write_csv_gzip_stable <- function(x, path) {
  connection <- gzfile(path, open = "wb", compression = 9L)
  on.exit(close(connection), add = TRUE)
  utils::write.csv(x, connection, row.names = FALSE, quote = TRUE, na = "")
  close(connection)
  on.exit(NULL, add = FALSE)
}

append_context <- function(x, scenario_row) {
  cbind(
    scenario_row[rep(1L, nrow(x)), c(
      "seed_id", "scenario_id", "scenario_index", "confounding_axis",
      "nominal_strength"
    ), drop = FALSE],
    x,
    stringsAsFactors = FALSE
  )
}

validate_complete_output <- function(output_root, expected_config_sha256) {
  run_path <- file.path(output_root, "run_manifest.csv")
  files_path <- file.path(output_root, "file_manifest.csv")
  if (!file.exists(run_path) || !file.exists(files_path)) return(FALSE)
  run <- tryCatch(
    read.csv(run_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  files <- tryCatch(
    read.csv(files_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (
    is.null(run) || nrow(run) != 1L ||
      !identical(run$status[[1L]], "completed") ||
      !identical(run$config_sha256[[1L]], expected_config_sha256) ||
      is.null(files) ||
      !all(c("relative_path", "bytes", "sha256") %in% names(files))
  ) {
    return(FALSE)
  }
  paths <- file.path(output_root, files$relative_path)
  if (!all(file.exists(paths))) return(FALSE)
  if (!all(as.numeric(file.info(paths)$size) == files$bytes)) return(FALSE)
  observed_hashes <- vapply(paths, sha256_file, character(1))
  identical(unname(observed_hashes), unname(files$sha256))
}

parse_cli <- function(args) {
  output <- list(force = FALSE)
  for (arg in args) {
    if (identical(arg, "--force")) {
      output$force <- TRUE
    } else if (arg %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript scripts/robustness/",
        "generate_partial_confounding_manifests.R [--force]\n",
        sep = ""
      )
      quit(save = "no", status = 0L)
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }
  output
}

generate_partial_confounding_manifests <- function(force = FALSE) {
  output_root <- file.path(
    repo_root, "results", "robustness", "06_partial_confounding"
  )
  seed_manifest_path <- file.path(
    repo_root, "results", "robustness", "00_provenance", "seed_manifest.csv"
  )
  declared_design_path <- file.path(
    repo_root, "results", "robustness", "00_provenance",
    "partial_confounding_design.csv"
  )
  simulation_core_path <- file.path(script_dir, "simulation_core.R")
  confounding_core_path <- file.path(script_dir, "partial_confounding_core.R")
  required_inputs <- c(
    seed_manifest_path, declared_design_path, simulation_core_path,
    confounding_core_path, script_path
  )
  missing_inputs <- required_inputs[!file.exists(required_inputs)]
  if (length(missing_inputs)) {
    stop("Required input(s) missing: ", paste(missing_inputs, collapse = ", "))
  }

  seed_manifest <- read.csv(
    seed_manifest_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  seed_rows <- seed_manifest[
    seed_manifest$workstream == "partial_confounding", , drop = FALSE
  ]
  required_seed_fields <- c(
    "seed_id", "component", "seed", "prespecified_before_pilot"
  )
  missing_seed_fields <- setdiff(required_seed_fields, names(seed_rows))
  if (length(missing_seed_fields)) {
    stop(
      "Seed manifest lacks field(s): ",
      paste(missing_seed_fields, collapse = ", "), call. = FALSE
    )
  }
  expected_ids <- sprintf("CONF%02d", seq_len(20L))
  seed_rows <- seed_rows[match(expected_ids, seed_rows$seed_id), , drop = FALSE]
  if (
    nrow(seed_rows) != 20L || anyNA(seed_rows$seed_id) ||
      !identical(seed_rows$seed_id, expected_ids) ||
      anyDuplicated(seed_rows$seed_id) ||
      !all(seed_rows$component == "paired_phenotype_allocation") ||
      !all(as.logical(seed_rows$prespecified_before_pilot))
  ) {
    stop("The prespecified CONF01--CONF20 seed ledger failed validation.",
         call. = FALSE)
  }

  declared_design <- read.csv(
    declared_design_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  if (
    nrow(declared_design) != 360L ||
      !setequal(unique(declared_design$seed_id), expected_ids) ||
      !setequal(unique(declared_design$confounding_axis),
                .partial_confounding_axes) ||
      !setequal(unique(declared_design$nominal_strength),
                c(0, .partial_confounding_positive_strengths))
  ) {
    stop("The prespecified 360-row pre-dedup confounding grid is invalid.",
         call. = FALSE)
  }

  configuration <- list(
    bundle_format_version = "1.0.0",
    workstream = "06_partial_confounding",
    source_files_sha256 = list(
      seed_manifest = sha256_file(seed_manifest_path),
      declared_design = sha256_file(declared_design_path),
      simulation_core = sha256_file(simulation_core_path),
      partial_confounding_core = sha256_file(confounding_core_path),
      generator = sha256_file(script_path)
    ),
    design = list(
      master_seeds = 20L,
      scenarios_per_seed = 16L,
      total_scenarios = 320L,
      zero_strength_deduplicated = TRUE,
      axes = as.list(.partial_confounding_axes),
      strengths = as.list(c(0, .partial_confounding_positive_strengths)),
      n_features = 1000L,
      n_injections = 500L,
      n_study = 450L,
      n_qc = 50L,
      n_case = 225L,
      n_control = 225L,
      responsive_fraction = 0.10,
      responsive_count = 100L,
      effect_magnitude_log = 0.35,
      effect_parameterization =
        "centered study phenotype (-0.5,+0.5), zero for QC; case-control contrast is beta",
      maximum_batch_case_fractions = as.list(c(1, 1, 0.5, 0, 0)),
      order_alignment =
        "fixed-priority random allocation progressively swapped toward early cases",
      run_order_quantiles = "within-plate study-order quartiles",
      design_formula = "~ phenotype + plate + within_plate_position_scaled",
      phenotype_supplied_to_primary_correction_methods = FALSE,
      correction_methods_run_by_this_script = 0L
    )
  )
  config_sha256 <- sha256_object(configuration)

  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  if (validate_complete_output(output_root, config_sha256)) {
    message(
      "Existing partial-confounding manifests passed configuration, size, and ",
      "SHA-256 validation; reusing all 320 scenario bundles."
    )
    return(invisible(output_root))
  }
  generated_markers <- c(
    "scenario_manifest.csv", "phenotype_allocation.csv.gz", "run_manifest.csv",
    "file_manifest.csv", "bundles"
  )
  generated_markers <- file.path(output_root, generated_markers)
  if (any(file.exists(generated_markers)) && !force) {
    stop(
      "Existing partial-confounding output is incomplete or has a different ",
      "configuration. Re-run with --force to archive and replace only this ",
      "workstream's generated files.",
      call. = FALSE
    )
  }

  started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  temporary_root <- file.path(
    output_root, paste0(".generation-tmp-", Sys.getpid())
  )
  if (dir.exists(temporary_root)) {
    stop("Temporary output already exists: ", temporary_root, call. = FALSE)
  }
  dir.create(file.path(temporary_root, "bundles"), recursive = TRUE,
             showWarnings = FALSE)
  committed <- FALSE
  on.exit({
    if (!committed && dir.exists(temporary_root)) {
      unlink(temporary_root, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)

  scenario_manifest_rows <- list()
  allocation_rows <- list()
  feature_truth_rows <- list()
  seed_ledger_rows <- list()
  diagnostic_rows <- list()
  plate_correlation_rows <- list()
  batch_count_rows <- list()
  quantile_count_rows <- list()
  validation_rows <- list()

  for (seed_index in seq_len(nrow(seed_rows))) {
    seed_id <- seed_rows$seed_id[[seed_index]]
    master_seed <- as.integer(seed_rows$seed[[seed_index]])
    component_seeds <- expand_partial_confounding_seed(master_seed)
    seed_offsets <- component_seeds - master_seed
    seed_ledger_rows[[length(seed_ledger_rows) + 1L]] <- data.frame(
      seed_id = seed_id,
      component = c("master_seed", names(component_seeds)),
      seed = c(master_seed, unname(component_seeds)),
      offset_from_master = c(0L, unname(seed_offsets)),
      derivation = c(
        "prespecified 00_provenance/seed_manifest.csv",
        rep("master_seed + fixed component offset", length(component_seeds))
      ),
      scenario_invariant = TRUE,
      prespecified_before_pilot = TRUE,
      stringsAsFactors = FALSE
    )

    generation_seeds <- component_seeds[.canonical_required_seed_components]
    base <- generate_canonical_simulation(
      generation_seeds,
      include_artifact_matrices = TRUE
    )
    feature_truth <- partial_confounding_feature_truth(
      rownames(base$raw_intensity),
      component_seeds[["responsive_features"]],
      responsive_fraction = 0.10,
      effect_magnitude_log = 0.35
    )
    feature_truth_with_seed <- cbind(
      data.frame(seed_id = seed_id, stringsAsFactors = FALSE),
      feature_truth
    )
    feature_truth_rows[[length(feature_truth_rows) + 1L]] <-
      feature_truth_with_seed
    allocation_priority <- partial_confounding_allocation_priority(
      base$metadata, component_seeds[["allocation_priority"]]
    )
    technical_artifact <-
      base$artifact_matrices$batch_shift_by_injection +
      base$artifact_matrices$drift_by_injection
    base_hashes <- c(
      clean_baseline_log_sha256 = sha256_object(base$artifact_matrices$clean_log),
      batch_artifact_log_sha256 = sha256_object(
        base$artifact_matrices$batch_shift_by_injection
      ),
      drift_artifact_log_sha256 = sha256_object(
        base$artifact_matrices$drift_by_injection
      ),
      technical_artifact_log_sha256 = sha256_object(technical_artifact),
      metadata_sha256 = sha256_object(base$metadata),
      responsive_feature_truth_sha256 = sha256_object(feature_truth),
      allocation_priority_sha256 = sha256_object(allocation_priority),
      method_rng_policy_sha256 = sha256_object(
        component_seeds[grepl("^method_", names(component_seeds))]
      )
    )

    grid <- partial_confounding_scenario_grid(seed_id)
    seed_artifact_hashes <- character(nrow(grid))
    for (scenario_index in seq_len(nrow(grid))) {
      scenario <- grid[scenario_index, , drop = FALSE]
      allocation <- allocate_partial_confounding_phenotype(
        base$metadata,
        allocation_priority,
        scenario$confounding_axis[[1L]],
        scenario$nominal_strength[[1L]]
      )
      target_case_counts <- attr(allocation, "target_case_counts")
      use_batch_alignment <- scenario$confounding_axis[[1L]] %in%
        c("batch", "joint")
      ideal_case_counts <- if (use_batch_alignment) {
        attr(
          partial_confounding_batch_case_counts(
            scenario$nominal_strength[[1L]]
          ),
          "ideal"
        )
      } else {
        rep(45, 5L)
      }
      constructed <- construct_partial_confounding_matrices(
        base, allocation, feature_truth, include_intensities = TRUE
      )
      checks <- validate_partial_confounding_scenario(
        base, allocation, feature_truth, constructed
      )
      diagnostics <- partial_confounding_diagnostics(allocation)
      phenotype_allocation_hash <- sha256_object(
        allocation[, c("sample_id", "phenotype", "phenotype_numeric")]
      )
      clean_hash <- sha256_object(constructed$clean_ground_truth)
      raw_hash <- sha256_object(constructed$raw_intensity)
      technical_hash <- sha256_object(constructed$technical_artifact_log)
      seed_artifact_hashes[[scenario_index]] <- technical_hash
      scenario_config_sha256 <- sha256_object(list(
        global_config_sha256 = config_sha256,
        seed_id = seed_id,
        scenario = scenario,
        component_seeds = component_seeds,
        phenotype_allocation_sha256 = phenotype_allocation_hash,
        responsive_feature_truth_sha256 =
          base_hashes[["responsive_feature_truth_sha256"]]
      ))

      bundle <- list(
        bundle_format_version = "1.0.0",
        workstream = "06_partial_confounding",
        scenario = scenario,
        global_config_sha256 = config_sha256,
        scenario_config_sha256 = scenario_config_sha256,
        source_files_sha256 = configuration$source_files_sha256,
        master_seed = master_seed,
        component_seeds = component_seeds,
        generation_seeds = generation_seeds,
        correction_policy = list(
          phenotype_supplied = FALSE,
          correction_input_fields = names(constructed$correction_metadata),
          forbidden_evaluation_fields =
            .partial_confounding_forbidden_correction_fields,
          note = paste(
            "Phenotype and response truth are separate evaluation-only objects;",
            "they must never be passed to phenotype-blind correction or tuning."
          )
        ),
        correction_input_metadata = constructed$correction_metadata,
        evaluation_truth = list(
          phenotype_allocation = allocation,
          feature_truth = feature_truth
        ),
        realized_diagnostics = diagnostics$summary,
        hashes = as.list(c(
          base_hashes,
          phenotype_allocation_sha256 = phenotype_allocation_hash,
          clean_ground_truth_sha256 = clean_hash,
          raw_intensity_sha256 = raw_hash,
          scenario_technical_artifact_sha256 = technical_hash
        )),
        reconstruction = list(
          simulation_core = relative_to(simulation_core_path, repo_root),
          partial_confounding_core = relative_to(
            confounding_core_path, repo_root
          ),
          function_sequence = c(
            "generate_canonical_simulation(include_artifact_matrices=TRUE)",
            "partial_confounding_feature_truth()",
            "partial_confounding_allocation_priority()",
            "allocate_partial_confounding_phenotype()",
            "construct_partial_confounding_matrices()"
          ),
          large_matrices_persisted = FALSE
        )
      )
      bundle_path <- file.path(
        temporary_root, "bundles", paste0(scenario$scenario_id[[1L]], ".rds")
      )
      saveRDS(bundle, bundle_path, version = 3, compress = "gzip")
      bundle_sha256 <- sha256_file(bundle_path)

      allocation_output <- allocation
      allocation_output$target_batch_case_count <- unname(
        target_case_counts[match(allocation_output$plate,
                                 names(target_case_counts))]
      )
      allocation_output$target_batch_case_fraction <-
        allocation_output$target_batch_case_count / 90
      allocation_rows[[length(allocation_rows) + 1L]] <-
        append_context(allocation_output, scenario)

      diagnostics_output <- diagnostics$summary
      for (plate_index in seq_len(5L)) {
        diagnostics_output[[paste0("target_case_count_P", plate_index)]] <-
          unname(target_case_counts[[plate_index]])
        diagnostics_output[[paste0("ideal_case_count_P", plate_index)]] <-
          unname(ideal_case_counts[[plate_index]])
        diagnostics_output[[paste0("realized_case_fraction_P", plate_index)]] <-
          unname(target_case_counts[[plate_index]]) / 90
      }
      diagnostic_rows[[length(diagnostic_rows) + 1L]] <-
        append_context(diagnostics_output, scenario)
      plate_correlation_rows[[length(plate_correlation_rows) + 1L]] <-
        append_context(diagnostics$within_plate_correlations, scenario)
      batch_count_rows[[length(batch_count_rows) + 1L]] <-
        append_context(diagnostics$batch_group_counts, scenario)
      quantile_count_rows[[length(quantile_count_rows) + 1L]] <-
        append_context(diagnostics$run_order_quantile_counts, scenario)

      checks <- append_context(checks, scenario)
      validation_rows[[length(validation_rows) + 1L]] <- checks
      scenario_manifest_rows[[length(scenario_manifest_rows) + 1L]] <- cbind(
        scenario,
        diagnostics_output,
        data.frame(
          n_qc = 50L,
          n_responsive_features = 100L,
          effect_magnitude_log = 0.35,
          phenotype_supplied_to_correction = FALSE,
          correction_input_fields = paste(
            names(constructed$correction_metadata), collapse = ";"
          ),
          config_sha256 = config_sha256,
          scenario_config_sha256 = scenario_config_sha256,
          clean_baseline_log_sha256 =
            base_hashes[["clean_baseline_log_sha256"]],
          batch_artifact_log_sha256 =
            base_hashes[["batch_artifact_log_sha256"]],
          drift_artifact_log_sha256 =
            base_hashes[["drift_artifact_log_sha256"]],
          technical_artifact_log_sha256 =
            base_hashes[["technical_artifact_log_sha256"]],
          metadata_sha256 = base_hashes[["metadata_sha256"]],
          responsive_feature_truth_sha256 =
            base_hashes[["responsive_feature_truth_sha256"]],
          allocation_priority_sha256 =
            base_hashes[["allocation_priority_sha256"]],
          method_rng_policy_sha256 =
            base_hashes[["method_rng_policy_sha256"]],
          phenotype_allocation_sha256 = phenotype_allocation_hash,
          clean_ground_truth_sha256 = clean_hash,
          raw_intensity_sha256 = raw_hash,
          bundle_relative_path = file.path(
            "bundles", paste0(scenario$scenario_id[[1L]], ".rds")
          ),
          bundle_sha256 = bundle_sha256,
          validation_passed = all(checks$passed),
          stringsAsFactors = FALSE
        ),
        stringsAsFactors = FALSE
      )
      rm(constructed, bundle)
    }

    paired_hash_check <- length(unique(seed_artifact_hashes)) == 1L &&
      identical(
        unique(seed_artifact_hashes),
        base_hashes[["technical_artifact_log_sha256"]]
      )
    validation_rows[[length(validation_rows) + 1L]] <- data.frame(
      seed_id = seed_id,
      scenario_id = paste0(seed_id, "_paired_seed_invariant"),
      scenario_index = NA_integer_,
      confounding_axis = "all",
      nominal_strength = NA_real_,
      check = "technical_artifact_hash_identical_across_16_scenarios",
      observed = paste(unique(seed_artifact_hashes), collapse = ";"),
      passed = paired_hash_check,
      stringsAsFactors = FALSE
    )
    if (!paired_hash_check) {
      stop(seed_id, " did not preserve an identical technical-artifact matrix.")
    }
    message(seed_id, ": generated and validated 16 paired scenario bundles.")
    rm(base, feature_truth, technical_artifact)
    invisible(gc(FALSE))
  }

  scenario_manifest <- do.call(rbind, scenario_manifest_rows)
  phenotype_allocation <- do.call(rbind, allocation_rows)
  feature_truth_manifest <- do.call(rbind, feature_truth_rows)
  seed_ledger <- do.call(rbind, seed_ledger_rows)
  design_diagnostics <- do.call(rbind, diagnostic_rows)
  within_plate_correlations <- do.call(rbind, plate_correlation_rows)
  batch_group_counts <- do.call(rbind, batch_count_rows)
  run_order_quantile_counts <- do.call(rbind, quantile_count_rows)
  validation_checks <- do.call(rbind, validation_rows)
  rownames(scenario_manifest) <- NULL
  rownames(phenotype_allocation) <- NULL
  rownames(feature_truth_manifest) <- NULL
  rownames(seed_ledger) <- NULL
  rownames(design_diagnostics) <- NULL
  rownames(within_plate_correlations) <- NULL
  rownames(batch_group_counts) <- NULL
  rownames(run_order_quantile_counts) <- NULL
  rownames(validation_checks) <- NULL

  global_validation <- data.frame(
    check = c(
      "20_master_seeds",
      "16_deduplicated_scenarios_per_seed",
      "320_total_scenarios",
      "160000_injection_allocation_rows",
      "all_scenario_validations_pass",
      "one_artifact_hash_per_master_seed",
      "phenotype_never_in_correction_schema",
      "all_designs_full_rank",
      "all_design_condition_numbers_finite",
      "all_bundle_files_exist"
    ),
    observed = c(
      length(unique(scenario_manifest$seed_id)),
      paste(sort(unique(table(scenario_manifest$seed_id))), collapse = ","),
      nrow(scenario_manifest),
      nrow(phenotype_allocation),
      sum(!validation_checks$passed),
      paste(sort(unique(tapply(
        scenario_manifest$technical_artifact_log_sha256,
        scenario_manifest$seed_id,
        function(x) length(unique(x))
      ))), collapse = ","),
      sum(scenario_manifest$phenotype_supplied_to_correction),
      sum(!design_diagnostics$design_full_rank),
      sum(!is.finite(design_diagnostics$design_condition_number)),
      sum(!file.exists(file.path(
        temporary_root, scenario_manifest$bundle_relative_path
      )))
    ),
    expected = c("20", "16", "320", "160000", "0", "1", "0", "0", "0", "0"),
    passed = c(
      length(unique(scenario_manifest$seed_id)) == 20L,
      all(table(scenario_manifest$seed_id) == 16L),
      nrow(scenario_manifest) == 320L,
      nrow(phenotype_allocation) == 160000L,
      all(validation_checks$passed),
      all(tapply(
        scenario_manifest$technical_artifact_log_sha256,
        scenario_manifest$seed_id,
        function(x) length(unique(x)) == 1L
      )),
      !any(scenario_manifest$phenotype_supplied_to_correction),
      all(design_diagnostics$design_full_rank),
      all(is.finite(design_diagnostics$design_condition_number)),
      all(file.exists(file.path(
        temporary_root, scenario_manifest$bundle_relative_path
      )))
    ),
    stringsAsFactors = FALSE
  )
  if (!all(global_validation$passed)) {
    stop(
      "Global validation failed: ",
      paste(global_validation$check[!global_validation$passed], collapse = ", "),
      call. = FALSE
    )
  }

  correction_input_schema <- data.frame(
    field = c(
      .partial_confounding_correction_fields,
      "phenotype", "phenotype_numeric", "responsive",
      "effect_direction", "phenotype_effect_log"
    ),
    supplied_to_primary_correction = c(
      rep(TRUE, length(.partial_confounding_correction_fields)),
      rep(FALSE, 5L)
    ),
    role = c(
      "alignment identifier", "supplied categorical batch", "supplied order",
      "supplied global order", "QC identity/training role",
      rep("evaluation-only biological truth", 5L)
    ),
    stringsAsFactors = FALSE
  )

  write_csv_stable(scenario_manifest,
                   file.path(temporary_root, "scenario_manifest.csv"))
  write_csv_gzip_stable(
    phenotype_allocation,
    file.path(temporary_root, "phenotype_allocation.csv.gz")
  )
  write_csv_stable(feature_truth_manifest,
                   file.path(temporary_root, "responsive_feature_truth.csv"))
  write_csv_stable(seed_ledger,
                   file.path(temporary_root, "seed_ledger.csv"))
  write_csv_stable(design_diagnostics,
                   file.path(temporary_root, "design_diagnostics.csv"))
  write_csv_stable(
    within_plate_correlations,
    file.path(temporary_root, "within_plate_order_correlations.csv")
  )
  write_csv_stable(batch_group_counts,
                   file.path(temporary_root, "batch_group_counts.csv"))
  write_csv_stable(
    run_order_quantile_counts,
    file.path(temporary_root, "run_order_quantile_counts.csv")
  )
  write_csv_stable(validation_checks,
                   file.path(temporary_root, "validation_checks.csv"))
  write_csv_stable(global_validation,
                   file.path(temporary_root, "global_validation.csv"))
  write_csv_stable(
    correction_input_schema,
    file.path(temporary_root, "correction_input_schema.csv")
  )
  jsonlite::write_json(
    configuration,
    file.path(temporary_root, "design_configuration.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA
  )

  readme <- c(
    "# Paired partial-confounding simulation design",
    "",
    paste(
      "This directory contains the validated construction layer for 20",
      "prespecified master seeds and 320 paired scenarios. No correction method",
      "was run by this generator."
    ),
    "",
    "## Design",
    "",
    paste(
      "Each seed has one deduplicated strength-zero scenario plus batch,",
      "run-order, and joint scenarios at strengths 0.25, 0.50, 0.75, 0.90,",
      "and 1.00 (16 scenarios per seed). All scenarios retain 225 cases and",
      "225 controls among 450 study injections; 50 pooled QCs are phenotype-neutral."
    ),
    paste(
      "Exactly 100 of 1,000 features are responsive, with 50 positive and 50",
      "negative case-control log contrasts of magnitude 0.35. Biological effects",
      "are added to the clean log matrix before the fixed batch and drift artifacts."
    ),
    paste(
      "The maximum batch pattern is (1, 1, 0.5, 0, 0) case fraction across",
      "P1--P5. Run-order alignment uses deterministic swaps from a fixed random",
      "allocation toward cases occurring early within each plate."
    ),
    "",
    "## Key files",
    "",
    "- `scenario_manifest.csv`: one row per scenario with realized confounding, design diagnostics, seeds/configuration hashes, matrix object hashes, and bundle hashes.",
    "- `phenotype_allocation.csv.gz`: exact injection-level allocation for every scenario.",
    "- `batch_group_counts.csv` and `run_order_quantile_counts.csv`: exact balance diagnostics.",
    "- `design_diagnostics.csv`: Cramer's V, absolute phenotype/order correlation, rank, and condition number.",
    "- `bundles/*.rds`: lightweight reconstruction bundles; raw and clean matrices are intentionally not duplicated 320 times.",
    "- `correction_input_schema.csv`: explicit separation of correction inputs from evaluation-only phenotype truth.",
    "- `validation_checks.csv` and `global_validation.csv`: machine-readable invariant checks.",
    "",
    "## Phenotype-blind correction contract",
    "",
    paste(
      "Primary correction methods may receive only sample ID, supplied plate,",
      "acquisition order, and sample type/QC role. Phenotype labels and responsive-",
      "feature truth are separate evaluation objects and are never correction or",
      "tuning inputs."
    ),
    "",
    "## Reproduction",
    "",
    "Run `Rscript scripts/robustness/generate_partial_confounding_manifests.R` from the repository root. A matching completed output is reused only after file-size and SHA-256 validation."
  )
  writeLines(readme, file.path(temporary_root, "README.md"), useBytes = TRUE)

  completed_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  run_manifest <- data.frame(
    run_id = paste0("PARTIAL_CONFOUNDING_DESIGN_", format(Sys.time(), "%Y%m%d")),
    workstream = "06_partial_confounding",
    status = "completed",
    started_at = started_at,
    completed_at = completed_at,
    config_sha256 = config_sha256,
    master_seeds = 20L,
    scenarios_per_seed = 16L,
    scenario_units = 320L,
    correction_methods_run = 0L,
    validation_failures = sum(!validation_checks$passed) +
      sum(!global_validation$passed),
    note = paste(
      "Strength-zero axes deduplicated from the prespecified 360-row grid;",
      "all exact allocations, realized confounding diagnostics, reconstruction",
      "seeds, and hashes persisted."
    ),
    stringsAsFactors = FALSE
  )
  write_csv_stable(run_manifest,
                   file.path(temporary_root, "run_manifest.csv"))

  persisted_files <- list.files(
    temporary_root, recursive = TRUE, full.names = TRUE,
    all.files = FALSE, include.dirs = FALSE
  )
  persisted_files <- persisted_files[
    basename(persisted_files) != "file_manifest.csv"
  ]
  persisted_files <- sort(persisted_files)
  file_manifest <- data.frame(
    relative_path = vapply(
      persisted_files, relative_to, character(1), root = temporary_root
    ),
    bytes = as.numeric(file.info(persisted_files)$size),
    sha256 = vapply(persisted_files, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  write_csv_stable(file_manifest,
                   file.path(temporary_root, "file_manifest.csv"))

  if (any(file.exists(generated_markers))) {
    archive_root <- file.path(
      output_root, "archive",
      paste0("generation-", format(Sys.time(), "%Y%m%dT%H%M%S"))
    )
    dir.create(archive_root, recursive = TRUE, showWarnings = FALSE)
    for (path in generated_markers[file.exists(generated_markers)]) {
      if (!file.rename(path, file.path(archive_root, basename(path)))) {
        stop("Could not archive prior generated output: ", path, call. = FALSE)
      }
    }
  }

  top_level <- list.files(
    temporary_root, full.names = TRUE, all.files = TRUE, no.. = TRUE
  )
  for (path in top_level) {
    destination <- file.path(output_root, basename(path))
    if (file.exists(destination)) {
      stop("Refusing to overwrite unexpected output: ", destination,
           call. = FALSE)
    }
    if (!file.rename(path, destination)) {
      stop("Could not commit generated output: ", destination, call. = FALSE)
    }
  }
  if (!unlink(temporary_root, recursive = TRUE, force = TRUE) == 0L) {
    warning("Could not remove empty temporary directory: ", temporary_root)
  }
  committed <- TRUE

  if (!validate_complete_output(output_root, config_sha256)) {
    stop("Committed partial-confounding output failed SHA-256 validation.",
         call. = FALSE)
  }
  message(
    "Completed 320/320 lightweight partial-confounding bundles; all ",
    nrow(validation_checks) + nrow(global_validation),
    " scenario/global invariant checks passed."
  )
  invisible(output_root)
}

if (sys.nframe() == 0L) {
  setwd(repo_root)
  cli <- parse_cli(commandArgs(trailingOnly = TRUE))
  generate_partial_confounding_manifests(force = cli$force)
}
