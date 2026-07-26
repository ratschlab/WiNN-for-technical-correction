#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("scripts/robustness/generate_simulation_bundles.R", mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
source(file.path(dirname(script_path), "simulation_core.R"), local = FALSE)
source(file.path(dirname(script_path), "canonical_cache.R"), local = FALSE)

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

portable_matrix_sha256 <- function(x, digits = 9L) {
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x))) {
    stop("Portable simulation hashes require a finite numeric matrix.", call. = FALSE)
  }
  canonical_sha256_object(round(x, digits = digits))
}

write_csv_stable <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, quote = TRUE, na = "")
}

read_seed_manifest <- function(path) {
  if (!file.exists(path)) stop("Seed manifest is missing: ", path, call. = FALSE)
  manifest <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c(
    "workstream", "seed_id", "display_realization", "component", "seed",
    "prespecified_before_pilot"
  )
  missing <- setdiff(required, names(manifest))
  if (length(missing)) {
    stop("Seed manifest lacks column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  manifest
}

read_bundle_hash_manifest <- function(path) {
  if (!file.exists(path)) stop("Simulation hash manifest is missing: ", path, call. = FALSE)
  manifest <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c(
    "seed_id", "object", "features", "injections",
    "reference_object_sha256", "portable_round9_sha256", "portable_digits"
  )
  missing <- setdiff(required, names(manifest))
  if (length(missing)) {
    stop("Simulation hash manifest lacks column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (nrow(manifest) != 20L || anyDuplicated(manifest[c("seed_id", "object")])) {
    stop("Simulation hash manifest must contain one raw and one truth row for SIM01-SIM10.", call. = FALSE)
  }
  manifest
}

validate_simulation_identity <- function(seed_id, simulation, hash_manifest, require_exact = FALSE) {
  objects <- list(
    raw_intensity = simulation$raw_intensity,
    clean_ground_truth = simulation$clean_ground_truth
  )
  checks <- lapply(names(objects), function(object_name) {
    expected <- hash_manifest[
      hash_manifest$seed_id == seed_id & hash_manifest$object == object_name,
      ,
      drop = FALSE
    ]
    if (nrow(expected) != 1L) {
      stop("Missing unique simulation hash record for ", seed_id, "/", object_name, call. = FALSE)
    }
    x <- objects[[object_name]]
    observed_object <- sha256_object(x)
    observed_portable <- portable_matrix_sha256(x, expected$portable_digits[[1L]])
    passed <- nrow(x) == expected$features[[1L]] &&
      ncol(x) == expected$injections[[1L]] &&
      identical(observed_portable, expected$portable_round9_sha256[[1L]]) &&
      (!require_exact || identical(observed_object, expected$reference_object_sha256[[1L]]))
    data.frame(
      object = object_name,
      observed_object_sha256 = observed_object,
      reference_object_sha256 = expected$reference_object_sha256[[1L]],
      observed_portable_round9_sha256 = observed_portable,
      expected_portable_round9_sha256 = expected$portable_round9_sha256[[1L]],
      portable_digits = expected$portable_digits[[1L]],
      exact_object_hash_required = require_exact,
      passed = passed,
      stringsAsFactors = FALSE
    )
  })
  checks <- do.call(rbind, checks)
  if (!all(checks$passed)) {
    stop(seed_id, ": generated matrices do not match the frozen simulation identity.", call. = FALSE)
  }
  checks
}

simulation_ledger_for <- function(seed_manifest, seed_id) {
  rows <- seed_manifest[
    seed_manifest$workstream == "canonical_repeated_simulation" &
      seed_manifest$seed_id == seed_id,
    ,
    drop = FALSE
  ]
  if (!nrow(rows)) stop("No seed ledger rows found for ", seed_id, call. = FALSE)
  if (anyDuplicated(rows$component)) {
    stop("Duplicate component rows in seed ledger for ", seed_id, call. = FALSE)
  }
  rows$seed <- as.integer(rows$seed)
  rows
}

validate_existing_bundle <- function(bundle_dir, expected_config_hash) {
  file_manifest_path <- file.path(bundle_dir, "bundle_files.csv")
  provenance_path <- file.path(bundle_dir, "bundle_provenance.json")
  if (!file.exists(file_manifest_path) || !file.exists(provenance_path)) return(FALSE)

  provenance <- tryCatch(
    jsonlite::read_json(provenance_path, simplifyVector = TRUE),
    error = function(e) NULL
  )
  if (is.null(provenance) || !identical(provenance$config_sha256, expected_config_hash)) {
    return(FALSE)
  }

  files <- tryCatch(
    read.csv(file_manifest_path, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(files) || !all(c("relative_path", "bytes", "sha256") %in% names(files))) {
    return(FALSE)
  }
  paths <- file.path(bundle_dir, files$relative_path)
  if (!all(file.exists(paths))) return(FALSE)
  if (!all(as.numeric(file.info(paths)$size) == files$bytes)) return(FALSE)
  observed <- vapply(paths, sha256_file, character(1))
  identical(unname(observed), unname(files$sha256))
}

verify_committed_primary <- function(output_root, hash_manifest) {
  matrix_path <- file.path(output_root, "SIM01", "raw_intensity.rds")
  truth_path <- file.path(output_root, "SIM01", "clean_ground_truth.rds")
  manifest_path <- file.path(repo_root, "analysis", "config", "dataset_manifest.csv")
  if (!file.exists(matrix_path) || !file.exists(truth_path) || !file.exists(manifest_path)) {
    stop("The committed SIM01 bundle or dataset manifest is missing.", call. = FALSE)
  }
  frozen_inputs <- read.csv(manifest_path, stringsAsFactors = FALSE)
  expected_hash <- frozen_inputs$matrix_sha256[frozen_inputs$dataset_key == "simulation"]
  observed_hash <- sha256_object(readRDS(matrix_path))
  if (length(expected_hash) != 1L || !identical(observed_hash, expected_hash)) {
    stop("The committed SIM01 matrix does not match the frozen input manifest.", call. = FALSE)
  }
  validate_simulation_identity(
    "SIM01",
    list(raw_intensity = readRDS(matrix_path), clean_ground_truth = readRDS(truth_path)),
    hash_manifest,
    require_exact = TRUE
  )
  message("SIM01: committed canonical bundle matches the frozen input hash; reusing it.")
  invisible(matrix_path)
}

write_bundle <- function(
  seed_id,
  seed_manifest,
  seed_manifest_path,
  hash_manifest,
  hash_manifest_path,
  output_root,
  force = FALSE
) {
  if (identical(seed_id, "SIM01")) {
    return(verify_committed_primary(output_root, hash_manifest))
  }

  ledger <- simulation_ledger_for(seed_manifest, seed_id)
  generation_rows <- ledger[
    ledger$component %in% .canonical_required_seed_components,
    ,
    drop = FALSE
  ]
  missing_generation <- setdiff(.canonical_required_seed_components, generation_rows$component)
  if (length(missing_generation)) {
    stop(
      seed_id, " lacks generation seed component(s): ",
      paste(missing_generation, collapse = ", "),
      call. = FALSE
    )
  }
  generation_seeds <- setNames(
    generation_rows$seed,
    generation_rows$component
  )[.canonical_required_seed_components]

  hidden_row <- ledger[ledger$component == "hidden_reference", , drop = FALSE]
  if (nrow(hidden_row) != 1L) {
    stop(seed_id, " must have exactly one hidden_reference seed.", call. = FALSE)
  }

  core_path <- file.path(dirname(script_path), "simulation_core.R")
  config_record <- list(
    bundle_format_version = "1.1.0",
    seed_id = seed_id,
    seed_manifest_sha256 = sha256_file(seed_manifest_path),
    simulation_hash_manifest_sha256 = sha256_file(hash_manifest_path),
    simulation_core_sha256 = sha256_file(core_path),
    bundle_writer_sha256 = sha256_file(script_path),
    seed_ledger = ledger[, c("component", "seed"), drop = FALSE],
    save_format = list(matrix = "RDS version 3, gzip", table = "quoted CSV")
  )
  config_hash <- sha256_object(config_record)
  bundle_dir <- file.path(output_root, seed_id)

  if (dir.exists(bundle_dir) && validate_existing_bundle(bundle_dir, config_hash)) {
    message(seed_id, ": existing bundle passed hashes and config validation; reusing it.")
    return(invisible(bundle_dir))
  }
  if (dir.exists(bundle_dir) && !force) {
    stop(
      seed_id, ": an existing bundle failed validation or has a different config. ",
      "Re-run with --force to archive and replace it.",
      call. = FALSE
    )
  }

  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  temporary_dir <- file.path(output_root, paste0(".", seed_id, ".tmp-", Sys.getpid()))
  if (dir.exists(temporary_dir)) {
    stop("Temporary bundle directory already exists: ", temporary_dir, call. = FALSE)
  }
  dir.create(temporary_dir, recursive = FALSE, showWarnings = FALSE)
  committed <- FALSE
  on.exit({
    if (!committed && dir.exists(temporary_dir)) {
      unlink(temporary_dir, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)

  simulation <- generate_canonical_simulation(
    generation_seeds,
    include_artifact_matrices = FALSE
  )
  identity_checks <- validate_simulation_identity(seed_id, simulation, hash_manifest)
  hidden_assignment <- select_canonical_hidden_references(
    simulation$metadata,
    seed = hidden_row$seed,
    hidden_per_plate = 2L
  )

  standardized_metadata <- simulation$metadata
  standardized_metadata$batch <- standardized_metadata$plate
  standardized_metadata$within_batch_order <- standardized_metadata$order_in_plate
  standardized_metadata$class <- ifelse(
    standardized_metadata$sample_type == "control", "QC", "Study"
  )
  standardized_metadata$is_qc <- standardized_metadata$sample_type == "control"
  standardized_metadata$is_hidden_reference <- standardized_metadata$sample_id %in%
    hidden_assignment$sample_id[hidden_assignment$is_hidden_reference]
  standardized_metadata$reference_role <- "not_a_control"
  standardized_metadata$reference_role[standardized_metadata$is_qc] <- "training_control"
  standardized_metadata$reference_role[standardized_metadata$is_hidden_reference] <- "hidden_test"

  saveRDS(
    simulation$raw_intensity,
    file.path(temporary_dir, "raw_intensity.rds"),
    version = 3,
    compress = "gzip"
  )
  saveRDS(
    simulation$clean_ground_truth,
    file.path(temporary_dir, "clean_ground_truth.rds"),
    version = 3,
    compress = "gzip"
  )
  write_csv_stable(standardized_metadata, file.path(temporary_dir, "sample_metadata.csv"))
  write_csv_stable(simulation$feature_truth, file.path(temporary_dir, "feature_truth.csv"))
  write_csv_stable(
    simulation$feature_plate_truth,
    file.path(temporary_dir, "feature_plate_truth.csv")
  )
  write_csv_stable(
    hidden_assignment,
    file.path(temporary_dir, "hidden_reference_assignment.csv")
  )
  write_csv_stable(ledger, file.path(temporary_dir, "seed_ledger.csv"))
  jsonlite::write_json(
    simulation$design_parameters,
    file.path(temporary_dir, "design_parameters.json"),
    auto_unbox = TRUE,
    pretty = TRUE,
    digits = NA
  )

  provenance <- list(
    bundle_format_version = "1.1.0",
    seed_id = seed_id,
    display_realization = isTRUE(ledger$display_realization[[1L]]),
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    config_sha256 = config_hash,
    source_seed_manifest = "analysis/config/simulation_seed_ledger.csv",
    source_seed_manifest_sha256 = sha256_file(seed_manifest_path),
    simulation_hash_manifest = "analysis/config/simulation_bundle_hashes.csv",
    simulation_hash_manifest_sha256 = sha256_file(hash_manifest_path),
    simulation_core = "scripts/robustness/simulation_core.R",
    simulation_core_sha256 = sha256_file(core_path),
    bundle_writer = "scripts/robustness/generate_simulation_bundles.R",
    bundle_writer_sha256 = sha256_file(script_path),
    dimensions = list(
      features = nrow(simulation$raw_intensity),
      injections = ncol(simulation$raw_intensity),
      plates = length(unique(simulation$metadata$plate)),
      controls = sum(simulation$metadata$sample_type == "control"),
      study_injections = sum(simulation$metadata$sample_type == "sample"),
      hidden_references = sum(hidden_assignment$is_hidden_reference)
    ),
    object_sha256 = list(
      raw_intensity = sha256_object(simulation$raw_intensity),
      clean_ground_truth = sha256_object(simulation$clean_ground_truth),
      metadata_core = sha256_object(simulation$metadata),
      feature_truth = sha256_object(simulation$feature_truth),
      feature_plate_truth = sha256_object(simulation$feature_plate_truth),
      hidden_reference_assignment = sha256_object(hidden_assignment)
    ),
    matrix_identity = split(identity_checks, seq_len(nrow(identity_checks))),
    R = list(version = R.version.string, platform = R.version$platform)
  )
  jsonlite::write_json(
    provenance,
    file.path(temporary_dir, "bundle_provenance.json"),
    auto_unbox = TRUE,
    pretty = TRUE,
    digits = NA
  )

  persisted_files <- sort(list.files(temporary_dir, full.names = TRUE, recursive = FALSE))
  file_manifest <- data.frame(
    relative_path = basename(persisted_files),
    bytes = as.numeric(file.info(persisted_files)$size),
    sha256 = vapply(persisted_files, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  write_csv_stable(file_manifest, file.path(temporary_dir, "bundle_files.csv"))

  if (dir.exists(bundle_dir)) {
    archive_dir <- file.path(
      output_root,
      paste0(seed_id, ".replaced-", format(Sys.time(), "%Y%m%dT%H%M%S"))
    )
    if (!file.rename(bundle_dir, archive_dir)) {
      stop("Could not archive existing bundle to: ", archive_dir, call. = FALSE)
    }
    message(seed_id, ": moved the previous invalid bundle to ", archive_dir)
  }
  if (!file.rename(temporary_dir, bundle_dir)) {
    stop("Could not atomically commit bundle to: ", bundle_dir, call. = FALSE)
  }
  committed <- TRUE

  if (!validate_existing_bundle(bundle_dir, config_hash)) {
    stop(seed_id, ": committed bundle failed post-write validation.", call. = FALSE)
  }
  message(
    seed_id, ": wrote canonical bundle (",
    nrow(simulation$raw_intensity), " features x ",
    ncol(simulation$raw_intensity), " injections; ",
    sum(hidden_assignment$is_hidden_reference), " hidden references)."
  )
  invisible(bundle_dir)
}

parse_cli <- function(args) {
  options <- list(seed_ids = character(), all = FALSE, force = FALSE)
  for (arg in args) {
    if (identical(arg, "--all")) {
      options$all <- TRUE
    } else if (identical(arg, "--force")) {
      options$force <- TRUE
    } else if (startsWith(arg, "--seed-id=")) {
      options$seed_ids <- c(options$seed_ids, sub("^--seed-id=", "", arg))
    } else if (grepl("^SIM[0-9]{2}$", arg)) {
      options$seed_ids <- c(options$seed_ids, arg)
    } else if (arg %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript scripts/robustness/generate_simulation_bundles.R ",
        "[--seed-id=SIM01 | SIM01 ... | --all] [--force]\n",
        sep = ""
      )
      quit(save = "no", status = 0L)
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }
  options
}

if (sys.nframe() == 0L) {
  setwd(repo_root)
  cli <- parse_cli(commandArgs(trailingOnly = TRUE))
  seed_manifest_path <- file.path(repo_root, "analysis", "config", "simulation_seed_ledger.csv")
  seed_manifest <- read_seed_manifest(seed_manifest_path)
  hash_manifest_path <- file.path(repo_root, "analysis", "config", "simulation_bundle_hashes.csv")
  hash_manifest <- read_bundle_hash_manifest(hash_manifest_path)
  available_ids <- unique(seed_manifest$seed_id[
    seed_manifest$workstream == "canonical_repeated_simulation"
  ])
  available_ids <- intersect(sprintf("SIM%02d", 1:10), available_ids)
  seed_ids <- if (cli$all) available_ids else unique(cli$seed_ids)
  if (!length(seed_ids)) seed_ids <- "SIM01"
  unknown <- setdiff(seed_ids, available_ids)
  if (length(unknown)) {
    stop("Unknown canonical simulation seed ID(s): ", paste(unknown, collapse = ", "))
  }

  output_root <- file.path(repo_root, "data", "simulated", "canonical")
  for (seed_id in seed_ids) {
    write_bundle(
      seed_id = seed_id,
      seed_manifest = seed_manifest,
      seed_manifest_path = seed_manifest_path,
      hash_manifest = hash_manifest,
      hash_manifest_path = hash_manifest_path,
      output_root = output_root,
      force = cli$force
    )
  }
}
