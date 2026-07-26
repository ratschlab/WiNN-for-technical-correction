#!/usr/bin/env Rscript

release_root <- Sys.getenv("WINN_RELEASE_ROOT", unset = "")
if (!nzchar(release_root) || !dir.exists(release_root)) {
  stop("WINN_RELEASE_ROOT must identify the release directory.", call. = FALSE)
}
source(file.path(release_root, "analysis", "release_helpers.R"), local = FALSE)

audit_rows <- list()
add_audit <- function(family, dataset_key, task_id, method_id, check, passed, detail) {
  audit_rows[[length(audit_rows) + 1L]] <<- data.frame(
    family = family,
    dataset_key = dataset_key,
    task_id = task_id,
    method_id = method_id,
    check = check,
    passed = isTRUE(passed),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

as_flag <- function(value) {
  if (is.logical(value)) return(isTRUE(value))
  identical(tolower(as.character(value)), "true")
}

control_ids <- function(details, metadata) {
  controls <- details$control_samples
  if (is.null(controls) || !length(controls)) return(character())
  controls <- as.integer(controls)
  if (anyNA(controls) || any(controls < 1L | controls > nrow(metadata))) return(NA_character_)
  as.character(metadata$sample_id[controls])
}

audit_winn_controls <- function(
  family, dataset_key, task_id, method_id, details_path, metadata,
  training_ids, hidden_ids
) {
  details <- readRDS(details_path)
  observed <- control_ids(details, metadata)
  if (identical(method_id, "winn_fixed_no_qc")) {
    passed <- length(observed) == 0L
    detail <- paste0("fixed mode control count=", length(observed))
  } else {
    passed <- !anyNA(observed) && setequal(observed, training_ids) &&
      !length(intersect(observed, hidden_ids))
    detail <- paste0(
      "recorded controls=", length(observed),
      "; training references=", length(training_ids),
      "; withheld overlap=", length(intersect(observed, hidden_ids))
    )
  }
  add_audit(
    family, dataset_key, task_id, method_id,
    "recorded_winn_controls_match_training_role", passed, detail
  )
}

selection <- read.csv(
  file.path(release_root, "analysis", "config", "endpoint_free_selection_manifest.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
forbidden_fields <- c(
  "hidden_qc_values_available", "biological_labels_available",
  "replicate_identities_available", "preservation_endpoints_available",
  "batch_r2_available", "final_rank_available"
)
selection_ok <- all(forbidden_fields %in% names(selection)) &&
  all(vapply(forbidden_fields, function(field) {
    !any(vapply(selection[[field]], as_flag, logical(1)))
  }, logical(1)))
add_audit(
  "selection", "all", "selection_manifest", "all",
  "evaluation_fields_unavailable_during_selection", selection_ok,
  paste(length(forbidden_fields), "forbidden availability fields audited")
)

runner_text <- paste(deparse(body(release_run_method)), collapse = "\n")
winn_runner_text <- paste(deparse(body(release_run_winn)), collapse = "\n")
winn_runner_compact <- gsub("[[:space:]]+", " ", winn_runner_text)
runner_role_ok <-
  grepl("prepared$training_ids", runner_text, fixed = TRUE) &&
  !grepl("prepared$hidden_ids", runner_text, fixed = TRUE) &&
  !grepl("prepared$external_ids", runner_text, fixed = TRUE) &&
  grepl("auto_batch = TRUE", runner_text, fixed = TRUE) &&
  grepl(
    "batch = if (isTRUE(auto_batch)) NULL else prepared$meta$batch",
    winn_runner_compact, fixed = TRUE
  )
add_audit(
  "implementation", "all", "release_run_method", "all",
  "method_runner_exposes_training_roles_only_and_omits_batch_in_auto_batch_mode",
  runner_role_ok,
  "The frozen runner passes training reference IDs only; hidden/external IDs are absent, and auto-batch WiNN receives batch=NULL."
)

dataset_keys <- c(
  "simulation", "mtbls79", "clinical_fiams", "batchcorr_set1",
  "sacurine", "waveica"
)
winn_modes <- c("winn_auto_qc", "winn_auto_batch_qc", "winn_fixed_no_qc")
all_methods <- c(
  "raw", "combat", "qc_rlsc", "qc_rfsc", "tiger", "serrf",
  winn_modes
)
training_reference_methods <- c(
  "qc_rlsc", "qc_rfsc", "tiger", "serrf",
  "winn_auto_qc", "winn_auto_batch_qc"
)

for (dataset_key in dataset_keys) {
  prepared <- release_load_dataset(dataset_key)
  roles_ok <- !length(intersect(
    prepared$training_ids, c(prepared$hidden_ids, prepared$external_ids)
  )) && !length(intersect(prepared$hidden_ids, prepared$external_ids)) &&
    all(c(
      prepared$training_ids, prepared$hidden_ids, prepared$external_ids
    ) %in% prepared$meta$sample_id)
  add_audit(
    "primary", dataset_key, dataset_key, "all",
    "training_and_hidden_reference_roles_are_disjoint", roles_ok,
    paste0(
      "training=", length(prepared$training_ids),
      "; hidden=", length(prepared$hidden_ids)
    )
  )
  for (method_id in all_methods) {
    configuration <- selection[
      selection$dataset_key == dataset_key & selection$method_id == method_id,
      , drop = FALSE
    ]
    allowed_ids <- if (method_id %in% training_reference_methods) {
      prepared$training_ids
    } else character()
    protocol_ok <-
      nrow(configuration) == 1L &&
      !length(intersect(allowed_ids, c(prepared$hidden_ids, prepared$external_ids))) &&
      !as_flag(configuration$hidden_qc_values_available) &&
      !as_flag(configuration$biological_labels_available) &&
      !as_flag(configuration$replicate_identities_available)
    add_audit(
      "primary", dataset_key, dataset_key, method_id,
      "frozen_primary_control_protocol_excludes_hidden_and_external_references",
      protocol_ok,
      paste0(
        "allowed training references=", length(allowed_ids),
        "; hidden/external overlap=",
        length(intersect(allowed_ids, c(prepared$hidden_ids, prepared$external_ids)))
      )
    )
  }
  for (method_id in winn_modes) {
    audit_winn_controls(
      "primary", dataset_key, dataset_key, method_id,
      file.path(
        release_root, "results", "primary", dataset_key, method_id,
        "method_details.rds"
      ),
      prepared$meta, prepared$training_ids,
      c(prepared$hidden_ids, prepared$external_ids)
    )
  }
}

for (seed_id in sprintf("SIM%02d", 1:10)) {
  prepared <- release_load_simulation(seed_id)
  roles_ok <- !length(intersect(prepared$training_ids, prepared$hidden_ids)) &&
    all(c(prepared$training_ids, prepared$hidden_ids) %in% prepared$meta$sample_id)
  add_audit(
    "simulation_seed", "simulation", seed_id, "all",
    "training_and_hidden_reference_roles_are_disjoint", roles_ok,
    paste0(
      "training=", length(prepared$training_ids),
      "; hidden=", length(prepared$hidden_ids)
    )
  )
  for (method_id in winn_modes) {
    audit_winn_controls(
      "simulation_seed", "simulation", seed_id, method_id,
      file.path(
        release_root, "results", "simulation_seeds", seed_id, method_id,
        "method_details.rds"
      ),
      prepared$meta, prepared$training_ids, prepared$hidden_ids
    )
  }
}

reference_tasks <- read.csv(
  file.path(release_root, "analysis", "config", "reference_split_task_manifest.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
metadata_cache <- new.env(parent = emptyenv())
for (index in seq_len(nrow(reference_tasks))) {
  task <- reference_tasks[index, , drop = FALSE]
  dataset_key <- as.character(task$dataset_key)
  split_id <- as.character(task$split_id)
  method_id <- as.character(task$method_id)
  directory <- file.path(release_root, as.character(task$output_dir))
  assignment <- read.csv(
    file.path(directory, "reference_assignment.csv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  training_ids <- as.character(assignment$sample_id[assignment$assignment == "training"])
  hidden_ids <- as.character(assignment$sample_id[assignment$assignment == "hidden"])
  external_ids <- as.character(assignment$sample_id[assignment$assignment == "external"])
  role_ok <- !anyDuplicated(assignment$sample_id) &&
    !length(intersect(training_ids, c(hidden_ids, external_ids))) &&
    !length(intersect(hidden_ids, external_ids))

  manifest <- jsonlite::read_json(
    file.path(directory, "run_manifest.json"), simplifyVector = TRUE
  )
  manifest_ok <- identical(
    as.character(manifest$training_qc_ids_sha256),
    release_sha256_object(sort(training_ids))
  ) && identical(
    as.character(manifest$hidden_qc_ids_sha256),
    release_sha256_object(sort(hidden_ids))
  ) && identical(as.integer(manifest$n_training_qc), length(training_ids)) &&
    identical(as.integer(manifest$n_hidden_qc), length(hidden_ids)) &&
    identical(as.integer(manifest$n_external_reference), length(external_ids)) &&
    !as_flag(manifest$hidden_qc_values_available_during_selection) &&
    !as_flag(manifest$biological_labels_available_during_selection) &&
    !as_flag(manifest$replicate_identities_available_during_selection)

  detail <- paste0(
    "training=", length(training_ids), "; hidden=", length(hidden_ids),
    "; external=", length(external_ids), "; role_disjoint=", role_ok,
    "; manifest_blinding=", manifest_ok
  )
  add_audit(
    "reference_split", dataset_key, as.character(task$task_id), method_id,
    "assignment_roles_and_blinding_manifest_match", role_ok && manifest_ok, detail
  )

  if (method_id %in% winn_modes) {
    if (!exists(dataset_key, envir = metadata_cache, inherits = FALSE)) {
      assign(dataset_key, release_load_dataset(dataset_key)$meta, envir = metadata_cache)
    }
    audit_winn_controls(
      "reference_split", dataset_key, as.character(task$task_id), method_id,
      file.path(directory, "selected_details.rds"),
      get(dataset_key, envir = metadata_cache, inherits = FALSE),
      training_ids, c(hidden_ids, external_ids)
    )
  }
}

audit <- do.call(rbind, audit_rows)
output_dir <- file.path(release_root, "results", "final", "validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  audit,
  file.path(output_dir, "reference_separation_audit.csv"),
  row.names = FALSE, quote = TRUE
)
summary <- aggregate(
  passed ~ family + check, data = audit,
  FUN = function(x) c(passed = sum(x), total = length(x))
)
summary <- data.frame(
  family = summary$family,
  check = summary$check,
  passed = summary$passed[, "passed"],
  total = summary$passed[, "total"],
  stringsAsFactors = FALSE
)
write.csv(
  summary,
  file.path(output_dir, "reference_separation_summary.csv"),
  row.names = FALSE, quote = TRUE
)
if (!all(audit$passed)) {
  stop(
    "Reference-separation audit failed for ",
    paste(audit$task_id[!audit$passed], collapse = ", "),
    call. = FALSE
  )
}
message("Reference-separation audit passed (", nrow(audit), " checks).")
