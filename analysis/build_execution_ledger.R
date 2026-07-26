#!/usr/bin/env Rscript

release_root <- Sys.getenv("WINN_RELEASE_ROOT", unset = "")
if (!nzchar(release_root) || !dir.exists(release_root)) {
  stop("WINN_RELEASE_ROOT must identify the release directory.", call. = FALSE)
}

status_dir <- file.path(release_root, "scheduler", "status")
status_paths <- list.files(status_dir, pattern = "[.]tsv$", full.names = TRUE)
if (!length(status_paths)) stop("No scheduler status files were found.", call. = FALSE)

events <- do.call(rbind, lapply(status_paths, function(path) {
  value <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(value)) return(NULL)
  names(value)[names(value) == "started"] <- "event_time"
  value$status_file <- basename(path)
  value
}))
events$event_time <- as.character(events$event_time)
events$event_order <- seq_len(nrow(events))
events <- events[order(events$event_time, events$event_order), ]

attempts <- do.call(rbind, lapply(split(events, events$status_file), function(value) {
  value <- value[order(value$event_time, value$event_order), , drop = FALSE]
  final_status <- as.character(tail(value$status, 1L))
  data.frame(
    task_id = as.character(value$task_id[[1L]]),
    job_id = as.character(value$job_id[[1L]]),
    array_index = as.integer(value$array_index[[1L]]),
    started_at = as.character(value$event_time[[1L]]),
    finished_at = if (final_status == "running") NA_character_ else as.character(tail(value$event_time, 1L)),
    status = final_status,
    status_file = as.character(value$status_file[[1L]]),
    stringsAsFactors = FALSE
  )
}))
attempts$status_source <- "runner event"
open_attempts <- which(attempts$status == "running")
if (length(open_attempts) && nzchar(Sys.which("sacct"))) {
  for (index in open_attempts) {
    job_spec <- paste0(attempts$job_id[index], "_", attempts$array_index[index])
    accounting <- tryCatch(
      system2(
        "sacct",
        c("-X", "-j", job_spec, "--format=State,End,ExitCode", "-n", "-P"),
        stdout = TRUE, stderr = FALSE
      ),
      error = function(e) character()
    )
    accounting <- accounting[nzchar(trimws(accounting))]
    if (!length(accounting)) next
    fields <- strsplit(accounting[[1L]], "|", fixed = TRUE)[[1L]]
    if (length(fields) < 3L) next
    state <- toupper(sub("[ +].*$", "", trimws(fields[[1L]])))
    if (state %in% c("RUNNING", "PENDING", "CONFIGURING", "COMPLETING")) next
    exit_code <- trimws(fields[[3L]])
    attempts$status[index] <- switch(
      state,
      COMPLETED = "completed",
      CANCELLED = "canceled",
      TIMEOUT = "failed_timeout",
      OUT_OF_MEMORY = "failed_out_of_memory",
      paste0("failed_", tolower(state), "_", exit_code)
    )
    end_time <- trimws(fields[[2L]])
    attempts$finished_at[index] <- if (nzchar(end_time) && end_time != "Unknown") {
      end_time
    } else {
      NA_character_
    }
    attempts$status_source[index] <- "Slurm accounting reconciliation"
  }
}
attempts <- attempts[order(attempts$task_id, attempts$started_at, attempts$job_id), , drop = FALSE]
attempts$attempt_index <- ave(
  seq_len(nrow(attempts)), attempts$task_id,
  FUN = function(index) seq_along(index)
)
attempts$total_attempts <- ave(
  seq_len(nrow(attempts)), attempts$task_id,
  FUN = length
)
attempts$retry_count <- attempts$total_attempts - 1L

read_config <- function(name) {
  read.csv(
    file.path(release_root, "analysis", "config", name),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}
primary <- read_config("primary_run_matrix.csv")
reference <- read_config("reference_split_task_manifest.csv")
simulation <- read_config("simulation_task_manifest.csv")
partial <- read_config("partial_confounding_task_manifest.csv")
simulation_configuration <- if ("public_reproduction_decision" %in% names(simulation)) {
  paste0(
    "original_release=", simulation$decision,
    "; public_reproduction=", simulation$public_reproduction_decision
  )
} else {
  as.character(simulation$decision)
}

lookup <- rbind(
  data.frame(
    task_id = primary$task_id, output_path = primary$output_dir,
    dataset_key = primary$dataset_key, method_id = primary$method_id,
    seed_id = as.character(primary$seed), split_id = NA_character_,
    configuration = "frozen primary configuration", family = "primary",
    stringsAsFactors = FALSE
  ),
  data.frame(
    task_id = reference$task_id, output_path = reference$output_dir,
    dataset_key = reference$dataset_key, method_id = reference$method_id,
    seed_id = NA_character_, split_id = reference$split_id,
    configuration = reference$selection_rule, family = "reference_split",
    stringsAsFactors = FALSE
  ),
  data.frame(
    task_id = simulation$task_id, output_path = simulation$output_dir,
    dataset_key = "simulation", method_id = simulation$method_id,
    seed_id = simulation$seed_id, split_id = NA_character_,
    configuration = simulation_configuration, family = "simulation_method",
    stringsAsFactors = FALSE
  ),
  data.frame(
    task_id = partial$task_id,
    output_path = file.path("results", "partial_confounding", partial$scenario_id),
    dataset_key = "partial_confounding", method_id = partial$method_id,
    seed_id = partial$seed_id, split_id = partial$scenario_id,
    configuration = paste(partial$confounding_type, partial$level, sep = "="),
    family = "partial_confounding", stringsAsFactors = FALSE
  )
)

dataset_keys <- unique(primary$dataset_key)
simulation_seeds <- unique(simulation$seed_id)
reference_evaluations <- unique(reference[c("dataset_key", "split_id")])
lookup <- rbind(
  lookup,
  data.frame(
    task_id = paste0("ABLATION_", dataset_keys),
    output_path = file.path("results", "ablations", dataset_keys),
    dataset_key = dataset_keys, method_id = NA_character_, seed_id = NA_character_,
    split_id = NA_character_, configuration = "nine-variant WiNN ablation bundle",
    family = "primary_ablation", stringsAsFactors = FALSE
  ),
  data.frame(
    task_id = paste0("SIMABL_", simulation_seeds),
    output_path = file.path("results", "simulation_seed_ablations", simulation_seeds),
    dataset_key = "simulation", method_id = NA_character_, seed_id = simulation_seeds,
    split_id = NA_character_, configuration = "simulation WiNN ablation bundle",
    family = "simulation_ablation", stringsAsFactors = FALSE
  ),
  data.frame(
    task_id = paste0("EVAL_primary_", dataset_keys),
    output_path = file.path("results", "evaluation", "primary", dataset_keys),
    dataset_key = dataset_keys, method_id = NA_character_, seed_id = NA_character_,
    split_id = NA_character_, configuration = "frozen primary evaluation",
    family = "primary_evaluation", stringsAsFactors = FALSE
  ),
  data.frame(
    task_id = paste0("EVAL_ablations_", dataset_keys),
    output_path = file.path("results", "evaluation", "ablations", dataset_keys),
    dataset_key = dataset_keys, method_id = NA_character_, seed_id = NA_character_,
    split_id = NA_character_, configuration = "frozen ablation evaluation",
    family = "ablations_evaluation", stringsAsFactors = FALSE
  ),
  data.frame(
    task_id = paste0("SIMEVAL_", simulation_seeds),
    output_path = file.path("results", "evaluation", "simulation_seeds", simulation_seeds),
    dataset_key = "simulation", method_id = NA_character_, seed_id = simulation_seeds,
    split_id = NA_character_, configuration = "frozen repeated-simulation evaluation",
    family = "simulation_evaluation", stringsAsFactors = FALSE
  ),
  data.frame(
    task_id = paste0("REFEVAL_", reference_evaluations$dataset_key, "_", reference_evaluations$split_id),
    output_path = file.path(
      "results", "evaluation", "reference_splits",
      reference_evaluations$dataset_key, reference_evaluations$split_id
    ),
    dataset_key = reference_evaluations$dataset_key, method_id = NA_character_,
    seed_id = NA_character_, split_id = reference_evaluations$split_id,
    configuration = "frozen shared reference-split evaluation",
    family = "reference_evaluation", stringsAsFactors = FALSE
  )
)

extra_ids <- setdiff(unique(attempts$task_id), lookup$task_id)
extra_lookup <- lapply(extra_ids, function(task_id) {
  if (grepl("^ABLATION_", task_id)) {
    dataset <- sub("^ABLATION_", "", task_id)
    values <- c(file.path("results", "ablations", dataset), NA, "primary_ablation")
  } else if (grepl("^SIMABL_", task_id)) {
    seed_id <- sub("^SIMABL_", "", task_id)
    values <- c(file.path("results", "simulation_seed_ablations", seed_id), NA, "simulation_ablation")
  } else if (grepl("^SIMEVAL_", task_id)) {
    seed_id <- sub("^SIMEVAL_", "", task_id)
    values <- c(file.path("results", "evaluation", "simulation_seeds", seed_id), NA, "simulation_evaluation")
  } else if (grepl("^REFEVAL_", task_id)) {
    key <- sub("^REFEVAL_", "", task_id)
    split_id <- sub("^.*_(SPLIT[0-9]+)$", "\\1", key)
    dataset <- sub(paste0("_", split_id, "$"), "", key)
    values <- c(file.path("results", "evaluation", "reference_splits", dataset, split_id), NA, "reference_evaluation")
  } else if (grepl("^EVAL_", task_id)) {
    key <- sub("^EVAL_", "", task_id)
    family <- sub("_.*$", "", key)
    dataset <- sub(paste0("^", family, "_"), "", key)
    values <- c(file.path("results", "evaluation", family, dataset), NA, paste0(family, "_evaluation"))
  } else {
    values <- c(NA, NA, "unclassified")
  }
  data.frame(
    task_id = task_id, output_path = values[[1L]], dataset_key = NA_character_,
    method_id = values[[2L]], seed_id = NA_character_, split_id = NA_character_,
    configuration = "unclassified scheduler record",
    family = values[[3L]], stringsAsFactors = FALSE
  )
})
if (length(extra_lookup)) lookup <- rbind(lookup, do.call(rbind, extra_lookup))
lookup <- lookup[!duplicated(lookup$task_id), , drop = FALSE]

ledger <- merge(attempts, lookup, by = "task_id", all.x = TRUE, sort = FALSE)
ledger$requested_cpus <- 1L
ledger$requested_memory <- ""
ledger$requested_time <- ""
fast_primary <- ledger$family %in% c("primary", "simulation_method") &
  !is.na(ledger$method_id) & ledger$method_id != "tiger"
ledger$requested_memory[fast_primary] <- "24G"
ledger$requested_time[fast_primary] <- "12:00:00"
tiger <- ledger$method_id == "tiger" & !is.na(ledger$method_id)
ledger$requested_cpus[tiger] <- 8L
ledger$requested_memory[tiger] <- "8G per CPU"
ledger$requested_time[tiger] <- "48:00:00"
fast_reference <- ledger$family == "reference_split" & !tiger
ledger$requested_memory[fast_reference] <- "32G"
ledger$requested_time[fast_reference] <- "36:00:00"
ledger$requested_memory[ledger$family %in% c("primary_ablation", "partial_confounding")] <- "32G"
ledger$requested_time[ledger$family %in% c("primary_ablation", "partial_confounding")] <- "24:00:00"
ledger$requested_memory[ledger$family == "simulation_ablation"] <- "24G"
ledger$requested_time[ledger$family == "simulation_ablation"] <- "24:00:00"
evaluation <- grepl("evaluation$", ledger$family)
ledger$requested_memory[evaluation] <- "48G"
ledger$requested_time[evaluation] <- "36:00:00"

completion_paths <- file.path(release_root, ledger$output_path, "complete.json")
completion_present <- !is.na(ledger$output_path) & file.exists(completion_paths)
ledger$validation_status <- ifelse(
  ledger$status == "completed" & completion_present,
  "complete_marker_present",
  ifelse(ledger$status == "completed", "output_validation_pending", "not_complete")
)
ledger <- ledger[order(ledger$task_id, ledger$attempt_index), , drop = FALSE]

output_dir <- file.path(release_root, "results", "final", "validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(events, file.path(output_dir, "scheduler_event_history.csv"), row.names = FALSE)
write.csv(ledger, file.path(output_dir, "execution_ledger.csv"), row.names = FALSE)

latest_attempt <- ledger[order(ledger$task_id, ledger$attempt_index), , drop = FALSE]
latest_attempt <- latest_attempt[!duplicated(latest_attempt$task_id, fromLast = TRUE), , drop = FALSE]
latest_fields <- latest_attempt[c(
  "task_id", "job_id", "array_index", "started_at", "finished_at", "status",
  "total_attempts", "retry_count", "requested_cpus", "requested_memory",
  "requested_time"
)]
names(latest_fields)[names(latest_fields) == "status"] <- "latest_scheduler_status"
logical_tasks <- merge(lookup, latest_fields, by = "task_id", all.x = TRUE, sort = FALSE)
logical_tasks <- logical_tasks[!duplicated(logical_tasks$task_id), , drop = FALSE]
logical_completion_paths <- file.path(release_root, logical_tasks$output_path, "complete.json")
logical_tasks$complete_marker_present <-
  !is.na(logical_tasks$output_path) & file.exists(logical_completion_paths)
logical_tasks$execution_provenance <- ifelse(
  !is.na(logical_tasks$job_id), "scheduler attempt recorded",
  ifelse(
    logical_tasks$complete_marker_present,
    "validated existing artifact; no scheduler attempt recorded in this release ledger",
    "no completed artifact or scheduler attempt recorded"
  )
)
logical_tasks$validation_status <- ifelse(
  logical_tasks$complete_marker_present, "complete_marker_present", "missing_complete_marker"
)
logical_tasks <- logical_tasks[order(logical_tasks$family, logical_tasks$task_id), , drop = FALSE]
write.csv(logical_tasks, file.path(output_dir, "logical_task_ledger.csv"), row.names = FALSE)
message(
  "Execution ledger written with ", nrow(ledger), " scheduler task-attempt rows and ",
  nrow(logical_tasks), " prespecified logical tasks."
)
