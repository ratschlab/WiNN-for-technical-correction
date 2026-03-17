if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render the notebooks.")
}

args <- commandArgs(trailingOnly = TRUE)
download_public_data <- "--download-public-data" %in% args

render_dir <- file.path("notebooks", "rendered")
dir.create(render_dir, showWarnings = FALSE, recursive = TRUE)

rmarkdown::render(
  input = file.path("notebooks", "01_simulation_benchmark_data.Rmd"),
  output_dir = render_dir,
  quiet = TRUE,
  envir = new.env(parent = globalenv())
)

rmarkdown::render(
  input = file.path("notebooks", "02_public_mtbls79_preprocessing.Rmd"),
  output_dir = render_dir,
  quiet = TRUE,
  params = list(download_if_missing = download_public_data),
  envir = new.env(parent = globalenv())
)
