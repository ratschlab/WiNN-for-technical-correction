assert_frozen_winn <- function(
    release_root = Sys.getenv("WINN_RELEASE_ROOT", unset = ""),
    expected_version = "0.1.4",
    expected_sha256 = "71a0964cee2778b2e5789d20621147e074c7945e813cf76af2ceeb696104aae1") {
  if (!nzchar(release_root)) {
    stop("WINN_RELEASE_ROOT is not set.", call. = FALSE)
  }
  release_root <- normalizePath(release_root, mustWork = TRUE)
  frozen_library <- normalizePath(file.path(release_root, "Rlib"), mustWork = TRUE)
  archive <- file.path(release_root, "package", "winn_0.1.4.tar.gz")
  if (!file.exists(archive)) {
    stop("Frozen WiNN source archive is missing: ", archive, call. = FALSE)
  }
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package digest is required to validate the frozen WiNN archive.", call. = FALSE)
  }
  actual_sha256 <- digest::digest(file = archive, algo = "sha256")
  if (!identical(actual_sha256, expected_sha256)) {
    stop("Frozen WiNN source archive SHA-256 mismatch.", call. = FALSE)
  }
  package_path <- normalizePath(find.package("winn"), mustWork = TRUE)
  if (!startsWith(package_path, paste0(frozen_library, .Platform$file.sep))) {
    stop("WiNN was not loaded from the frozen release library: ", package_path,
         call. = FALSE)
  }
  actual_version <- as.character(utils::packageVersion("winn"))
  if (!identical(actual_version, expected_version)) {
    stop("Expected WiNN ", expected_version, "; loaded ", actual_version, ".",
         call. = FALSE)
  }
  invisible(list(
    version = actual_version,
    package_path = package_path,
    source_archive = normalizePath(archive),
    source_archive_sha256 = actual_sha256
  ))
}
