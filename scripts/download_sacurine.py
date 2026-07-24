#!/usr/bin/env python3
"""Download the pinned, small Sacurine/W4M00001 benchmark artifacts."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import mimetypes
import os
from pathlib import Path
import shutil
import tempfile
import urllib.request


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIR = REPO_ROOT / "data" / "public" / "sacurine" / "source"
MANIFEST_DIR = REPO_ROOT / "data" / "public" / "sacurine" / "manifests"
SOURCE_REPOSITORY = "https://github.com/SciDoPhenIA/phenomis"
SOURCE_COMMIT = "1e83ce6997a8d16b89ce5f0f899a1570004ebc0e"
RAW_BASE = f"https://raw.githubusercontent.com/SciDoPhenIA/phenomis/{SOURCE_COMMIT}"

FILES = {
    "Galaxy1_dataMatrix.tabular": {
        "repository_path": "inst/extdata/W4M00001_Sacurine-statistics/Galaxy1_dataMatrix.tabular",
        "git_blob_sha1": "49354577dbabd424bd8544cce5373b8a6501aeb9",
        "role": "uncorrected_or_minimally_processed_feature_intensity_matrix",
        "used": True,
    },
    "Galaxy2_sampleMetadata.tabular": {
        "repository_path": "inst/extdata/W4M00001_Sacurine-statistics/Galaxy2_sampleMetadata.tabular",
        "git_blob_sha1": "4834c36cb5f1a8fcf7a4103dec735c6478f23ed2",
        "role": "injection_metadata",
        "used": True,
    },
    "Galaxy3_variableMetadata.tabular": {
        "repository_path": "inst/extdata/W4M00001_Sacurine-statistics/Galaxy3_variableMetadata.tabular",
        "git_blob_sha1": "c3f77b61351009ffa254de6edcbe02e7ca4178df",
        "role": "feature_annotation",
        "used": True,
    },
    "phenomis_DESCRIPTION": {
        "repository_path": "DESCRIPTION",
        "git_blob_sha1": "423935909669896e79b442d841f3b37171509a32",
        "role": "package_provenance_and_license",
        "used": True,
    },
    "phenomis_vignette.Rmd": {
        "repository_path": "vignettes/phenomis.Rmd",
        "git_blob_sha1": "9ad87c38ffbe4883e5b0f67df3a6f2820b666157",
        "role": "processing_provenance",
        "used": True,
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_blob_sha1(path: Path) -> str:
    payload = path.read_bytes()
    digest = hashlib.sha1()
    digest.update(f"blob {len(payload)}\0".encode("ascii"))
    digest.update(payload)
    return digest.hexdigest()


def fetch(url: str, target: Path, expected_blob: str, refresh: bool) -> str:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and not refresh and git_blob_sha1(target) == expected_blob:
        return "reused_checksum_match"

    request = urllib.request.Request(
        url, headers={"User-Agent": "metabolomics-winn-sacurine-audit/1.0"}
    )
    fd, temporary_name = tempfile.mkstemp(prefix=target.name + ".", dir=target.parent)
    os.close(fd)
    temporary_path = Path(temporary_name)
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            with temporary_path.open("wb") as handle:
                shutil.copyfileobj(response, handle)
        observed_blob = git_blob_sha1(temporary_path)
        if observed_blob != expected_blob:
            raise RuntimeError(
                f"Git blob checksum mismatch for {target.name}: "
                f"expected {expected_blob}, observed {observed_blob}"
            )
        temporary_path.replace(target)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()
    return "downloaded"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true")
    args = parser.parse_args()

    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    downloaded_at = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
    rows = []

    for local_name, metadata in FILES.items():
        repository_path = metadata["repository_path"]
        url = f"{RAW_BASE}/{repository_path}"
        target = SOURCE_DIR / local_name
        disposition = fetch(
            url, target, str(metadata["git_blob_sha1"]), args.refresh
        )
        rows.append(
            {
                "source_url": url,
                "source_repository": SOURCE_REPOSITORY,
                "source_commit": SOURCE_COMMIT,
                "source_repository_path": repository_path,
                "source_git_blob_sha1": metadata["git_blob_sha1"],
                "local_path": str(target.relative_to(REPO_ROOT)),
                "file_size": target.stat().st_size,
                "sha256": sha256(target),
                "download_date_utc": downloaded_at,
                "file_type": mimetypes.guess_type(target.name)[0] or "text/tab-separated-values",
                "license": "Artistic-2.0 (phenomis DESCRIPTION)",
                "dataset": "Sacurine / W4M00001 / MTBLS404 negative-ion LC-MS",
                "role": metadata["role"],
                "used": str(bool(metadata["used"])).lower(),
                "disposition": disposition,
            }
        )

    manifest_path = MANIFEST_DIR / "download_manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Downloaded/reused {len(rows)} small Sacurine files.")
    print(f"Pinned source commit: {SOURCE_COMMIT}")
    print(f"Manifest: {manifest_path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
