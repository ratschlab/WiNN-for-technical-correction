#!/usr/bin/env python3
"""Download the pinned WaveICA adenocarcinoma matrix and small supporting source files."""

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
SOURCE_DIR = REPO_ROOT / "data" / "public" / "waveica_adenocarcinoma" / "source"
MANIFEST_DIR = REPO_ROOT / "data" / "public" / "waveica_adenocarcinoma" / "manifests"
SOURCE_REPOSITORY = "https://github.com/dengkuistat/WaveICA_2.0"
SOURCE_COMMIT = "56fff2e6b0410b5957c6ea83bc658241df222f41"
RAW_BASE = f"https://raw.githubusercontent.com/dengkuistat/WaveICA_2.0/{SOURCE_COMMIT}"

FILES = {
    "Amide_data.rda": ("data/Amide_data.rda", "91bf470baeaf8359303ac2bac6cb7ec246e97520", "feature_intensity_matrix_and_injection_metadata", True),
    "Amide_data.Rd": ("man/Amide_data.Rd", "07331b0d05573bfc6ac3b1854d0a8b5e95fc9e8f", "dataset_documentation", True),
    "WaveICA_2.0.R": ("R/WaveICA_2.0.R", "bdd8b4c056bd79de1e1de110090511478fafbe27", "published_method_and_input_provenance", True),
    "unbiased_stICA.R": ("R/unbiased_stICA.R", "f3216c3029147bc51521081a3fb528b90fddee8f", "supporting_package_source", False),
    "README.md": ("README.md", "3e59c5f021401e4e04342e3c28629890ec6544b4", "repository_documentation", True),
    "DESCRIPTION": ("DESCRIPTION", "39781ad4c7a7fc7168124375ddb84e7458d30797", "package_provenance", True),
    "LICENSE": ("LICENSE", "b1fd004ebf1aae5f1d7e2ca9c2200671ada79b0b", "repository_license", True),
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
        url, headers={"User-Agent": "metabolomics-winn-waveica-audit/1.0"}
    )
    fd, temporary_name = tempfile.mkstemp(prefix=target.name + ".", dir=target.parent)
    os.close(fd)
    temporary_path = Path(temporary_name)
    try:
        with urllib.request.urlopen(request, timeout=300) as response:
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

    for local_name, (repository_path, expected_blob, role, used) in FILES.items():
        url = f"{RAW_BASE}/{repository_path}"
        target = SOURCE_DIR / local_name
        disposition = fetch(url, target, expected_blob, args.refresh)
        rows.append(
            {
                "source_url": url,
                "source_repository": SOURCE_REPOSITORY,
                "source_commit": SOURCE_COMMIT,
                "source_repository_path": repository_path,
                "source_git_blob_sha1": expected_blob,
                "local_path": str(target.relative_to(REPO_ROOT)),
                "file_size": target.stat().st_size,
                "sha256": sha256(target),
                "download_date_utc": downloaded_at,
                "file_type": mimetypes.guess_type(target.name)[0] or "application/octet-stream",
                "license": "MIT (repository LICENSE); package DESCRIPTION contains the conflicting value 'No'",
                "dataset": "WaveICA 2.0 human plasma adenocarcinoma demonstration dataset",
                "role": role,
                "used": str(bool(used)).lower(),
                "disposition": disposition,
            }
        )

    manifest_path = MANIFEST_DIR / "download_manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Downloaded/reused {len(rows)} WaveICA matrix/source files.")
    print(f"Pinned source commit: {SOURCE_COMMIT}")
    print(f"Manifest: {manifest_path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
