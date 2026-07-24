#!/usr/bin/env python3
"""Download the small, uncorrected MTBLS79 processed workbook used here."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import os
from pathlib import Path
import shutil
import tempfile
import urllib.request


REPO_ROOT = Path(__file__).resolve().parents[1]
TARGET = REPO_ROOT / "data" / "public" / "raw" / "Dataset07__SFPM.xlsx"
MANIFEST = REPO_ROOT / "data" / "public" / "raw" / "download_manifest.csv"
URL = (
    "https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/"
    "MTBLS79/FILES/Dataset07__SFPM.xlsx"
)
EXPECTED_SHA256 = "ad27f28fa57dc334a7c2ae3e9c1173263969d1c7a7062719ddeb250e282175e7"
EXPECTED_SIZE = 4_060_472


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def valid(path: Path) -> bool:
    return (
        path.is_file()
        and path.stat().st_size == EXPECTED_SIZE
        and sha256(path) == EXPECTED_SHA256
    )


def download(refresh: bool) -> str:
    TARGET.parent.mkdir(parents=True, exist_ok=True)
    if valid(TARGET) and not refresh:
        return "reused_checksum_match"

    request = urllib.request.Request(URL, headers={"User-Agent": "winn-benchmark/1.0"})
    descriptor, temporary_name = tempfile.mkstemp(prefix=TARGET.name + ".", dir=TARGET.parent)
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with urllib.request.urlopen(request, timeout=180) as response:
            with temporary.open("wb") as handle:
                shutil.copyfileobj(response, handle)
        if not valid(temporary):
            raise RuntimeError("Downloaded MTBLS79 workbook failed size or SHA-256 validation.")
        temporary.replace(TARGET)
    finally:
        if temporary.exists():
            temporary.unlink()
    return "downloaded"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true")
    args = parser.parse_args()
    disposition = download(args.refresh)
    row = {
        "source_url": URL,
        "local_path": str(TARGET.relative_to(REPO_ROOT)),
        "file_size": TARGET.stat().st_size,
        "sha256": sha256(TARGET),
        "download_date_utc": dt.datetime.now(dt.timezone.utc)
        .replace(microsecond=0)
        .isoformat(),
        "file_type": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        "dataset": "MTBLS79",
        "role": "uncorrected_or_minimally_processed_feature_intensity_workbook",
        "used": "true",
        "disposition": disposition,
    }
    with MANIFEST.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)
    print(f"{disposition}: {TARGET.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
