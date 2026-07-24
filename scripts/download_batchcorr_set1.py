#!/usr/bin/env python3
"""Download the pinned BatchCorrMetabolomics Set 1 package data and documentation."""

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
SOURCE_DIR = REPO_ROOT / "data" / "public" / "batchcorr_set1" / "source"
COMMIT = "e0c7668140e206dcdae9afa602dd2e1b337ac4f6"
RAW_BASE = f"https://raw.githubusercontent.com/rwehrens/BatchCorrMetabolomics/{COMMIT}"
FILES = {
    "BC.RData": f"{RAW_BASE}/data/BC.RData",
    "BC.Rd": f"{RAW_BASE}/man/BC.Rd",
    "BC_demo.R": f"{RAW_BASE}/demo/BC.R",
    "DESCRIPTION": f"{RAW_BASE}/DESCRIPTION",
    "README.md": f"{RAW_BASE}/README.md",
}
PAPER_XML_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=pmc&id=4796354"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(url: str, target: Path, refresh: bool) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and not refresh:
        return
    request = urllib.request.Request(
        url, headers={"User-Agent": "metabolomics-winn-audit/1.0"}
    )
    fd, temporary_name = tempfile.mkstemp(prefix=target.name + ".", dir=target.parent)
    os.close(fd)
    temporary_path = Path(temporary_name)
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            with temporary_path.open("wb") as handle:
                shutil.copyfileobj(response, handle)
        temporary_path.replace(target)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true")
    args = parser.parse_args()
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    downloaded_at = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
    rows = []

    all_files = dict(FILES)
    all_files["PMC4796354_full_text.xml"] = PAPER_XML_URL
    for filename, url in all_files.items():
        target = SOURCE_DIR / filename
        fetch(url, target, args.refresh)
        rows.append(
            {
                "source_url": url,
                "local_path": str(target.relative_to(REPO_ROOT)),
                "file_size": target.stat().st_size,
                "sha256": sha256(target),
                "download_date_utc": downloaded_at,
                "file_type": mimetypes.guess_type(target.name)[0]
                or "application/octet-stream",
                "inferred_assay": "BatchCorrMetabolomics Set 1 LC-MS negative ion",
                "used": "true",
                "source_commit": COMMIT if url.startswith(RAW_BASE) else "",
            }
        )

    fields = list(rows[0])
    with (SOURCE_DIR / "download_manifest.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Downloaded/reused {len(rows)} small package and article files.")
    print(f"Pinned Git commit: {COMMIT}")


if __name__ == "__main__":
    main()
