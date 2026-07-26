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
from typing import Optional
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
EXPECTED_SHA256 = {
    "BC.RData": "a5d3918c69a902af3886c0141292a614cc20c01526a68930853157e5c9aec113",
    "BC.Rd": "c1d60f1765436187bffb23a86ce4254d0edbe0d967917cba478249f1ddcaf644",
    "BC_demo.R": "70f0abb3d056eb12ea283db367c6b05ee7e6b862c8c566c531a824636061e187",
    "DESCRIPTION": "3aee180104116f2d61226b5aed8a232616aae2504aebd501ae32887d462ba830",
    "README.md": "8c41f92e73731a7caecc088fbbeadca7783502d0bc9deca7465b44ef7e8e7d95",
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


def fetch(url: str, target: Path, refresh: bool, expected_sha256: Optional[str]) -> str:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and not refresh:
        if expected_sha256 is None or sha256(target) == expected_sha256:
            return "reused_checksum_match" if expected_sha256 else "reused_recorded_file"
        raise RuntimeError(f"Existing file failed SHA-256 validation: {target}")
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
        if expected_sha256 is not None and sha256(temporary_path) != expected_sha256:
            raise RuntimeError(f"Downloaded file failed SHA-256 validation: {target.name}")
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
    downloaded_at = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
    rows = []

    all_files = dict(FILES)
    all_files["PMC4796354_full_text.xml"] = PAPER_XML_URL
    for filename, url in all_files.items():
        target = SOURCE_DIR / filename
        disposition = fetch(url, target, args.refresh, EXPECTED_SHA256.get(filename))
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
                "license": "GPL (>= 2) (BatchCorrMetabolomics DESCRIPTION)",
                "used": "true",
                "source_commit": COMMIT if url.startswith(RAW_BASE) else "",
                "disposition": disposition,
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
