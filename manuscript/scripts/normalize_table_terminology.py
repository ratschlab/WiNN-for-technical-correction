#!/usr/bin/env python3
"""Apply manuscript display terminology to generated TeX, leaving data keys intact.

Run after detailed table generation. The condensed-table assembler also calls this
step, so the published tables and supplementary source tables use the same terms.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def normalize(text):
    replacements = {
        "clean-reference": "ground-truth", "clean simulated matrix": "ground-truth simulated matrix",
        "corrected-to-clean": "corrected-to-ground-truth",
        "injected-effect ratio": "generating-effect ratio", "signed injected coefficient": "signed generating coefficient",
        "identity-collapsed": "sample-collapsed", "identity-cluster": "sample-cluster",
        "repeated-identity datasets": "datasets with repeated samples",
        "replicate/sample-identity agreement": "replicate/sample agreement",
        " & Residual GAM ": " & Residual drift ", " & GAM (": " & Drift (",
        "(GAM deviance explained) (proportion)": "(GAM deviance explained, proportion)",
        "imposed-effect ratio": "generating-effect ratio", "signed imposed coefficient": "signed generating coefficient",
        " & Imposed ratio & ": " & Generating ratio & ",
        "Clean-truth": "Ground-truth", "clean-truth": "ground-truth",
        "clean truth": "ground truth", "Clean simulation truth": "Simulation ground truth",
        "the clean matrix": "the ground-truth matrix",
        "Held-out reference CV": "Held-out QC CV", "held-out reference CV": "held-out QC CV",
        "Hidden-reference CV": "Held-out QC CV", "Hidden-ref.\\ CV": "Held-out QC CV",
        "GAM is residual run-order deviance explained": "Drift is residual run-order drift (GAM deviance explained)",
        "GAM is residual acquisition-order deviance explained": "Drift is residual run-order drift (GAM deviance explained)",
        "GAM is ": "Drift is ",
        "residual run-order deviance explained (proportion)": "residual run-order drift (GAM deviance explained, proportion)",
        "Residual GAM deviance": "Residual run-order drift", "residual GAM deviance": "residual run-order drift",
        " & GAM & ": " & Drift & ", "CV, GAM, batch": "CV, drift, batch",
        "categorical batch weighted-PC $R^2$": "variance explained by batch (weighted-PC $R^2$)",
        "batch $R^2$ is categorical weighted-PC $R^2$": "batch $R^2$ is variance explained by batch, calculated as variance-weighted PC $R^2$",
        "Batch $R^2$ is categorical weighted-PC $R^2$": "Batch $R^2$ is variance explained by batch, calculated as variance-weighted PC $R^2$",
        "within-accession": "within-plant-line", "accession-cluster": "plant-line-cluster",
        "repeated accession identity": "repeated measurements of plant lines",
        "neutral deposited group": "biological group", "deposited groups": "biological groups",
        "deposited biological endpoint": "biological endpoint",
        "deposited phenotype-group endpoint": "biological group endpoint",
        "Cross-batch identity discrimination": "Cross-batch sample matching",
        "Cross-batch identity retrieval": "Cross-batch sample matching",
        "biological-identity bootstrap": "biological-sample bootstrap",
        "within every resampled identity": "within every resampled biological sample",
        "within-identity feature-profile": "within-sample feature-profile",
        "346 identities": "346 plant lines", "267 identities": "267 biological samples",
        "repeated biological identities": "repeated biological samples",
        " & Clean Pearson & Clean ratio & Injected ratio & ": " & Truth $r$ & Truth ratio & Generating ratio & ",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text.replace("(GAM deviance explained) (proportion)", "(GAM deviance explained, proportion)")


def main():
    changed = 0
    for path in sorted((ROOT / "generated/supplementary_tables").glob("*.tex")):
        old = path.read_text()
        new = normalize(old)
        if new != old:
            path.write_text(new)
            changed += 1
    print(f"Updated display terminology in {changed} generated tables; numeric data unchanged.")


if __name__ == "__main__":
    main()
