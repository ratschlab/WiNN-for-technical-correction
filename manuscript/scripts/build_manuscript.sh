#!/usr/bin/env bash
# Build both documents with references resolved. Requires pdfLaTeX and BibTeX.
set -euo pipefail
cd "$(dirname "$0")/.."
build_dir=build/publication
mkdir -p "$build_dir"
for document in main supplementary; do
  pdflatex -interaction=nonstopmode -halt-on-error -output-directory="$build_dir" "$document.tex" > "$build_dir/${document}_build_stdout.txt" 2>&1
  bibtex "$build_dir/$document" > "$build_dir/${document}_bibtex_stdout.txt" 2>&1
  for pass in 1 2; do
    pdflatex -interaction=nonstopmode -halt-on-error -output-directory="$build_dir" "$document.tex" > "$build_dir/${document}_build_stdout.txt" 2>&1
  done
  cp "$build_dir/$document.pdf" "$document.pdf"
  echo "Built $document.pdf"
done
