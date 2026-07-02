#!/usr/bin/env bash
# fix_renv_rpaths.sh
#
# Workaround for: Positron's R kernel process does not inherit LD_LIBRARY_PATH
# from the r-conda conda environment, so compiled R packages that depend on
# conda-provided shared libraries (arrow -> libarrow_substrait.so,
# sf -> libproj.so) fail to load with "cannot open shared object file" even
# though the .so files are present on disk.
#
# This script finds the current renv cache locations of arrow.so and sf.so
# (their paths are keyed by a content hash that changes on reinstall) and
# uses patchelf to bake the conda environment's lib directory into each
# file's RPATH, so the dynamic linker finds the dependencies without relying
# on LD_LIBRARY_PATH at all.
#
# Run this after any renv::restore() / renv::hydrate() that reinstalls
# arrow or sf. It is idempotent -- safe to run repeatedly.

set -euo pipefail

CONDA_LIB="/home/mtalbot/software/miniconda3/envs/r-conda/lib"
RENV_CACHE="$HOME/.cache/R/renv/cache"

if ! command -v patchelf >/dev/null 2>&1; then
  echo "patchelf not found on PATH. Install it with:"
  echo "  conda install -c conda-forge patchelf -n r-conda"
  exit 1
fi

if [ ! -d "$CONDA_LIB" ]; then
  echo "Conda lib directory not found: $CONDA_LIB"
  exit 1
fi

patched_any=0

for target in arrow.so sf.so; do
  # Find all matches -- there may be more than one if multiple versions
  # are cached, so patch every match found.
  while IFS= read -r -d '' so_path; do
    echo "Patching RPATH: $so_path"
    patchelf --set-rpath "$CONDA_LIB" "$so_path"
    patched_any=1
  done < <(find "$RENV_CACHE" -type f -name "$target" -print0 2>/dev/null)
done

if [ "$patched_any" -eq 0 ]; then
  echo "No arrow.so or sf.so found under $RENV_CACHE."
  echo "Has renv::restore() / renv::hydrate() been run yet?"
  exit 1
fi

echo "Done. Verify in R with: library(arrow); library(sf)"