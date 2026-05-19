#!/bin/bash
set -e

cd "$(dirname "$0")/../src"

echo "[INFO] Running SM2 correctness workflow..."
make clean
make correctness
