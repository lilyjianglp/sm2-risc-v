#!/bin/bash

cd "$(dirname "$0")/.."

echo "[INFO] Counting source lines..."

find src scripts \
  -type f \( -name "*.c" -o -name "*.h" -o -name "*.S" -o -name "*.s" -o -name "*.sh" -o -name "Makefile" \) \
  -print0 | xargs -0 wc -l
