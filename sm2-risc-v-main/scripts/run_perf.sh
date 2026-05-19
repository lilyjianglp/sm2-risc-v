#!/bin/bash
set -e

cd "$(dirname "$0")/../src"

echo "[INFO] Enable perf user access..."
sudo sysctl -w kernel.perf_user_access=2

echo "[INFO] Running perf benchmark..."
make clean
make all
make perf-all
