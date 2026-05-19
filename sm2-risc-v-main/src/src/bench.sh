#!/usr/bin/env bash
set -euo pipefail

# ====== config: adapt to your environment ======
ITER="${ITER:-20}"

SYSROOT="${SYSROOT:-/home/wen/riscv-sysroot}"
QEMU="${QEMU:-qemu-riscv64}"

BIN_C="${BIN_C:-./test_sm2_kex_c}"
BIN_RV="${BIN_RV:-./test_sm2_kex_rv}"
BIN_MONT="${BIN_MONT:-./test_sm2_kex_mont_rv}"

# If your outputs are in different folder, cd there or adjust BIN_*
# ==============================================

need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Missing command: $1"; exit 1; }; }
need_cmd awk
need_cmd grep
need_cmd bc

run_bin() {
  local bin="$1"
  # always run via qemu with sysroot
  "$QEMU" -L "$SYSROOT" "$bin"
}

# Extract the "mean=" value from a given KEY line.
# Works for:
#   "fp_mul (Precise) : mean=221, min=..."
#   "fp_mul (avg cycles/op) : mean=452, ..."
extract_mean_once() {
  local bin="$1"
  local key="$2"

  run_bin "$bin" \
    | awk -v k="$key" '
        index($0, k) {
          # find "mean=" then capture number until comma/space
          m = match($0, /mean=[0-9]+(\.[0-9]+)?/)
          if (m) {
            s = substr($0, RSTART+5, RLENGTH-5)
            print s
            exit
          }
        }
      '
}

# Average mean value across ITER runs; returns a float with 2 decimals.
avg_mean() {
  local bin="$1"
  local key="$2"
  local sum="0"

  for _ in $(seq 1 "$ITER"); do
    local v
    v="$(extract_mean_once "$bin" "$key" || true)"
    if [[ -z "${v:-}" ]]; then
      echo "ERROR: cannot find mean for key [$key] in $bin" >&2
      echo "Tip: run once to inspect: $QEMU -L $SYSROOT $bin" >&2
      exit 1
    fi
    sum="$(echo "$sum + $v" | bc -l)"
  done

  echo "$(echo "scale=2; $sum / $ITER" | bc -l)"
}

# Print integer if it is x.00 else print float
pretty_num() {
  local x="$1"
  if echo "$x" | grep -qE '^[0-9]+(\.0+)?$'; then
    # integer-ish
    printf "%.0f" "$x"
  else
    printf "%s" "$x"
  fi
}

# Speedup C / X (float)
speedup() {
  local c="$1"
  local x="$2"
  echo "$(echo "scale=2; $c / $x" | bc -l)"
}

echo "Running each version $ITER times..."
echo "QEMU=$QEMU  SYSROOT=$SYSROOT"
echo ""

echo "===== Basic Arithmetic ====="
echo ""
printf "%-12s %-10s %-10s %-10s\n" "Operation" "C" "RV ASM" "Mont"

C_reduce="$(avg_mean "$BIN_C"   "fp_reduce")"
RV_reduce="$(avg_mean "$BIN_RV" "fp_reduce")"
M_reduce="$(avg_mean "$BIN_MONT" "fp_reduce")"
printf "%-12s %-10s %-10s %-10s\n" "fp_reduce" \
  "$(pretty_num "$C_reduce")" "$(pretty_num "$RV_reduce")" "$(pretty_num "$M_reduce")"

# 兼容：你某些版本可能输出 "fp_mul (Precise)"，纯C可能输出 "fp_mul (avg cycles/op)"
# 所以用一个 key 列表做 fallback（不改你的代码也能适配）
avg_mean_fallback() {
  local bin="$1"; shift
  for key in "$@"; do
    if run_bin "$bin" | grep -q "$key"; then
      avg_mean "$bin" "$key"
      return 0
    fi
  done
  echo "ERROR: none of keys matched in $bin: $*" >&2
  exit 1
}

C_mul="$(avg_mean_fallback "$BIN_C" "fp_mul (Precise)" "fp_mul (avg cycles/op)" "fp_mul")"
RV_mul="$(avg_mean_fallback "$BIN_RV" "fp_mul (Precise)" "fp_mul (avg cycles/op)" "fp_mul")"
M_mul="$(avg_mean_fallback "$BIN_MONT" "fp_mul (Precise)" "fp_mul (avg cycles/op)" "fp_mul")"
printf "%-12s %-10s %-10s %-10s\n" "fp_mul" \
  "$(pretty_num "$C_mul")" "$(pretty_num "$RV_mul")" "$(pretty_num "$M_mul")"

C_sqr="$(avg_mean_fallback "$BIN_C" "fp_sqr (Precise)" "fp_sqr (avg cycles/op)" "fp_sqr")"
RV_sqr="$(avg_mean_fallback "$BIN_RV" "fp_sqr (Precise)" "fp_sqr (avg cycles/op)" "fp_sqr")"
M_sqr="$(avg_mean_fallback "$BIN_MONT" "fp_sqr (Precise)" "fp_sqr (avg cycles/op)" "fp_sqr")"
printf "%-12s %-10s %-10s %-10s\n" "fp_sqr" \
  "$(pretty_num "$C_sqr")" "$(pretty_num "$RV_sqr")" "$(pretty_num "$M_sqr")"

echo ""
echo "===== High-Level Ops ====="
echo ""
printf "%-12s %-10s %-10s %-10s\n" "Operation" "C" "RV ASM" "Mont"

C_scalar="$(avg_mean "$BIN_C"   "ScalarMul")"
RV_scalar="$(avg_mean "$BIN_RV" "ScalarMul")"
M_scalar="$(avg_mean "$BIN_MONT" "ScalarMul")"
printf "%-12s %-10s %-10s %-10s\n" "ScalarMul" \
  "$(pretty_num "$C_scalar")" "$(pretty_num "$RV_scalar")" "$(pretty_num "$M_scalar")"

C_kex="$(avg_mean_fallback "$BIN_C" "KEX Initiator (full verify)" "KEX Initiator" "KEX Verify")"
RV_kex="$(avg_mean_fallback "$BIN_RV" "KEX Initiator (full verify)" "KEX Initiator" "KEX Verify")"
M_kex="$(avg_mean_fallback "$BIN_MONT" "KEX Initiator (full verify)" "KEX Initiator" "KEX Verify")"
printf "%-12s %-10s %-10s %-10s\n" "KEX Verify" \
  "$(pretty_num "$C_kex")" "$(pretty_num "$RV_kex")" "$(pretty_num "$M_kex")"

echo ""
echo "===== Speedup (vs C) ====="
echo ""

printf "fp_mul     RV ASM x%s   Mont x%s\n" \
  "$(speedup "$C_mul" "$RV_mul")" \
  "$(speedup "$C_mul" "$M_mul")"

printf "ScalarMul  RV ASM x%s   Mont x%s\n" \
  "$(speedup "$C_scalar" "$RV_scalar")" \
  "$(speedup "$C_scalar" "$M_scalar")"

printf "KEX Verify RV ASM x%s   Mont x%s\n" \
  "$(speedup "$C_kex" "$RV_kex")" \
  "$(speedup "$C_kex" "$M_kex")"