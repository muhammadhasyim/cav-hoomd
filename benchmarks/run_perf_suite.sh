#!/usr/bin/env bash
# Short perf benchmarks for aging_weak_lambda production stack (~30 ps each).
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPT="$REPO/examples/05_advanced_run.py"
GSD="$REPO/square_lambda0.025_diffeq/prod-0.gsd"
OUT="$REPO/benchmarks/results"
export PYTHONPATH=/opt/hoomd-md:${PYTHONPATH:-}

mkdir -p "$OUT"

COMMON=(
  --molecular-bath bussi --cavity-bath langevin
  --lambda-coupling 0.03
  --temperature 100.0 --frequency 1560.0
  --runtime 30.0 --switch-time 200.0
  --input-gsd "$GSD" --frame -1
  --fixed-timestep --timestep 1.0
  --disable-gsd --console-output-period-ps 10.0
  --replicas 0 --seed 1
)

run_case() {
  local name="$1"
  shift
  local dir="$OUT/$name"
  mkdir -p "$dir"
  echo "=== Running $name ==="
  (
    cd "$dir"
    CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}" \
    python3 "$SCRIPT" "${COMMON[@]}" "$@" > run.log 2>&1
  ) || true
  grep -E '^\s+[0-9]+\s+[0-9]+\.[0-9]+\s+[0-9]+\.[0-9]+\s+[0-9]+\.[0-9]+' "$dir/run.log" | tail -1 \
    | awk -v n="$name" '{printf "%s tps=%s ns_per_day=%s sim_ps=%s\n", n, $2, $4, $3}' \
    >> "$OUT/summary.txt"
}

: > "$OUT/summary.txt"

# A: full production stack
run_case A_baseline \
  --coupling-type step --device GPU --gpu-id 0 \
  --enable-energy-tracker --energy-output-period-ps 1.0 \
  --enable-fkt --fkt-kmag 6.02 --fkt-wavevectors 50 \
  --fkt-ref-interval 200.0 --fkt-max-refs 13 --fkt-output-period-ps 1.0

# B: physics only (no trackers / HDF5)
run_case B_minimal \
  --coupling-type step --device GPU --gpu-id 0 \
  --disable-hdf5-output --disable-temp-tracker

# C: no energy tracker only
run_case C_no_energy \
  --coupling-type step --device GPU --gpu-id 0 \
  --enable-fkt --fkt-kmag 6.02 --fkt-wavevectors 50 \
  --fkt-ref-interval 200.0 --fkt-max-refs 13 --fkt-output-period-ps 1.0

# D: constant coupling (no StepVariant Python callback)
run_case D_constant_coupling \
  --coupling-type constant --device GPU --gpu-id 0 \
  --enable-energy-tracker --energy-output-period-ps 1.0 \
  --enable-fkt --fkt-kmag 6.02 --fkt-wavevectors 50 \
  --fkt-ref-interval 200.0 --fkt-max-refs 13 --fkt-output-period-ps 1.0

# E: CPU reference
run_case E_cpu \
  --coupling-type step --device CPU \
  --enable-energy-tracker --energy-output-period-ps 1.0 \
  --enable-fkt --fkt-kmag 6.02 --fkt-wavevectors 50 \
  --fkt-ref-interval 200.0 --fkt-max-refs 13 --fkt-output-period-ps 1.0

echo "=== Summary ==="
cat "$OUT/summary.txt"
