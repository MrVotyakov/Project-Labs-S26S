#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CASSANDRA_EXE="$PROJECT_ROOT/Cassandra/bin/cassandra.exe"
LIB_DIR="$PROJECT_ROOT/Cassandra/Libraries/locals/gfortran/local/lib"
MONITOR="$PROJECT_ROOT/MC/data/monitor_gemc_progress.py"

cd "$SCRIPT_DIR"

export DYLD_LIBRARY_PATH="$LIB_DIR${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"

if [[ ! -x "$CASSANDRA_EXE" ]]; then
  echo "Cassandra executable not found or not executable: $CASSANDRA_EXE" >&2
  exit 1
fi

if [[ ! -f "$MONITOR" ]]; then
  echo "Progress monitor not found: $MONITOR" >&2
  exit 1
fi

if [[ "${1:-}" != "" && "${1:-}" != "--vle-only" ]]; then
  echo "Usage: $0 [--vle-only]" >&2
  exit 2
fi

run_stage() {
  local input_file="$1"
  local run_prefix="$2"

  if [[ ! -f "$input_file" ]]; then
    echo "Input file not found: $SCRIPT_DIR/$input_file" >&2
    exit 1
  fi

  echo
  echo "Starting Cassandra stage: $input_file"
  echo "Progress prefix: $run_prefix"

  "$CASSANDRA_EXE" "$input_file" &
  local cassandra_pid=$!
  echo "Cassandra PID: $cassandra_pid"

  python3 "$MONITOR" --prefix "$run_prefix" --pid "$cassandra_pid" --interval 2 &
  local monitor_pid=$!

  local status=0
  wait "$cassandra_pid" || status=$?
  wait "$monitor_pid" 2>/dev/null || true

  if [[ "$status" -ne 0 ]]; then
    echo "Cassandra stage failed with status $status: $input_file" >&2
    exit "$status"
  fi

  if [[ -f "${run_prefix}.log" ]] && grep -Eq "Fatal Error|\\*+ERROR\\*+" "${run_prefix}.log"; then
    echo "Cassandra reported a fatal error: ${run_prefix}.log" >&2
    exit 1
  fi

  echo "Cassandra stage finished: $input_file"
}

require_liquid_like_box1() {
  local property_file="$1"
  local min_density="450.0"

  if [[ ! -f "$property_file" ]]; then
    echo "Liquid-phase property file was not created: $property_file" >&2
    exit 1
  fi

  awk -v min_density="$min_density" '
    $1 !~ /^#/ && NF >= 5 {
      step = $1
      pressure = $3
      density = $5
      found = 1
    }
    END {
      if (!found) {
        print "No readable property rows in " FILENAME > "/dev/stderr"
        exit 1
      }

      printf "Pre-equilibration box1 final record: MC_STEP %s P=%g bar rho=%g kg/m3\n", step, pressure, density

      if (density < min_density) {
        printf "Stopping before VLE stage: box1 density %g kg/m3 is below liquid-like guard %g kg/m3.\n", density, min_density > "/dev/stderr"
        exit 1
      }
    }
  ' "$property_file"
}

if [[ "${1:-}" == "--vle-only" ]]; then
  if [[ ! -f "azeotrope_gemc_npt_3000_preswap_eq.chk" ]]; then
    echo "Pre-equilibration checkpoint is missing." >&2
    exit 1
  fi

  echo "Using existing pre-equilibration checkpoint for VLE stage."
  run_stage "azeotrope_gemc_npt_600_50k.inp" "azeotrope_gemc_npt_3000_50k"
  exit 0
fi

run_stage "azeotrope_gemc_npt_3000_preswap_eq.inp" "azeotrope_gemc_npt_3000_preswap_eq"

if [[ ! -f "azeotrope_gemc_npt_3000_preswap_eq.chk" ]]; then
  echo "Pre-equilibration checkpoint was not created." >&2
  exit 1
fi

require_liquid_like_box1 "azeotrope_gemc_npt_3000_preswap_eq.box1.prp"

run_stage "azeotrope_gemc_npt_600_50k.inp" "azeotrope_gemc_npt_3000_50k"

echo
echo "Two-stage azeotrope run finished."
