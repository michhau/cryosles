#!/usr/bin/env bash
# Run the three supported SGS cases one after another and archive each result.
# The launcher is invoked in foreground mode, so a subsequent case never starts
# until the previous mpirun exits successfully.
set -euo pipefail

model_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
archive_dir=${SGS_ARCHIVE_DIR:-"$model_dir/sgs_comparison/runs"}
nsteps=${SGS_NSTEPS:-800000}

if [[ "$archive_dir" != /* ]]; then
  archive_dir="$model_dir/$archive_dir"
fi

if ! [[ "$nsteps" =~ ^[1-9][0-9]*$ ]]; then
  echo "SGS_NSTEPS must be a positive integer; got: $nsteps" >&2
  exit 2
fi

if pgrep -u "$(id -u)" -f '(^|/)program\.exe([[:space:]]|$)' >/dev/null; then
  echo "A CryosLES program.exe process is already running; refusing to mix runs." >&2
  exit 1
fi

declare -a cases=(static dynamic lagrangian_scale_dependent)
declare -A labels=(
  [static]=static_Smagorinski
  [dynamic]=dynamic_Smagorinski
  [lagrangian_scale_dependent]=lagrangian_scale_dependent
)

declare -a run_artifacts=(
  output_0
  restart_files_0
  PARAMETERS.py_0
  patches.in
  hetero_srf_props.in
  ke.txt
  nohup.out
  run.log
  total_time.dat
  np_global.dat
  t_sgsFluxUpperBoundary_0.txt
  tstar_0.txt
  ustar_0.txt
  uvpheight.log
  wheight.log
)

mkdir -p -- "$archive_dir"

# Fail before launching anything if an archive would be overwritten.
for case_name in "${cases[@]}"; do
  destination="$archive_dir/${labels[$case_name]}"
  if [[ -e "$destination" ]]; then
    echo "Archive destination already exists: $destination" >&2
    echo "Choose another SGS_ARCHIVE_DIR or move the existing archive first." >&2
    exit 1
  fi
done

reset_run_artifacts() {
  local artifact
  for artifact in "${run_artifacts[@]}"; do
    rm -rf -- "$model_dir/$artifact"
  done
}

archive_run() {
  local destination=$1
  local artifact

  mkdir -- "$destination"
  for artifact in "${run_artifacts[@]}"; do
    if [[ -e "$model_dir/$artifact" ]]; then
      cp -a -- "$model_dir/$artifact" "$destination/"
    fi
  done
  printf 'case=%s\nnsteps=%s\n' "$case_name" "$nsteps" > "$destination/run_metadata.txt"
}

cd -- "$model_dir"
for case_name in "${cases[@]}"; do
  destination="$archive_dir/${labels[$case_name]}"
  echo "=== Preparing $case_name (${labels[$case_name]}) ==="

  # Remove only per-run products. Source, executable, and build artifacts stay intact.
  reset_run_artifacts
  SGS_NSTEPS="$nsteps" ./sgs_comparison/prepare_sgs_case.sh "$case_name"

  echo "=== Running $case_name ==="
  CRYOSLES_FOREGROUND=1 ./a_launch_pc.sh > nohup.out 2>&1

  echo "=== Archiving $case_name to $destination ==="
  archive_run "$destination"
done

echo "All SGS cases completed successfully. Archives: $archive_dir"
