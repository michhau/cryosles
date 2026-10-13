#!/usr/bin/env bash
# Stage one reproducible SGS-comparison case in ../PARAMETERS.py_0.
# Run from this directory, for example:
#   ./sgs_comparison/prepare_sgs_case.sh static
set -euo pipefail

nsteps=${SGS_NSTEPS:-800000}
if ! [[ "$nsteps" =~ ^[1-9][0-9]*$ ]]; then
  echo "SGS_NSTEPS must be a positive integer; got: $nsteps" >&2
  exit 2
fi

case_name=${1:-}
case "$case_name" in
  static)
    model=1
    description="static Smagorinsky"
    ;;
  dynamic)
    model=2
    description="plane-averaged dynamic Smagorinsky"
    ;;
  lagrangian_scale_dependent)
    model=5
    description="Lagrangian scale-dependent dynamic Smagorinsky"
    ;;
  *)
    echo "Usage: $0 {static|dynamic|lagrangian_scale_dependent}" >&2
    exit 2
    ;;
esac

sed \
  -e 's/^L_X=.*/L_X= 12.566370614359172 # 4*pi m/' \
  -e 's/^L_Y=.*/L_Y= 12.566370614359172 # 4*pi m/' \
  -e 's/^L_Z=.*/L_Z= 12.566370614359172 # 4*pi m/' \
  -e "s/^NSTEPS=.*/NSTEPS= ${nsteps} # total integration steps/" \
  -e "s/^MODEL=.*/MODEL= ${model}  # ${description}/" \
  -e 's/^LBC_DZ0=.*/LBC_DZ0= 0.01 # homogeneous momentum roughness length (m)/' \
  -e 's/^C_CNT2=.*/C_CNT2= 10 # 40 Hz centreline output with DT=0.0025 s/' \
  20_PARAMETERS.py_0 > PARAMETERS.py_0

cp sgs_comparison/patches_homogeneous_32.in patches.in
cp sgs_comparison/hetero_srf_props_homogeneous_0.01m.in hetero_srf_props.in

echo "Prepared ${case_name}: ${description}"
echo "Grid/domain are compiled as 32^3 on 16 MPI ranks; recompile before running."
