#!/bin/bash
ulimit -c unlimited
ulimit -s unlimited

export OMP_NUM_THREADS=2

if [[ "${CRYOSLES_FOREGROUND:-0}" == "1" ]]; then
  exec mpirun -np 16 ./program.exe
fi

nohup mpirun -np 16 ./program.exe &
