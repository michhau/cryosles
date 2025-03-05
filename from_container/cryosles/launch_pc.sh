#!/bin/bash
ulimit -c unlimited
ulimit -s unlimited

export OMP_NUM_THREADS=2

nohup mpirun -np 16 ./program.exe &
