#!/bin/bash -l

#SBATCH --job-name="snowshine_with_omp_2_all"
#SBATCH --output=snowshine.%j.o
#SBATCH --error=snowshine.%j.e
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=2
#SBATCH --constraint=gpu
#SBATCH --time=24:00:00
##SBATCH --partition=debug
#SBATCH --account=s569
#SBATCH --mail-type=ALL
#SBATCH --mail-user=varun.sharma@epfl.ch

module load slurm
module load daint-gpu

module switch PrgEnv-cray PrgEnv-intel


#mkdir turbine_data
#mkdir wind_farm_forces
#mkdir wind_farm_setup

#cp ../default_turbine.txt .

 #module load Score-P/3.0-CrayIntel-2016.11 
 #module load Scalasca/2.3.1-CrayIntel-2016.11 
 #module load cray-libsci

 #export SCOREP_EXPERIMENT_DIRECTORY=scorep_sum_trace
 #export SCOREP_ENABLE_TRACING=true
 #export SCOREP_ENABLE_PROFILING=true
 #export SCOREP_TOTAL_MEMORY=1024M

ulimit -s unlimited
ulimit -c unlimited
export OMP_NUM_THREADS=2

#export KMP_AFFINITY=scatter

#numactl --interleave=all
#/usr/bin/time -p srun --distribution=block:block --exclusive --hint=nomultithread numactl --interleave=all ./program.exe
/usr/bin/time -p srun --hint=nomultithread ./program.exe
#srun ./program.exe

echo " ****** SIMULATION FINISHED ! ****** "
exit

