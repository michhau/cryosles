#!/bin/bash -l

#SBATCH --job-name="cryosles"
#SBATCH --account=s1329
#SBATCH --output=jobout_err/jobout_%j.out
#SBATCH --error=jobout_err/joberr_%j.err
#SBATCH --partition=normal
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=4
#SBATCH --constraint=mc
#SBATCH --time=24:00:00
##SBATCH --mail-type=ALL
##SBATCH --mail-user=michael.haugeneder@slf.ch

module load cray/23.12
module load PrgEnv-intel/8.5.0
#module load slurm
#module load daint-gpu

#module switch PrgEnv-cray PrgEnv-intel

 #module load Score-P/3.0-CrayIntel-2016.11 
 #module load Scalasca/2.3.1-CrayIntel-2016.11 
 #module load cray-libsci

ulimit -s unlimited
ulimit -c unlimited
export OMP_NUM_THREADS=4

#/usr/bin/time -p srun --distribution=block:block --exclusive --hint=nomultithread numactl --interleave=all ./program.exe
#/usr/bin/time -p srun --hint=nomultithread ./program.exe
/usr/bin/time -p srun --hint=nomultithread --cpu-bind=sockets ./program.exe

echo " ****** SIMULATION FINISHED ! ****** "
exit

