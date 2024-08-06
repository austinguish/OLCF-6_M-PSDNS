#!/bin/bash

#SBATCH -A STF006_FRONTIER
#SBATCH -J dns-pencil
#SBATCH -o stdout
#SBATCH -e stderr
#SBATCH -t 00:10:00
#SBATCH -p batch
#SBATCH -N 128
#SBATCH --threads-per-core=1
##SBATCH --threads-per-core=2

export CODE_PATH=/ccs/home/nichols/PROJECTS/OLCF6_BENCHMARK/hipfft-dns-benchmark
source $CODE_PATH/setUpModules_frontier.sh

module list

export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}

export HIP_LAUNCH_BLOCKING=1

exec=$CODE_PATH/DNS_PEN_GPU_p8.x

if [ ! -d "MPI_timings" ]; then
   mkdir MPI_timings
fi

## check which rocfft libs are being used
ldd $exec | grep rocfft

## add --unbuffered for I/O when debugging

## for --threads-per-core=1
export OMP_NUM_THREADS=7
srun -N128 -n1024 -c7 --ntasks-per-gpu=1 --gpu-bind=closest $exec


## for --threads-per-core=2
#export OMP_NUM_THREADS=14
#srun -N128 -n1024 -c14 --ntasks-per-gpu=1 --gpu-bind=closest $exec

