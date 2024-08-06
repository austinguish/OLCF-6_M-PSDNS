#!/bin/bash

#SBATCH -A STF006
#SBATCH -J benchmark
#SBATCH -o stdout
#SBATCH -e stderr
#SBATCH -t 00:20:00
#SBATCH -p batch
#SBATCH -N 4096
##SBATCH --threads-per-core=2
#SBATCH -C nvme

BASE_DIR=${PWD}
cd ${BASE_DIR}

if [ ! -d "MPI_timings" ]; then
   mkdir MPI_timings
fi

## set CODE_PATH to point to your build directory
export CODE_PATH=/ccs/home/nichols/PROJECTS/OLCF6_BENCHMARK/hipfft-dns-benchmark

## set your environment
source ${CODE_PATH}/setUpModules_frontier.sh
module list

export HIP_LAUNCH_BLOCKING=1

## this is just the name of the executable ... it is copied to NVMe in the next section
## NOTE:  do NOT include the path here
exec=DNS_PEN_GPU_p4.x

# SBCAST executable from filesystem to NVMe 
# NOTE: ``-C nvme`` is needed in SBATCH headers to use the NVMe drive
# NOTE: dlopen'd files will NOT be picked up by sbcast
sbcast --send-libs --exclude=NONE -pf ${CODE_PATH}/${exec} /mnt/bb/$USER/${exec}
if [ ! "$?" == "0" ]; then
    # CHECK EXIT CODE. When SBCAST fails, it may leave partial files on the compute nodes, and if you continue to launch srun,
    # your application may pick up partially complete shared library files, which would give you confusing errors.
    echo "SBCAST failed!"
    exit 1
fi

# Some applications may expect a generic name of certain libraries (ie, `libamdhip64.so` instead of `libamdhip64.so.5`)
# You may need to use an `srun` line like this to add sym-links on each node
# The 2 most common cases are `libhsa-runtime64.so` and `libamdhip64.so`
 srun -N ${SLURM_NNODES} -n ${SLURM_NNODES} --ntasks-per-node=1 --label -D /mnt/bb/$USER/${exec}_libs \
   bash -c "if [ -f libhsa-runtime64.so.1 ]; then ln -s libhsa-runtime64.so.1 libhsa-runtime64.so; fi;
            if [ -f libamdhip64.so.5 ]; then ln -s libamdhip64.so.5 libamdhip64.so; fi;
            if [ -f librocfft.so.0 ]; then ln -s librocfft.so.0 librocfft.so; fi;
            if [ -f libroctracer64.so.4 ]; then ln -s libroctracer64.so.4 libroctracer64.so; fi"

# Check to see if file exists
echo "*****ls -lh /mnt/bb/$USER*****"
ls -lh /mnt/bb/$USER/
echo "*****ls -lh /mnt/bb/$USER/${exec}_libs*****"
ls -lh /mnt/bb/$USER/${exec}_libs

# All required libraries are in 2 total directories - this can have significant benefit at scale
export LD_LIBRARY_PATH="/mnt/bb/$USER/${exec}_libs:$(pkg-config --variable=libdir libfabric)"

## check which rocfft libs are being used
echo "*****ldd /mnt/bb/$USER/${exec}*****"
ldd /mnt/bb/$USER/${exec}
echo "*************************************"

## WITHOUT -S0 (ie. with default of -S8)
export OMP_NUM_THREADS=7
srun -u -N4096 -n32768 -c7 --ntasks-per-gpu=1 --gpu-bind=closest /mnt/bb/$USER/$exec

#export OMP_NUM_THREADS=14
#srun -u -N4096 -n32768 -c14 --ntasks-per-gpu=1 --gpu-bind=closest /mnt/bb/$USER/$exec
