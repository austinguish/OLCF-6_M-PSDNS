#!/bin/bash

BASE_DIR=${PWD}
cd ${BASE_DIR}

## set CODE_PATH to point to your build directory
export CODE_PATH=/ccs/home/nichols/PROJECTS/OLCF6_BENCHMARK/hipfft-dns-benchmark

## set your environment
source ${CODE_PATH}/setUpModules_frontier.sh
module list

export HIP_LAUNCH_BLOCKING=1

## this is just the name of the executable ... it is copied to NVMe in the next section
## NOTE:  do NOT include the path here
exec=DNS_PEN_GPU_p4.x

#ldd ${CODE_PATH}/$exec
ls -ltr ${CODE_PATH}/$exec
