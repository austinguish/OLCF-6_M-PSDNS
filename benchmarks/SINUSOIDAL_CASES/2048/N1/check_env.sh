#!/bin/bash

export CODE_PATH=/ccs/home/nichols/PROJECTS/OLCF6_BENCHMARK/hipfft-dns-benchmark
source $CODE_PATH/setUpModules_frontier.sh

module list

export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}

export HIP_LAUNCH_BLOCKING=1

exec=$CODE_PATH/DNS_PEN_GPU_p4.x

ls -ltr $exec
