#!/bin/bash

## unload darshan to use cray-mpich/8.1.25
module unload darshan-runtime/3.4.0

## cce/15.0.0 is the current default on Frontier
#module load cce/15.0.0
module load cray-mpich/8.1.25
module load rocm/5.4.0
module load craype-accel-amd-gfx90a

## These must be set before running with GPU-Aware Cray-MPICH
export MPICH_GPU_SUPPORT_ENABLED=1

## env var to avoid issues with big jobs 
export ROCFFT_RTC_CACHE_PATH=/dev/null

