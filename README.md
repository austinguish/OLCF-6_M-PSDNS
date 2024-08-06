# OLCF-6 Benchmark M-PSDNS Code

## History of the Algorithm
This benchmark code is a streamlined, minimalist version of a code developed by Professor P.K Yeung of the Georgia Institute of Technology. 

## Direct Numerical Simulation of Turbulence
Fluid turbulence is generally characterized as nonlinear, unsteady, disorderly fluctuations occurring across a vast range of scales in three-dimensions.  The Minimalist Pseudo-Spectral DNS (M-PSDNS) code is designed on top of a purpose-built 3D FFT algorithm to investigate the fundamental behavior of turbulence at very high resolution.  As is common with 3D FFT computations, MPI communications dominate the runtime (~76%) followed by the actual FFT computations (~15%) and by the packing/unpacking operations on the MPI send/receive buffers (~5%).  It is not uncommon for these three aspects of a 3D FFT algorithm to consume in excess of 96% of the overall runtime of the simulations.  Consequently, the M-PSDNS algorithm is an excellent bellwether of system performance for all scientific fields that must address large 3D FFTs or intense, large-scale MPI communications.  

## Requirements
Fortran compiler (2003 compliant) with OpenMP Offloading, MPI (preferably GPU-Aware), hipfort, and either CUDA or ROCM

## Documentation
Read the "M-PSDNS_Documentation.pdf" file in the "documentation" directory for a full description of the code.

## Run Rules
Read the "M-PSDNS_RunRules.pdf" file in the "documentation" directory for list of rules for evaluating the M-PSDNS code.




