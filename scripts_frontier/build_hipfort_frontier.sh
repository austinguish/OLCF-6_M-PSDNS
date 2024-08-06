#!/bin/bash 
#set -eux -o pipefail

# For ease all of this should be run in its own directory
export BASE_DIR=${PWD}

echo ${BASE_DIR}

#exit

### cce/15 is the current default on Frontier
export MYCOMPILERS="CCE-15.0.0"
export VERSION="5.4.0"

## load rest of modules
module load git cmake rocm/${VERSION} craype-accel-amd-gfx90a

export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}

module list

export CC=`which cc`
export FC=`which ftn`
export CXX=`which CC`

export AR=`which ar`
export RANLIB=`which ranlib`


echo $CC
echo $CXX
echo $FC
echo $AR
echo $RANLIB

#exit


#Build hipfort
echo ""
echo "BUILDING HIPFORT ... "
echo ""

## local install ... CAN BE CUSTOMIZED
export MY_HIPFORT_DIR=${BASE_DIR}/hipfort_v${VERSION}_${MYCOMPILERS}

## local install ... CAN BE CUSTOMIZED ... MUST BE USED WHEN BUILDING HIPFFT DNS CODE
export HIPFORT_INSTALL="${BASE_DIR}/INSTALLATION/hipfort_v${VERSION}_${MYCOMPILERS}"

# setting stuff to build hipfort
export ROCM="${OLCF_ROCM_ROOT}"
export HIPFORT_HOME="${HIPFORT_INSTALL}"
export HIPFORT_ARCH="amdgcn"
export HIPFORT_ARCHGPU="amdgcn-gfx90a"
export COMPILER_FLAGS="-ffree -homp -N1023 -eT" # -eT is needed to preprocess, can also use -eZ (gives *i files)


if [ ! -d "${MY_HIPFORT_DIR}" ]; then
   git clone https://github.com/ROCmSoftwarePlatform/hipfort.git
   mv hipfort ${MY_HIPFORT_DIR}
   cd ${MY_HIPFORT_DIR}

   # get a specific version
   git checkout rocm-${VERSION}

   cd ${BASE_DIR}
else
   echo " Directory hipfort exists"
fi

#exit

cd ${MY_HIPFORT_DIR}

if [ ! -d "build" ]; then
   mkdir build
fi
cd build

# use cce/15 compilers
# NOTE: cmake adds the flags in CMAKE_Fortran_FLAGS_<RELEASE,TESTING,DEBUG> to the flags
#       in COMPILER_FLAGS to form the complete set of compiling flags (due to the 
#       instructions in CMakeLists.txt and all of the *.cmake files)
cmake ../ -DHIPFORT_COMPILER=ftn \
          -DHIPFORT_COMPILER_FLAGS="${COMPILER_FLAGS}" \
          -DHIPFORT_INSTALL_DIR="${HIPFORT_INSTALL}" \
          -DHIPFORT_AR="${AR}" \
          -DHIPFORT_RANLIB="${RANLIB}" \
          -DHIPFORT_BUILD_TYPE="RELEASE" \
          -DCMAKE_Fortran_FLAGS_RELEASE="-O3" \
          -DCMAKE_Fortran_FLAGS_TESTING="-O2" \
          -DCMAKE_Fortran_FLAGS_DEBUG="-g -O0"\
          -DCMAKE_INSTALL_PREFIX="${HIPFORT_INSTALL}"

#exit

make -j8 |& tee my_hipfort_build
make install |& tee my_hipfort_install


echo ""
echo "FINISHED BUILDING HIPFORT"
echo ""
