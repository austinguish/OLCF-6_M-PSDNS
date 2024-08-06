#!/bin/bash 

# For ease all of this should be run in its own directory
export BASE_DIR=${PWD}

echo ${BASE_DIR}

module load cuda/11.4.2
module load hip-cuda/5.1.0
module load git cmake

module list

#exit

export CC=`which xlc`
export FC=`which xlcuf`
export CXX=`which xlc`

export AR=`which ar`
export RANLIB=`which ranlib`


echo $CC
echo $CXX
echo $FC
echo $AR
echo $RANLIB

export MYCOMPILERS="IBM-XL"
export VERSION="5.1.0"

#exit


#Build hipfort
echo ""
echo "BUILDING HIPFORT ... "
echo ""

## local build ... CAN BE CUSTOMIZED 
export MY_HIPFORT_DIR=${BASE_DIR}/hipfort_v${VERSION}_${MYCOMPILERS}

## local install ... CAN BE CUSTOMIZED ... MUST BE USED WHEN BUILDING HIPFFT DNS CODE 
export HIPFORT_INSTALL="${BASE_DIR}/INSTALLATION/hipfort_v${VERSION}/${MYCOMPILERS}"

# setting stuff to build hipfort
export HIPFORT_HOME="${HIPFORT_INSTALL}"
export HIPFORT_ARCH="nvptx"
export HIPFORT_ARCHGPU="sm70"
export COMPILER_FLAGS="-qfree -qpreprocess -qsmp=omp -qoffload -qcuda"

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

# use ibm xlcuf compilers
# NOTE: cmake adds the flags in CMAKE_Fortran_FLAGS_<RELEASE,TESTING,DEBUG> to the flags
#       in COMPILER_FLAGS to form the complete set of compiling flags (due to the 
#       instructions in CMakeLists.txt and all of the *.cmake files)
cmake ../ -DHIPFORT_COMPILER=xlcuf \
          -DHIPFORT_INSTALL_DIR="${HIPFORT_INSTALL}" \
          -DHIPFORT_AR="${AR}" \
          -DHIPFORT_RANLIB="${RANLIB}" \
          -DHIPFORT_BUILD_TYPE="RELEASE" \
          -DHIPFORT_COMPILER_FLAGS="-qsmp=omp -qoffload -qcuda" \
          -DCMAKE_Fortran_FLAGS="-qfree -qpreprocess" \
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
