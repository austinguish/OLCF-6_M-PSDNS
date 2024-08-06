#!/bin/bash

source setUpModules_frontier.sh
export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}

module list

if [ ! -d ".srcmake" ]; then
   make -f makefile.frontier srcmake
fi

make -f makefile.frontier clean
make -f makefile.frontier

# check the executable in this environment
ldd DNS_PEN_GPU_p*.x
