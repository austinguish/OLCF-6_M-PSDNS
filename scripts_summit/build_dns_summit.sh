#!/bin/bash

source setUpModules_summit.sh

module list

if [ ! -d ".srcmake" ]; then
   make -f makefile.summit srcmake
fi

make -f makefile.summit clean
make -f makefile.summit

# check the executable in this environment
ldd DNS_PEN_GPU_p4.x
