#!/bin/bash


if [ ! -d ".srcmake" ]; then
   make -f makefile.nvhpc srcmake
fi

make -f makefile.nvhpc clean
make -f makefile.nvhpc

# check the executable in this environment
ldd DNS_PEN_GPU_p4.x
