#!/bin/bash

if [ ! -d ".srcmake" ]; then
   make -f makefile.lonestar srcmake
fi

make -f makefile.lonestar clean
make -f makefile.lonestar

# check the executable in this environment
ldd DNS_PEN_GPU_p*.x
