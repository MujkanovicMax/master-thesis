#!/bin/bash
INCLUDES="-I/home/m/Mujkanovic.Max/ma/rayli/include \
-isystem /usr/include/eigen3 \
-I/software/opt/bionic/x86_64/netcdf/4.7.0-gcc/include
"
LDFLAGS="\
	-L/software/opt/bionic/x86_64/netcdf/4.7.0-gcc/lib
	-L/home/m/Mujkanovic.Max/ma/rayli/build/lib \
	"
LIBS="-lpthread -lrayli_common -lnetcdf_c++4 -lnetcdf -lhdf5_hl -lhdf5 -lm -ldl -lz -lcurl"

CFLAGS="-std=c++17 -ggdb -O2" #-fsanitize=address -fsanitize=undefined"

g++ $CFLAGS -o panorama $INCLUDES alt_mc.cpp $LDFLAGS $LIBS
