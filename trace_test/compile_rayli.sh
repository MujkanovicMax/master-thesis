#!/bin/bash
INCLUDES="-I../rayli/include \
-isystem /usr/include/eigen3 \
-I../rayli/build/_deps/fmt-src/include/ \
-I../rayli/build/_deps/spdlog-src/include/ \
"

INCLUDES="$INCLUDES $(pkg-config netcdf-cxx4 --cflags)"
LDFLAGS="\
	-L../rayli/build/lib \
	"

LDFLAGS="$LDFLAGS $(pkg-config hdf5 --libs)"
LDFLAGS="$LDFLAGS $(pkg-config netcdf-cxx4 --libs)"
LIBS="-lrayli_common -lnetcdf_c++4 -lnetcdf -lhdf5_hl -lhdf5 -lm -ldl -lz -lcurl"

CFLAGS="-std=c++17 -ggdb -O3" #-fsanitize=address -fsanitize=undefined"

g++ $CFLAGS -o trace_optical_thickness $INCLUDES trace_optical_thickness.cpp $LDFLAGS $LIBS
