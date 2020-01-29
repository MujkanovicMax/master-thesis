#!/bin/bash
set -eu -o pipefail

#directories
LIBRAD=$WORK/maproject/libRadtran
WORKDIR=$(pwd)
ATM=$WORKDIR/stdatm/afglus.dat

#angles
UMUS=$( sed -n "1 p" input_params.txt )       
PHIS=$( sed -n "2 p" input_params.txt )

for umu in ${UMUS}
do
for phi in ${PHIS}
do
    JOBDIR=job_${umu}_${phi}
    if [ -e $JOBDIR/mc.rad.spc.nc ]; then rm $JOBDIR/mc.rad.spc.nc; fi
    cd $JOBDIR
    bash 04_convertToNetCDF.sh
    cd $WORKDIR
done
done

    				








