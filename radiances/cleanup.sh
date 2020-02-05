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
    echo $JOBDIR
    if [ -e $JOBDIR/mc.flx.spc.nc ]; then rm $JOBDIR/mc.flx.spc.nc; fi
    cd $JOBDIR
    if [ ! -e convert_flx2nc.sh ]; then ln -s ../convert_flx2nc.sh; fi
    bash convert_flx2nc.sh
    rm convert_flx2nc.sh
    cd $WORKDIR
done
done

    				








