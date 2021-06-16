#!/bin/bash
set -eu -o pipefail


UMUS=$1
PHIS=$2
LIBRAD=$3
SAVEDIR=$4


for umu in ${UMUS}
do
for phi in ${PHIS}
do
    
    cd $SAVEDIR
    JOBDIR=$SAVEDIR/job_mu${NM}/job_${umu}_${phi}
    if [ ! -e uvspec_opprop.inp ]; then cp $JOBDIR/uvspec.inp uvspec_opprop.inp; echo test_optical_properties >> uvspec_opprop.inp ;fi
    if [ ! -e test.optical_properties.nc ]; then $LIBRAD/bin/uvspec -f uvspec_opprop.inp; fi
    if [ ! -e mc.flx.spc.nc ]; then mv $JOBDIR/mc.flx.spc.nc ./; fi
    break

done
done


    

