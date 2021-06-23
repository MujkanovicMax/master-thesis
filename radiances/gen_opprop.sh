#!/bin/bash
set -eu -o pipefail


UMUS=$1
PHIS=$2
FNAME=$3
LIBRAD=$4
SAVEDIR=$5
CLOUDDIR=$6
WORKDIR=$7
nmu=$8

for umu in ${UMUS}
do
for phi in ${PHIS}
do
    
    cd $SAVEDIR
    JOBDIR=$SAVEDIR/job_mu${nmu}/job_${umu}_${phi}
    if [ ! -e uvspec_opprop.inp ]; then cp $JOBDIR/uvspec.inp uvspec_opprop.inp; echo test_optical_properties >> uvspec_opprop.inp ;fi
    if [ ! -e $FNAME ]; then ln -s $CLOUDDIR/$FNAME; fi
    if [ ! -e test.optical_properties.nc ]; then $LIBRAD/bin/uvspec -f uvspec_opprop.inp; fi
    if [ ! -e mc.flx.spc.nc ]; then mv $JOBDIR/mc.flx.spc.nc ./; fi
    break

done
done


#python3 $WORKDIR/radiances/check_opprop.py

    


