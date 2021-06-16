#!/bin/bash
set -eu -o pipefail


UMUS=$1
PHIS=$2
LIBRAD=$3
WORKDIR=$4
SAVEDIR=$5
PANDIR=$6

FLAG=0
END=$(($UMUS*$PHIS))

while [ $FLAG ! -eq $END ]
do
    FLAG=0
for umu in ${UMUS}
do
for phi in ${PHIS}
do
    cd $SAVEDIR
    JOBDIR=$SAVEDIR/job_mu${NM}/job_${umu}_${phi}
    cd $JOBDIR
    if [ ! -e mc.rad.spc.nc ]; then continue; fi
    if [ ! -e mc.flx.spc.nc ]; then continue; fi
    ((FLAG++))
done
done
done

#while :
#do
#    cd $PANDIR
#    if [ -e mc.rad.spc.nc ]; then break; fi
#done





