#!/bin/bash
set -eu -o pipefail


UMUS=$1
PHIS=$2
LIBRAD=$3
WORKDIR=$4
SAVEDIR=$5
PANDIR=$6
nmu=$7
nphi=$8
rep=$9

FLAG=0
END=$(($nmu*$nphi))

while [ $FLAG != $END ]
do
    FLAG=0

for umu in ${UMUS}
do
for phi in ${PHIS}
do
    cd $SAVEDIR
    JOBDIR=$SAVEDIR/job_mu${nmu}/job_${umu}_${phi}
    cd $JOBDIR
    if [ ! -e mc.rad.spc.nc ]; then continue; fi
    if [ ! -e mc.flx.spc.nc ]; then continue; fi
    ((FLAG=FLAG+1))
done
done
done


FLAG=0
while [ $FLAG != $REP ]
do
    FLAG=0
for i in $( seq 0 $rep )
do
    cd $PANDIR/panorama_$i
    if [ ! -e mc.rad.spc.nc ]; then continue; fi
    ((FLAG=FLAG+1))
done
done






