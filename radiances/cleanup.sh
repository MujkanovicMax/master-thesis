#!/bin/bash
set -eu -o pipefail

#directories
LIBRAD=$WORK/maproject/libRadtran
WORKDIR=$(pwd)
ATM=$WORKDIR/stdatm/afglus.dat

#angles
NUMU=16                              #NUMU defines number of sample angles chosen from gaussian quadrature. Disort uses 16 streams per default (for irradiance) 
bash 05_gauss_nodes.sh "$NUMU"
UMUS=$( head -n +1 gauss_nodes.txt )    #when not using gaussian quadrature angles use something like this: UMUS="-0.861136" #-0.339981 0.339981 0.861136   
PHIS=$( seq 0 36 360 )

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

    				








