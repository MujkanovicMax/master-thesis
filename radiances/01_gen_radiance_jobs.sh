#!/bin/bash
set -eu -o pipefail

UMUS=$1
PHIS=$2
LIBRAD=$3
WORKDIR=$4
ATM=$5
ALBEDO=$6
WAVELENGTH=$7
SAMPLEGRID=$8
FNAME=$9

for umu in ${UMUS}
do
for phi in ${PHIS}
do
    JOBDIR=job_${umu}_${phi}
    if [ ! -e $JOBDIR ]; then mkdir $JOBDIR; fi
    cp $LIBRAD/data/atmmod/afglus.dat $JOBDIR/
    cat > $JOBDIR/uvspec.inp << EOFJOB
data_files_path ${LIBRAD}/data
atmosphere_file ${ATM}
source solar ${LIBRAD}/data/solar_flux/atlas_plus_modtran
albedo $ALBEDO
wavelength $WAVELENGTH
umu $umu
phi $phi
wc_file 3D $FNAME
mc_sample_grid $SAMPLEGRID 
verbose
EOFJOB

cd $JOBDIR
$LIBRAD/bin/uvspec -f uvspec.inp > out.data
cd $WORKDIR


done
done

