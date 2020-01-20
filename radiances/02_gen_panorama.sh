#!/bin/bash
set -eu -o pipefail

UMUS=$1
PHIS=$2
LIBRAD=$3
WORKDIR=$4
ATM=$5
ALBEDO=$6
WAVELENGTH=$7
FNAME=$8
mc_panorama_view=$9
mc_sensorposition=${10}
mc_sample_grid=${11}
mc_backward=${12}

for umu in ${UMUS}
do
for phi in ${PHIS}
do
    JOBDIR=job_${umu}_${phi}
    if [ ! -e $JOBDIR ]; then mkdir $JOBDIR; fi
    cp $LIBRAD/data/atmmod/afglus.dat $JOBDIR/
    cat > $JOBDIR/uvspec_panorama.inp << EOFJOB
data_files_path ${LIBRAD}/data
atmosphere_file ${ATM}
source solar ${LIBRAD}/data/solar_flux/atlas_plus_modtran
albedo $ALBEDO
wavelength $WAVELENGTH
umu $umu
phi $phi
wc_file 3D $FNAME
mc_panorama_view $mc_panorama_view 
mc_sensorposition $mc_sensorposition
mc_sample_grid $mc_sample_grid 
mc_backward $mc_backward 
verbose
EOFJOB

cd $JOBDIR
$LIBRAD/bin/uvspec -f uvspec_panorama.inp > out.data
cd $WORKDIR


done
done












