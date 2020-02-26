#!/bin/bash
set -eu -o pipefail

UMUS=$1
PHIS=$2
SZA=$3
PHI0=$4
LIBRAD=$5
WORKDIR=$6
ATM=$7
ALBEDO=$8
WAVELENGTH=$9
SAMPLEGRID=${10}
FNAME=${11}
ZLEV=${12}

for umu in ${UMUS}
do
for phi in ${PHIS}
do
    JOBDIR=job_${umu}_${phi}
    if [ ! -e $JOBDIR ]; then mkdir $JOBDIR; fi
    cp $ATM $JOBDIR/
    cat > $JOBDIR/uvspec.inp << EOFJOB
data_files_path ${LIBRAD}/data
atmosphere_file ${ATM}
source solar ${LIBRAD}/data/solar_flux/atlas_plus_modtran
albedo $ALBEDO
wavelength $WAVELENGTH
umu $umu
phi $phi
sza $SZA
phi0 $PHI0
zout -999 1
rte_solver mystic
wc_file 3D $FNAME
no_scattering mol
no_absorption mol
wc_properties hu
mc_vroom on
mc_photons 10000
mc_std
mc_sample_grid $SAMPLEGRID 
verbose
EOFJOB

cd $JOBDIR
if [ ! -e $FNAME ]; then  ln -s ../$FNAME; fi
$LIBRAD/bin/uvspec -f uvspec.inp > out.data
if [ ! -e 04_convertToNetCDF.sh ]; then ln -s ../04_convertToNetCDF.sh; fi
bash 04_convertToNetCDF.sh 
if [ ! -e convert_flx2nc.sh ]; then ln -s ../convert_flx2nc.sh; fi
bash convert_flx2nc.sh 

cd $WORKDIR

done
done

