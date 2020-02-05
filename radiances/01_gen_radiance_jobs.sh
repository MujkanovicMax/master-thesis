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
ZLEV=${10}

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
zout -999 $ZLEV
rte_solver mystic
wc_file 3D $FNAME
wc_properties hu
mc_vroom on
mc_std
mc_sample_grid $SAMPLEGRID 
verbose
EOFJOB

cd $JOBDIR
if [ ! -e $FNAME ]; then  ln -s ../$FNAME; fi
$LIBRAD/bin/uvspec -f uvspec.inp 2>&1 > out.data
if [ ! -e 04_convertToNetCDF.sh ]; then ln -s ../04_convertToNetCDF.sh; fi
bash 04_convertToNetCDF.sh 
if [ ! -e convert_flx2nc.sh ]; then ln -s ../convert_flx2nc.sh; fi
bash convert_flx2nc.sh 

cd $WORKDIR

done
done

