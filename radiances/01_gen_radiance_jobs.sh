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
PHOTONS=${13}
NM=${14}
SAVEDIR=${15}

for umu in ${UMUS}
do
for phi in ${PHIS}
do
    JOBDIR=$SAVEDIR/job_mu${NM}/job_${umu}_${phi}
    echo $JOBDIR
    if [ ! -e $JOBDIR ]; then mkdir -p $JOBDIR; fi
    cp $ATM $JOBDIR/
    cat > $JOBDIR/uvspec.inp << EOFJOB
data_files_path ${LIBRAD}/data
atmosphere_file ${ATM}
source solar ${LIBRAD}/data/solar_flux/kurudz_1.0nm.dat
albedo $ALBEDO
wavelength $WAVELENGTH
umu $umu
phi $phi
sza $SZA
phi0 $PHI0
zout $ZLEV
rte_solver mystic
wc_file 3D $FNAME
no_scattering mol
no_absorption mol
wc_properties hu
mc_vroom on
mc_photons $PHOTONS
mc_std
mc_sample_grid $SAMPLEGRID 
verbose
EOFJOB

cat > $JOBDIR/slurm.job << EOF
#!/bin/bash -l
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=cluster,met-ws
#SBATCH --mem=25G
#SBATCH --time=8:00:00
#SBATCH --output=$JOBDIR/log.%j.out
#SBATCH --error=$JOBDIR/log.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=Mujkanovic.Max@physik.uni-muenchen.de
cd $JOBDIR

if [ ! -e $FNAME ]; then  ln -s $WORKDIR/$FNAME; fi
if [ ! -e mc.flx.spc ]; then $LIBRAD/bin/uvspec -f uvspec.inp > out.data || exit 1; fi
if [ ! -e 04_convertToNetCDF.sh ]; then ln -s $WORKDIR/04_convertToNetCDF.sh; fi
bash 04_convertToNetCDF.sh 
if [ ! -e convert_flx2nc.sh ]; then ln -s $WORKDIR/convert_flx2nc.sh; fi
bash convert_flx2nc.sh
rm -f mc*{.rad,.rad.std,.flx,.flx.std}
EOF

sbatch $JOBDIR/slurm.job


done
done


