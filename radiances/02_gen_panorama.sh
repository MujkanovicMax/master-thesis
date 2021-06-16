#!/bin/bash
set -eu -o pipefail

UMUS=$1
PHIS=$2
SZA=$3
PHI0=$4
PHOTONS=$5
LIBRAD=$6
WORKDIR=$7
ATM=$8
ALBEDO=$9
WAVELENGTH=${10}
FNAME=${11}
mc_panorama_view=${12}
mc_sensorposition=${13}
mc_sample_grid=${14}
mc_backward=${15}
SAVEDIR=${16}


JOBDIR=$SAVEDIR
if [ ! -e $JOBDIR ]; then mkdir -p $JOBDIR; fi
echo $JOBDIR
cp $ATM $JOBDIR/
cat > $JOBDIR/uvspec_panorama.inp << EOFJOB 
data_files_path ${LIBRAD}/data
atmosphere_file ${ATM}
source solar ${LIBRAD}/data/solar_flux/atlas_plus_modtran
albedo $ALBEDO
sza $SZA  #45
wavelength $WAVELENGTH
mol_abs_param reptran coarse
mc_panorama_alignment mu
umu 0
phi 0
phi0 $PHI0 #270
rte_solver mystic
mc_photons $PHOTONS #10000 #1000
mc_minphotons 1
wc_properties hu
wc_file 3D $FNAME
mc_panorama_view $mc_panorama_view 
mc_sensorposition $mc_sensorposition
mc_sample_grid $mc_sample_grid 
mc_backward $mc_backward 
mc_surface_reflectalways
#mc_maxscatters 1
mc_vroom on
no_scattering mol
no_absorption mol

#test_optical_properties

EOFJOB

cat > $JOBDIR/slurm.job << EOF
#!/bin/bash -l
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=cluster,met-ws
#SBATCH --mem=25G
#SBATCH --time=36:00:00
#SBATCH --output=$JOBDIR/log.%j.out
#SBATCH --error=$JOBDIR/log.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=Mujkanovic.Max@physik.uni-muenchen.de


cd $JOBDIR
if [ ! -e $FNAME ]; then  ln -s ../$FNAME; fi 
#if [ ! -e out.data ]; then
$LIBRAD/bin/uvspec -f uvspec_panorama.inp > out_panorama.data
#fi
if [ ! -e 04_convertToNetCDF.sh ]; then ln -s /home/m/Mujkanovic.Max/ma/radiances/04_convertToNetCDF.sh; fi
bash 04_convertToNetCDF.sh
cd $WORKDIR
EOF

sbatch $JOBDIR/slurm.job













