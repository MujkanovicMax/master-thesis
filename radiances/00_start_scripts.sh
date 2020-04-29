#!/bin/bash
set -eu -o pipefail


#directories
LIBRAD=$WORK/maproject/libRadtran
WORKDIR=$(pwd)
ATM=$WORKDIR/stdatm/afglus_mgl.dat

nm="32"     # 1 2 4 6 8 10 16 32
np="2 4 8 16 32 64"

for i in $nm
do
for j in $np
do

	#angles
	NUMU=$i                             #NUMU defines number of sample angles chosen from gaussian quadrature. Disort uses 16 streams per default (for irradiance) 
	python3 angles.py -Nmu $NUMU -1 1 > numus.txt 
	UMUS=$( head -n +1 numus.txt )    #when not using gaussian quadrature angles use something like this: UMUS="-0.861136" #-0.339981 0.339981 0.861136   
	NPHIS=$j
	python3 angles.py -Nphi $NPHIS > nphis.txt 
	PHIS=$( head -n +1 nphis.txt )
    SZA=45
	PHI0=270


	#layers/levels
	ZLEV=$( tac $ATM | awk -F " " '{print $1}' | head  -n -2 | tr  "\n" " " )
	NLAY=$( cat $ATM | tail -n +4  | wc -l )
	ZLEVno0=$( echo $ZLEV | sed "s/^[^ ]* /-999 /" )

	#cloud data
	FNAME="wc3D.dat"   
	NX=70              #NX NY define the grid for cloud boxes, e.g. NX=6 NY=6 means a 6x6 grid
	NY=1
	DX=0.1              #DX DY define the grid extent in km
	DY=1
	CLDX=$( seq 30 40 )        #CLDX CLDY define the gridboxes in which to put clouds (like coordinates;gridbox count starts from 1) 
	CLDY="1"
	CLDZ=$( seq 1 50 )      #CLDZ layernumber in which to put clouds
	CR=10
	LWC="0.05"


	#simulation data
	ALBEDO=0.2
	WAVELENGTH=500
	SAMPLEGRID="70 1 0.1 1" #should probably be the same as "NX NY DX DY" ?

	#camera data
	mc_panorama_view="-45 45 90 90"   #"phi1 phi2 theta1 theta2"; fov of camera in hor and vert direction; phi=0 -> south; theta=0 -> vertically down 
	mc_sensorposition="3000 10 2000"          #"x y z"; sensorpos in m
	mc_sample_grid="90 1"              #"nphi ntheta"; how many samples are taken in each direction; (phi1-phi1)/nphi = camera resolution in phi direction 
	mc_backward="0 0 89 0"             #"phi_start theta_start phi_end theta_end"; ? should match sample grid, e.g. 90 90 -> 0 0 89 89
    PHOTONS=1000000




	#script execution
	#bash 03_gen_cloud.sh "$FNAME" "$ATM" "$NX" "$NY" "$DX" "$DY" "$CLDX" "$CLDY" "$CLDZ" "$ZLEV" "$NLAY" "$CR" "$LWC"
	#echo Cloud gernerated
	#bash 01_gen_radiance_jobs.sh "$UMUS" "$PHIS" "$SZA" "$PHI0" "$LIBRAD" "$WORKDIR" "$ATM" "$ALBEDO" "$WAVELENGTH" "$SAMPLEGRID" "$FNAME" "$ZLEVno0" "$PHOTONS" "$NUMU"
	#echo Radiances generated
	#bash 02_gen_panorama.sh "$UMUS" "$PHIS" "$SZA" "$PHI0" "$PHOTONS" "$LIBRAD" "$WORKDIR" "$ATM" "$ALBEDO" "$WAVELENGTH" "$FNAME" "$mc_panorama_view" "$mc_sensorposition" "$mc_sample_grid" "$mc_backward"
	#echo Panorama generated
	
	#output used parameters
	#bash op_parameters.sh "$UMUS" "$PHIS" "$SZA" "$PHI0" "$ZLEV" "$NLAY" "$SAMPLEGRID" "$LIBRAD" "$WORKDIR" "$ATM"
	

	#combine radiances into netcdf
	python3 mergerads.py -Nmu $NUMU -1 1 -Nphi $NPHIS -radfn radiances_mu${NUMU%.*}_phi${NPHIS%.*}.nc



done
done
