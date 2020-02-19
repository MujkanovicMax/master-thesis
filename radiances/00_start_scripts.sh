#!/bin/bash
set -eu -o pipefail


#directories
LIBRAD=$WORK/maproject/libRadtran
WORKDIR=$(pwd)
ATM=$WORKDIR/stdatm/afglus.dat

nangles="2"     # 1 2 4 6 8 10 16 32

for i in $nangles
do
	#angles
	NUMU=$i                              #NUMU defines number of sample angles chosen from gaussian quadrature. Disort uses 16 streams per default (for irradiance) 
	bash 05_gauss_nodes.sh "$NUMU" numus.txt "-1" "1" 
	UMUS=$( head -n +1 numus.txt )    #when not using gaussian quadrature angles use something like this: UMUS="-0.861136" #-0.339981 0.339981 0.861136   
	NPHIS=$i
	bash 05_gauss_nodes.sh "$NPHIS" nphis.txt "0" "6.283185307"
	bash convertToAngle.sh nphis.txt
	PHIS=$( head cangles.txt )

	#layers/levels
	ZLEV=$( tac $ATM | awk -F " " '{print $1}' | head  -n -2 | tr  "\n" " " )
	NLAY=$( cat $ATM | tail -n +4  | wc -l )
	ZLEVno0=$( tac $ATM | awk -F " " '{print $1}' | head  -n -2 | tail -n -49  | tr  "\n" " " )

	#cloud data
	FNAME="wc3D.dat"   
	NX=3              #NX NY define the grid for cloud boxes, e.g. NX=6 NY=6 means a 6x6 grid
	NY=1
	DX=1              #DX DY define the grid extent in km
	DY=1
	CLDX="2"        #CLDX CLDY define the gridboxes in which to put clouds (like coordinates;gridbox count starts from 1) 
	CLDY="1"
	CLDZ="1"      #CLDZ layernumber in which to put clouds
	CR=10
	LWC="1e-30"


	#simulation data
	ALBEDO=0.2
	WAVELENGTH=500
	SAMPLEGRID="30 1 0.1 1" #should probably be the same as "NX NY DX DY" ?

	#camera data
	mc_panorama_view="-45 45 90 90"   #"phi1 phi2 theta1 theta2"; fov of camera in hor and vert direction; phi=0 -> south; theta=0 -> vertically down 
	mc_sensorposition="1500 10 2000"          #"x y z"; sensorpos in m
	mc_sample_grid="90 1"              #"nphi ntheta"; how many samples are taken in each direction; (phi1-phi1)/nphi = camera resolution in phi direction 
	mc_backward="0 0 89 0"             #"phi_start theta_start phi_end theta_end"; ? should match sample grid, e.g. 90 90 -> 0 0 89 89





	#script execution
	bash 03_gen_cloud.sh "$FNAME" "$ATM" "$NX" "$NY" "$DX" "$DY" "$CLDX" "$CLDY" "$CLDZ" "$ZLEV" "$NLAY" "$CR" "$LWC"
	echo Cloud gernerated
	#bash 01_gen_radiance_jobs.sh "$UMUS" "$PHIS" "$LIBRAD" "$WORKDIR" "$ATM" "$ALBEDO" "$WAVELENGTH" "$SAMPLEGRID" "$FNAME" "$ZLEVno0"
	#echo Radiances generated
	bash 02_gen_panorama.sh "$UMUS" "$PHIS" "$LIBRAD" "$WORKDIR" "$ATM" "$ALBEDO" "$WAVELENGTH" "$FNAME" "$mc_panorama_view" "$mc_sensorposition" "$mc_sample_grid" "$mc_backward"
	echo Panorama generated

	#output used parameters
	bash op_parameters.sh "$UMUS" "$PHIS" "$ZLEV" "$NLAY" "$LIBRAD" "$WORKDIR" "$ATM"
	
	

done


