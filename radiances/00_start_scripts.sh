#!/bin/bash
set -eu -o pipefail

#directories
LIBRAD=$WORK/libRadtran
WORKDIR=$(pwd)
ATM=$WORKDIR/stdatm/afglus.dat

#angles
UMUS=".1 .5 .9 1"
PHIS="0 180"

#cloud data
FNAME=wc3D.dat   
NX=6              #NX NY define the grid for cloud boxes, e.g. NX=6 NY=6 means a 6x6 grid
NY=6
DX=1              #DX DY define the grid extent in km
DY=1
CLDX="3 4"        #CLDX CLDY define the gridboxes in which to put clouds (like coordinates;gridbox count starts from 1) 
CLDY="3 4"
CLDZ="10 11"      #CLDZ layernumber in which to put clouds

#simulation data
ALBEDO=0.2
WAVELENGTH=550
SAMPLEGRID="6 6 1 1" #should probably be the same as "NX NY DX DY" ?

#camera data
mc_panorama_view="315 225 350 80"   #"phi1 phi2 theta1 theta2"; fov of camera in hor and vert direction; phi=0 -> south; theta=0 -> vertically down 
mc_sensorposition="0 0 13"          #"x y z"; sensorpos in km
mc_sample_grid="90 90"              #"nphi ntheta"; how many samples are taken in each direction; (phi1-phi1)/nphi = camera resolution in phi direction 
mc_backward="0 0 89 89"             #"phi_start theta_start phi_end theta_end"; ? should match sample grid, e.g. 90 90 -> 0 0 89 89


bash 03_gen_cloud.sh "$FNAME" "$ATM" "$NX" "$NY" "$DX" "$DY" "$CLDX" "$CLDY" "$CLDZ"
bash 01_gen_radiance_jobs.sh "$UMUS" "$PHIS" "$LIBRAD" "$WORKDIR" "$ATM" "$ALBEDO" "$WAVELENGTH" "$SAMPLEGRID" "$FNAME"
bash 02_gen_panorama.sh "$UMUS" "$PHIS" "$LIBRAD" "$WORKDIR" "$ATM" "$ALBEDO" "$WAVELENGTH" "$FNAME" "$mc_panorama_view" "$mc_sensorposition" "$mc_sample_grid" "$mc_backward"

