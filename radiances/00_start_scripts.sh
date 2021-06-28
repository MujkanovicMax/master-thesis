#!/bin/bash
set -eu -o pipefail


#directories
WORKDIR=$HOME/testgit/ma-project   ###dir of repo
LIBRAD=$WORKDIR/libRadtran
SAVEDIR=$WORKDIR/radiances/rad_test
PANDIR=$SAVEDIR/job_panorama
TRACEDIR=$WORKDIR/trace_test
CLOUDDIR=$WORKDIR/radiances/clouds
ATM=$WORKDIR/radiances/stdatm/afglus_checkerboard.dat  

#number of angles for radiancesi (can be lists)
nm="4"                                          # 1 2 4 6 8 10 16 32
np="4"

#dirs and filenames
PANFN="test_panorama.nc" 
cloudfn="wc3d_checkerboard_full.dat"
radfn="rad_test.nc"
outputfn="output2d_test.nc"

#simulation data
ALBEDO=0.2
WAVELENGTH=500
SAMPLEGRID="2 2 1 1"            #"nx ny dx dy"; nx ny number of boxes in x and y; dx dy box lengths in km

#angles
SZA=45                          #solar zenith angle
PHI0=270                        #where the sun is shining from (0 north, 90 west, 180 south, 270 east)
nsub=5                          #number of subdivisions for phase function calculation

#camera data
fov_phi1="-45"                  #
fov_phi2="45"                   #   angles for camera (rotation and fov)
fov_theta1="45"                 #
fov_theta2="135"                #
xpixel=120                      #   number of pixel in x and y
ypixel=120                      #

#positions
xloc=1                          #camera position in km
yloc=1
zloc=5

#photon numbers
PHOTONS=1000000                 # ???
cam_photons=100                 #number of photons for mystic panorama cam
REP=100                         #number of repetitions for mystoc panorama calculations


#main loop execution
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

	#layers/levels
	ZLEV=$( tac $ATM | awk -F " " '{print $1}' | head  -n -2 | tr  "\n" " " )
	NLAY=$( cat $ATM | tail -n +4  | wc -l )
	ZLEVno0=$( echo $ZLEV | sed "s/^[^ ]* /-999 /" )

	#cloud data
    ### only necessary if you make clouds via script(gen clouds script)########################################
    NX=2 #def 70              #NX NY define the grid for cloud boxes, e.g. NX=6 NY=6 means a 6x6 grid
	NY=2  #def 1
	DX=1      #def 0.1        #DX DY define the grid extent in km
	DY=1 #def 1
    CLDX="1 "         #$( seq 30 40 )        #CLDX CLDY define the gridboxes in which to put clouds (like coordinates;gridbox count starts from 1) 
	CLDY="1"
	CLDZ=$( seq 2 3 )                     #$( seq 1 50 )      #CLDZ layernumber in which to put clouds
	CR=10
	LWC="0.5"
    ########################################################################################

	#panorama data
	mc_panorama_view="$fov_phi1 $fov_phi2 $fov_theta1 $fov_theta2"   #"phi1 phi2 theta1 theta2"; fov of camera in hor and vert direction; phi=0 -> south; theta=0 -> vertically down 
    mc_sensorposition="$(($xloc*1000)) $(($yloc*1000)) $((zloc*1000))"          #"x y z"; sensorpos in m
	mc_sample_grid="$xpixel $ypixel"              #"nphi ntheta"; how many samples are taken in each direction; (phi1-phi1)/nphi = camera resolution in phi direction 
    mc_backward="0 0 $((xpixel-1)) $((ypixel-1))"             #"phi_start theta_start phi_end theta_end"; ? should match sample grid, e.g. 90 90 -> 0 0 89 89


	#script execution
	#bash 03_gen_cloud.sh "$cloudfn" "$ATM" "$NX" "$NY" "$DX" "$DY" "$CLDX" "$CLDY" "$CLDZ" "$ZLEV" "$NLAY" "$CR" "$LWC" "$CLOUDDIR"
    bash 01_gen_radiance_jobs.sh "$UMUS" "$PHIS" "$SZA" "$PHI0" "$LIBRAD" "$WORKDIR" "$ATM" "$ALBEDO" "$WAVELENGTH" "$SAMPLEGRID" "$cloudfn" "$ZLEVno0" "$PHOTONS" "$NUMU" "$SAVEDIR" "$CLOUDDIR"
    bash 02_gen_panorama.sh "$UMUS" "$PHIS" "$SZA" "$PHI0" "$cam_photons" "$LIBRAD" "$WORKDIR" "$ATM" "$ALBEDO" "$WAVELENGTH" "$cloudfn" "$mc_panorama_view" \
        "$mc_sensorposition" "$mc_sample_grid" "$mc_backward" "$PANDIR" "$CLOUDDIR" "$REP"

	#output used parameters
	bash op_parameters.sh "$UMUS" "$PHIS" "$SZA" "$PHI0" "$ZLEV" "$NLAY" "$SAMPLEGRID" "$LIBRAD" "$WORKDIR" "$ATM"
    
    #waiting for file completion/error checking
    bash checkforrad.sh "$UMUS" "$PHIS" "$LIBRAD" "$WORKDIR" "$SAVEDIR" "$PANDIR" "$i" "$j" 
    bash checkforpano.sh "$PANDIR" "$REP" 

    #generate optical properties and flx file
    bash gen_opprop.sh "$UMUS" "$PHIS" "$cloudfn" "$LIBRAD" "$SAVEDIR" "$CLOUDDIR" "$WORKDIR "$WORKDIR"" "$i"
    
	#combine radiances into netcdf
    python3 mergerads.py -Nmu $NUMU -1 1 -Nphi $NPHIS -radfn $radfn -loc $SAVEDIR -panfn $PANFN -rep $REP -locpan $PANDIR 
        
    #generate config
    bash gen_config.sh "$radfn" "$SAVEDIR/test.optical_properties.nc" "$SAVEDIR/mc.flx.spc.nc" "$outputfn" "$DX" "$DY" "$ALBEDO" "$xpixel" "$ypixel" "$fov_phi1" "$fov_phi2" \
        "$fov_theta1" "$fov_theta2" "$xloc" "$yloc" "$zloc" "$nsub" "$TRACEDIR"
    
    cd $TRACEDIR
    ./trace_optical_thickness
    cd -
    ((ncview $outputfn)&)
    #((ncview $PANDIR/mc.rad.spc.nc)&)






done
done
