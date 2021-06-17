#!/bin/bash
set -eu -o pipefail


#directories
WORKDIR=$HOME/testgit/ma-project   ###dir of repo
LIBRAD=$WORKDIR/libRadtran
SAVEDIR=$WORKDIR/radiances/rad_test
PANDIR=$SAVEDIR/job_panorama_cam_below_clouds
TRACEDIR=$WORKDIR/trace_test
CLOUDDIR=$WORKDIR/radiances/clouds
ATM=$WORKDIR/radiances/stdatm/afglus_checkerboard.dat  

nm="4"     # 1 2 4 6 8 10 16 32
np="4"

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
    SZA=45 #def: 45
	PHI0=270 #def: 270


	#layers/levels
	ZLEV=$( tac $ATM | awk -F " " '{print $1}' | head  -n -2 | tr  "\n" " " )
	NLAY=$( cat $ATM | tail -n +4  | wc -l )
	ZLEVno0=$( echo $ZLEV | sed "s/^[^ ]* /-999 /" )

	#cloud data
	FNAME="wc3d_checkerboard_full.dat"  #def: "wc3D.dat" 
	
    ### only necessary if you make clouds via script(gen clouds script)########################################
    NX=2 #def 70              #NX NY define the grid for cloud boxes, e.g. NX=6 NY=6 means a 6x6 grid
	NY=2  #def 1
	DX=1      #def 0.1        #DX DY define the grid extent in km
	DY=1 #def 1
    CLDX="1 "         #$( seq 30 40 )        #CLDX CLDY define the gridboxes in which to put clouds (like coordinates;gridbox count starts from 1) 
	CLDY="1"
	CLDZ=$( seq 3 52 )                     #$( seq 1 50 )      #CLDZ layernumber in which to put clouds
	CR=10
	LWC="0.05"
    ########################################################################################

	#simulation data
	ALBEDO=0.2
	WAVELENGTH=500  #def: 500
    SAMPLEGRID="2 2 1 1" #def 70 1 0.1 1 #should probably be the same as "NX NY DX DY" ?

	#camera data
    fov_phi1="-45"
    fov_phi2="45"
    fov_theta1="45"
    fov_theta2="135"
    xloc=1  #in km
    yloc=1
    zloc=5
    xpixel=120
    ypixel=120
	mc_panorama_view="$fov_phi1 $fov_phi2 $fov_theta1 $fov_theta2"   #"phi1 phi2 theta1 theta2"; fov of camera in hor and vert direction; phi=0 -> south; theta=0 -> vertically down 
    mc_sensorposition="$(($xloc/1000)) $(($yloc/1000)) $((zloc/1000))"          #"x y z"; sensorpos in m
	mc_sample_grid="$xpixel $ypixel"              #"nphi ntheta"; how many samples are taken in each direction; (phi1-phi1)/nphi = camera resolution in phi direction 
    mc_backward="0 0 $((xpixel-1)) $((ypixel-1))"             #"phi_start theta_start phi_end theta_end"; ? should match sample grid, e.g. 90 90 -> 0 0 89 89
    cam_photons=1000
    PHOTONS=1000000




	#script execution
	#bash 03_gen_cloud.sh "$FNAME" "$ATM" "$NX" "$NY" "$DX" "$DY" "$CLDX" "$CLDY" "$CLDZ" "$ZLEV" "$NLAY" "$CR" "$LWC"
    bash 01_gen_radiance_jobs.sh "$UMUS" "$PHIS" "$SZA" "$PHI0" "$LIBRAD" "$WORKDIR" "$ATM" "$ALBEDO" "$WAVELENGTH" "$SAMPLEGRID" "$FNAME" "$ZLEVno0" "$PHOTONS" "$NUMU" "$SAVEDIR"
    bash 02_gen_panorama.sh "$UMUS" "$PHIS" "$SZA" "$PHI0" "$cam_photons" "$LIBRAD" "$WORKDIR" "$ATM" "$ALBEDO" "$WAVELENGTH" "$FNAME" "$mc_panorama_view" \
        "$mc_sensorposition" "$mc_sample_grid" "$mc_backward" "$PANDIR"


	#output used parameters
	bash op_parameters.sh "$UMUS" "$PHIS" "$SZA" "$PHI0" "$ZLEV" "$NLAY" "$SAMPLEGRID" "$LIBRAD" "$WORKDIR" "$ATM"
	
    #waiting for file completion/error checking
    bash checkforfiles.sh "$UMUS" "$PHIS" "$LIBRAD" "$WORKDIR" "$SAVEDIR" "$PANDIR" "$nm" "$np"

    #generate optical properties and flx file
    bash gen_opprop.sh "$UMUS" "$PHIS" "$LIBRAD" "$SAVEDIR"

	#combine radiances into netcdf
	radfn="rad_testcases/rad_cloud_above_cam__mu${NUMU%.*}_phi${NPHIS%.*}.nc"
    python3 mergerads.py -Nmu $NUMU -1 1 -Nphi $NPHIS -radfn $radfn

   
        
    #generate config
    outputfn="$TRACEDIR/output2d_checkerboard.nc"
    nsub=5
    bash gen_config.sh "$radfn" "$SAVEDIR/test.optical_properties.nc" "$SAVEDIR/mc.flx.spc.nc" "$outputfn" "$DX" "$DY" "$ALBEDO" "$xpixel" "$ypixel" "$fov_phi1" "$fov_phi2" 
                        "$fov_theta1" "$fov_theta2" "$xloc" "$yloc" "$zloc" "$nsub" "$TRACEDIR"
   
    $TRACEDIR/trace_optical_thickness
    ((ncview $outputfn)&)
    ((ncview $PANDIR/mc.rad.spc.nc)&)






done
done
