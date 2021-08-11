#!/bin/bash

dir="checkerboard_perf"
divs="50 25 10 5 2 1"
mus="32 16 8 4 2"
#phis="4 8"
nsubs="2 5 8 10"
arr=($nsubs)
i=0



for nsub in $nsubs
do
    

    outfn="nsubproblem/output2d_plane_nsub${nsub}.nc"
    
    sed -i 's~^output_filename.*$~output_filename = '"$outfn"'~' config.txt
    sed -i 's~^nsub.*$~nsub = '"$nsub"'~' config.txt

    ./trace_optical_thickness



done


#for mu in $mus
#do
#    nsub="${arr[i]}"
#    ((i=i+1))
#
#for div in $divs
#do
#    
#    radfn="${dir}/rad_checkerboard_32x32_meaned_to_${mu}x${mu}_zdiv${div}.nc"
#    flxfn="${dir}/flx_checkerboard_32x32_zdiv${div}.nc"
#    opfn="${dir}/op_checkerboard_32x32_zdiv${div}.nc"
#    outfn="${dir}/output2d_rad_checkerboard_32x32_meaned_to_${mu}x${mu}_zdiv${div}nc"
#
#    sed -i 's~^rad_filename.*$~rad_filename = '"$radfn"'~' config.txt
#    sed -i 's~^flx_filename.*$~flx_filename = '"$flxfn"'~' config.txt
#    sed -i 's~^op_filename.*$~op_filename = '"$opfn"'~' config.txt
#    sed -i 's~^output_filename.*$~output_filename = '"$outfn"'~' config.txt
#    sed -i 's~^nsub.*$~nsub = '"$nsub"'~' config.txt
#
#    ./trace_optical_thickness
#
#    ((i=i+1))
#
#done
#i=0
#done
#
#
#for div in $divs
#do
#    
#    radfn="${dir}/rad_checkerboard_32x32_meaned_to_2x1_zdiv${div}.nc"
#    flxfn="${dir}/flx_checkerboard_32x32_zdiv${div}.nc"
#    opfn="${dir}/op_checkerboard_32x32_zdiv${div}.nc"
#    outfn="${dir}/output2d_rad_checkerboard_32x32_meaned_to_2x1_zdiv${div}nc"
#
#    nsub="25"
#
#    sed -i 's~^rad_filename.*$~rad_filename = '"$radfn"'~' config.txt
#    sed -i 's~^flx_filename.*$~flx_filename = '"$flxfn"'~' config.txt
#    sed -i 's~^op_filename.*$~op_filename = '"$opfn"'~' config.txt
#    sed -i 's~^output_filename.*$~output_filename = '"$outfn"'~' config.txt
#    sed -i 's~^nsub.*$~nsub = '"$nsub"'~' config.txt
#
#    ./trace_optical_thickness
#
#
#
#done
#
#for div in $divs
#do
#   
#    radfn="${dir}/irr_checkerboard_32x32_zdiv${div}.nc"
#    flxfn="${dir}/flx_checkerboard_32x32_zdiv${div}.nc"
#    opfn="${dir}/op_checkerboard_32x32_zdiv${div}.nc"
#    outfn="${dir}/output2d_irr_checkerboard_32x32_meaned_to_2x1_zdiv${div}.nc"
#
#    nsub="25"
#
#    sed -i 's~^rad_filename.*$~rad_filename = '"$radfn"'~' config.txt
#    sed -i 's~^flx_filename.*$~flx_filename = '"$flxfn"'~' config.txt
#    sed -i 's~^op_filename.*$~op_filename = '"$opfn"'~' config.txt
#    sed -i 's~^output_filename.*$~output_filename = '"$outfn"'~' config.txt
#    sed -i 's~^nsub.*$~nsub = '"$nsub"'~' config.txt
#
#    ./trace_optical_thickness
#
#
#
#done


