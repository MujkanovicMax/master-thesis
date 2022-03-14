#!/bin/bash

dir="/home/m/Mujkanovic.Max/ma/trace_test/weightingtest"
#divs="50 25 10 5 2 1"
#mus="32 16 8 4 2"
#phis="4 8"
#nsubs="0 1 2 3 4 5 6 7 8 9 10 100 200 300" #11 12 13 14 15 16 17 18 19 20"
#arr=($nsubs)
#i=0
#inputs="rad_checkerboard_32x32_meaned_to_32x32.nc" #rad_checkerboard_32x32_meaned_to_16x16.nc rad_checkerboard_32x32_meaned_to_2x1.nc"

#for inp in $inputs
#do
#for nsub in $nsubs
#do
#    
#    inpfn="${dir}/rad_checkerboard_32x32_meaned_to_2x1.nc"
#    outfn="/home/m/Mujkanovic.Max/ma/trace_test/newtests/nsub/output_checkerboard_2x1_nsub${nsub}.nc"
#    
#    sed -i 's~^rad_filename.*$~rad_filename = '"$inpfn"'~' config.txt
#    sed -i 's~^output_filename.*$~output_filename = '"$outfn"'~' config.txt
#    sed -i 's~^nsub.*$~nsub = '"$nsub"'~' config.txt
#
#    ./trace_optical_thickness
#
#
#
#done
#done
modes="0 1 2 3 4 5"
for mode in $modes
do
    
    inpfn="${dir}/rad_checkerboard_32x32_meaned_to_32x32_zdiv50.nc"
    outfn="${dir}/output/output_rad_checkerboard_32x32_meaned_to_32x32_zdiv50_mode${mode}.nc"
    flxfn="${dir}/flx_checkerboard_32x32_zdiv50.nc"
    opfn="${dir}/op_checkerboard_32x32_zdiv50.nc"
    
    sed -i 's~^rad_filename.*$~rad_filename = '"$inpfn"'~' config.txt
    sed -i 's~^output_filename.*$~output_filename = '"$outfn"'~' config.txt
    sed -i 's~^mode.*$~mode = '"$mode"'~' config.txt
    sed -i 's~^flx_filename.*$~flx_filename = '"$flxfn"'~' config.txt
    sed -i 's~^op_filename.*$~op_filename = '"$opfn"'~' config.txt

    ./trace_optical_thickness
done


#radlist="rad_checkerboard_32x32_meaned_to_32x32_zdiv1.nc rad_checkerboard_32x32_meaned_to_2x1_zdiv1.nc irr_checkerboard_32x32_zdiv1.nc" 
#oplist="op_checkerboard_32x32_zdiv1.nc"
#flxlist="flx_checkerboard_32x32_zdiv1.nc"
#nsubs="0 12 12"
#arrnsubs=($nsubs)
#i=0
#for rad in $radlist
#do
#
#    outfn="problemtest/output2d_${rad}"
#    radfn="checkerboard_perf/${rad}"
#    opfn="checkerboard_perf/${oplist}"
#    flxfn="checkerboard_perf/${flxlist}"
#    nsub=${arrnsubs[i]}
#
#    sed -i 's~^rad_filename.*$~rad_filename = '"$radfn"'~' config.txt
#    sed -i 's~^flx_filename.*$~flx_filename = '"$flxfn"'~' config.txt
#    sed -i 's~^op_filename.*$~op_filename = '"$opfn"'~' config.txt
#    sed -i 's~^output_filename.*$~output_filename = '"$outfn"'~' config.txt
#    sed -i 's~^nsub.*$~nsub = '"$nsub"'~' config.txt
#
#    ./trace_optical_thickness
#    ((i=i+1))
#
#done

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


