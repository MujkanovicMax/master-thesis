#!/bin/bash
set -eu -o pipefail


radfn=$1
opfn=$2
flxfn=$3
outfn=$4
dx=$5
dy=$6
albedo=$7
xpixel=$8
ypixel=$9
fov_phi1=${10}
fov_phi2=${11}
fov_theta1=${12}
fov_theta2=${13}
xloc=${14}
yloc=${15}
zloc=${16}
nsub=${17}
TRACEDIR=${18}

cd $TRACEDIR
cat > config.txt << EOFJOB
rad_filename = $radfn
flx_filename = $flxfn
op_filename = $opfn
output_filename = $outfn
dx = $dx
dy = $dy
albedo = $albedo
xpixel = $xpixel
ypixel = $ypixel
fov_phi1 = $fov_phi1
fov_phi2 = $fov_phi2
fov_theta1 = $fov_theta1
fov_theta2 = $fov_theta2
xloc = $xloc
yloc = $yloc
zloc = $zloc
nsub = $nsub
EOFJOB


