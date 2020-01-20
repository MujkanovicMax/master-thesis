#!/bin/bash

set -eu -o pipefail

LIBRAD=$WORK/libRadtran
WORKDIR=$(pwd)
CLOUDDIR=$WORK/clouds
ATM=$WORKDIR/stdatm/afglus.dat
UMUS=".1 .5 .9 1"
PHIS="0 180"

function gen_cld {
        FNAME=$1
        NX=$2
        NY=$3
        DX=$4
        DY=$5
        CLDX=$6
        CLDY=$7
        CLDZ=$8
        ZLEV=$( tac $ATM | awk -F " " '{print $1}' | head  -n -2 | tr  "\n" " " )
        NLAY=$( cat $ATM | tail -n +4  | wc -l )

        cat > $FNAME << EOF
$NX $NY $NLAY 3
$DX $DY $ZLEV
EOF
for k in $(seq $NLAY)
do
        echo 1 1 $k 1e-30 10 >> $FNAME
done
for i in $CLDX
do
for j in $CLDY
do
for k in $CLDZ
do
        echo $i $j $k 0.1 10 >> $FNAME
done
done
done
}


