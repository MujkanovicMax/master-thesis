#!/bin/bash 

set -eu -o pipefail

function gen_cld {
        FNAME=$1
	ATM=$2
        NX=$3
        NY=$4
        DX=$5
        DY=$6
        CLDX=$7
        CLDY=$8
        CLDZ=$9
        ZLEV=${10}
        NLAY=${11}

        cat > $FNAME << EOF
$NX $NY $NLAY 3
$DX $DY $ZLEV
EOF
for k in $(seq $NLAY)
do
        echo 1 1 $k 1e-30 10 >> $FNAME
done
echo $CLDX
echo $CLDY
echo $CLDZ
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


gen_cld  "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}" "${11}"

#example call
#gen_cld wc3D.dat 6 6 1 1 "3 4" "3 4" "10 11"

