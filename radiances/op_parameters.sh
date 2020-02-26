#!/bin/bash
set -eu -o pipefail

cat > input_params.txt << EOF
$1    #UMUS
$2    #PHIS
$3    #SZA
$4    #PHI0
$5    #ZLEV
$6    #NLAY
$7    #SAMPLEGRID
$8    
$9    #Dirs
${10}
EOF



