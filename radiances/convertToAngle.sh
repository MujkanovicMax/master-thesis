#!/bin/bash

function convertToAngle {

	cat > convert.py << EOF
import numpy as np

bog = np.loadtxt("$1", max_rows=1)
angles=bog/np.pi * 180

np.savetxt("cangles.txt",(angles,),fmt="%f")
EOF
python convert.py
rm convert.py
}

convertToAngle "$1"
