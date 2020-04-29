function linspace {

    cat > linspace.py << EOF
import numpy as np

angles = np.linspace($3,$4,$1+1)[:-1]
weights = np.repeat(2*np.pi/$1,$1)
np.savetxt("$2",(angles,weights),fmt="%f")

EOF

python3 linspace.py
rm linspace.py
}

linspace "$1" "$2" "$3" "$4"

