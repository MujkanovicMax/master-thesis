function gauss_nodes {

	cat > gauss_nodes.py << EOF
import numpy as np

gq = np.polynomial.legendre.leggauss($1)

nodes = gq[0]
weights = gq[1]

#print(type($4),$3)
nodes = ($4 - $3)/2. * nodes + ($3 + $4)/2.
weights = weights * ($4 - $3)/2.

np.savetxt("$2", (nodes,weights), fmt="%f")
EOF
python gauss_nodes.py
rm gauss_nodes.py
}

gauss_nodes "$1" "$2" "$3" "$4"
