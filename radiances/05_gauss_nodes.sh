function gauss_nodes {

	cat > gauss_nodes.py << EOF
import numpy as np

gq = np.polynomial.legendre.leggauss($1)

nodes = gq[0]
weights = gq[1]

#nodes = (nodes+1)/2.
#nodes = (nodes-1)/2.

np.savetxt("gauss_nodes.txt", (nodes,weights), fmt="%f")
EOF
python gauss_nodes.py
rm gauss_nodes.py
}

gauss_nodes $1
