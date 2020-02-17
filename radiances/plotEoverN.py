import numpy as np
import matplotlib.pyplot as plt
from calc_irr import * 

def gauss_nodes(n,min,max):
    
    gq = np.polynomial.legendre.leggauss(n)

    nodes = gq[0]
    weights = gq[1]

    
    nodes = (max - min)/2. * nodes + (min + max)/2.
    weights = weights * (max - min)/2.

    return nodes, weights

nangles = np.array([2,4,6,8,10,16,32])
Eups = np.zeros(nangles.shape[0])
Edowns = np.zeros(nangles.shape[0])
Eavgds = np.zeros(nangles.shape[0])
Eavgus = np.zeros(nangles.shape[0])

for i,na in enumerate(nangles):
    print("/n ("+str(na)+") /n")
    UMUS,wumu = gauss_nodes(na,-1,1)
    UMUS=UMUS.astype(str)
    for k in  range(len(UMUS)):
        UMUS[k] = np.format_float_positional(float(UMUS[k]), unique=False, precision=6)
    phis_rad,wphi = gauss_nodes(na,0,6.283185307)
    phis_rad = phis_rad.astype(str)
    for k in range(len(phis_rad)):
        phis_rad[k] = np.format_float_positional(float(phis_rad[k]), unique=False, precision=6)
    PHIS = phis_rad.astype(float) *180 /np.pi
    PHIS = PHIS.astype(str)
    for k in range(len(PHIS)):
        PHIS[k] = np.format_float_positional(float(PHIS[k]), unique=False, precision=6)
        
    
    eu,ed,ea = calc_Es(UMUS.astype(str),PHIS.astype(str),wumu,wphi,"mc.rad.spc.nc","radiance")
    Eups[i] = eu["radiance"][0,0,0,0]
    Edowns[i] = ed["radiance"][0,0,0,0]
    Eavgds[i] = ea["Edown"][0,0,0,0]
    Eavgus[i] = ea["Eup"][0,0,0,0]

fig,ax = plt.subplots(1,2)
ax[0].plot(nangles, Eups, ".-b", label="Eup")
ax[0].plot(nangles, Eavgus, ".-k", label="Eavg (from data)")
ax[0].legend()
ax[0].set_xlabel("Number of angles")
ax[0].set_ylabel("E")
ax[0].set_xticks(nangles)

ax[1].plot(nangles, Edowns, ".-b", label="Edown")
ax[1].plot(nangles, Eavgds, ".-k", label="Eavg (from data)")
ax[1].legend()
ax[1].set_xlabel("Number of angles")
ax[1].set_ylabel("E")
ax[1].set_xticks(nangles)
plt.tight_layout()
plt.show()
        


