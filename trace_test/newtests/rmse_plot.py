import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from calc_rmse import *

nsubs = np.arange(0,11)
#nsubs = np.append(nsubs, [100,200,300])
rmse = np.zeros(nsubs.shape)
bias = np.zeros(nsubs.shape)
t = "/project/meteo/work/Mujkanovic.Max/maproject/trace_test/newtests/panorama_checkerboard_mystic.nc"
for i,nsub in enumerate(nsubs):
    obs = "/project/meteo/work/Mujkanovic.Max/maproject/trace_test/newtests/nsub/output_checkerboard_32x32_nsub" + str(nsub) + ".nc"
    rmse[i],bias[i] = calc_rmse(obs, t) 

fig,ax = plt.subplots()
ax.plot(nsubs, rmse)
ax.set_xlabel("Number of subdivisions")
ax.set_ylabel("RMSE")
ax.set_title("RMSE of checkerboard clouds over number of subdivions (res 32x32)")
locs,labels = plt.yticks()
plt.yticks(np.arange(0,np.nanmax(rmse)+0.1,0.1))
#plt.xscale("symlog")
#plt.xlim(0,300)
ax.grid(True)
fig.tight_layout()
fig.savefig("RMSE_32x32_nsub.pdf")
