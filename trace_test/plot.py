import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

a = xr.open_dataset("output.nc")
b = xr.open_dataset("../radiances/job_panorama/mc.rad.spc.nc")

fig, ax = plt.subplots()

ax.plot(np.arange(90),a["image"][0,:],".-b",label="local estimates")
ax.plot(np.arange(90),b["radiance"][0,:,0,0,0,0], ".-r", label="mystic panorama")
ax.legend()
plt.show()
