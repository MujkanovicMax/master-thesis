import numpy as np
import xarray as xr

a = xr.open_dataset("output.nc")
b = xr.open_dataset("../radiances/job_panorama/mc.rad.spc.nc")

diff = a["image"][0,:].values - b["radiance"][0,:,0,0,0,0].values
RMSE = np.sqrt(np.mean(diff * diff))

with open("rmse/rmse.txt","a") as f:
    print(RMSE,file=f)


