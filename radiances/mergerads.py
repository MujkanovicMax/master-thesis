import xarray as xr
import numpy as np

D = xr.open_mfdataset("job_*_*/mc.rad.spc.nc")

D.to_netcdf("radiances.nc")

