import xarray as xr
import numpy as np

phis = np.loadtxt("input_params.txt", dtype=str, skiprows=1, max_rows=1)
wphis = np.loadtxt("nphis.txt",skiprows=1, max_rows=1)
wmus = np.loadtxt("numus.txt",skiprows=1, max_rows=1)
Dphi = [xr.open_mfdataset("job_*_{}/mc.rad.spc.nc".format(j),concat_dim="mu") for j in phis]
D = xr.concat(Dphi, dim="phi")
D["wphi"] = wphis
D["wmu"] = wmus
D["wmu"] = xr.DataArray(D.wmu.data, dims=("mu",))
D["wphi"] = xr.DataArray(D.wphi.data, dims=("phi",))
D.to_netcdf("radiances.nc")
D.close()

sza = np.loadtxt("input_params.txt",dtype=float, skiprows=2, max_rows=1)
mu0 = np.cos(sza/180. * np.pi)
E = xr.open_mfdataset("job_*_*/mc.flx.spc.nc",concat_dim="files")
Em = E.mean(dim = "files") 
Em.attrs["mu0"] = -mu0
Em.to_netcdf("Edir.nc")
E.close()

