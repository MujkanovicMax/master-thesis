import xarray as xr
import numpy as np

ds = xr.open_dataset("test.optical_properties.nc")

da = ds["caoth3d_0_wc_nlyr"]

s = 20
e = 51

print(da[s:s+2].values)
arr = np.mean(da[s:s+2])
i = s+2
while i+2 < e:
    print(i)
    print(da[i:i+2].values)
    arr = xr.concat([arr, np.mean(da[i:i+2])],dim="caoth3d_0_wc_nlyr")
    i = i+2
print(da[i:e+1].values)
arr = xr.concat([arr,np.mean(da[i:e+1])],dim="caoth3d_0_wc_nlyr") 
            
