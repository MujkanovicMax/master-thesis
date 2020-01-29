import numpy as np
import netCDF4

def int_netcdf(x,y,wx,wy,fname,var):

    sum_arr = netCDF4.Dataset("job_back/job_" + x[0] + "_" + y[0] + "/"+ fname, "a")
    sum_arr[var][:] = sum_arr[var][:] * 0
    I=0

    for ix,m in enumerate(x):
        for iy,p in enumerate(y):
            fpath = "job_back/job_" + str(m) + "_" + str(p) + "/"+ fname
            arr = netCDF4.Dataset(fpath, "a")
            sum_arr[var][:] += float(m) * arr[var][:] * wx[ix] * wy[iy]
            arr.close()
            I+=1
            print(str(I),end="\r",flush=True)

    return sum_arr


UMUS = np.loadtxt("input_params.txt",dtype=str, max_rows=1)
PHIS = np.loadtxt("input_params.txt",dtype=str,skiprows=1, max_rows=1)
wumu = np.loadtxt("numus.txt", skiprows=1, max_rows=1)
wphi = np.loadtxt("nphis.txt", skiprows=1, max_rows=1)


UMUS_up = UMUS[np.where(UMUS.astype(float) >= 0)[0]]
UMUS_down = UMUS[np.where(UMUS.astype(float) <= 0)[0]]
wumu_down = wumu[0:len(UMUS_down)]
wumu_up = wumu[len(UMUS_down)-1:]

Eup = int_netcdf(UMUS_up,PHIS,wumu_up,wphi,"mc.rad.spc.nc","radiance")

Edown = int_netcdf(UMUS_down,PHIS,wumu_down,wphi,"mc.rad.spc.nc","radiance")







