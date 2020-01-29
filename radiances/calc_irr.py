import numpy as np
import netCDF4

def int_netcdf(x,y,wx,wy,fname,var):

   # start = netCDF4.Dataset("job_" + str(x[0]) + "_" + str(y[0]) + "/"+ fname, "r")
   # sum_arr = start[var]*wx[0]*wy[0]
   # start.close() 

    sum_arr = netCDF4.Dataset("job_" + x[0] + "_" + y[0] + "/"+ fname, "a")
    sum_arr[var][:] = sum_arr[var][:] * 0
    I=0

    for ix,m in enumerate(x):
        for iy,p in enumerate(y):
            fpath = "job_" + str(m) + "_" + str(p) + "/"+ fname
            arr = netCDF4.Dataset(fpath, "a")
            sum_arr[var][:] += arr[var] * wx[ix] * wy[iy]
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

print(UMUS_up,PHIS)

Eup = int_netcdf(UMUS_up,PHIS,wumu_up,wphi,"mc.rad.spc.nc","radiance")
Edown = int_netcdf(UMUS_down,PHIS,wumu_down,wphi,"mc.rad.spc.nc","radiance")



#start_down = netCDF4.Dataset("job_" + str(UMUS_down[0]) + "_" + str(PHIS[0]) + "/mc.rad.spc.nc", "r")
#start_up = netCDF4.Dataset("job_" + str(UMUS_up[0]) + "_" + str(PHIS[0]) + "/mc.rad.spc.nc", "r")
#
#sum_arr_down = start_down["radiance"] * weights_down[0]*wphi[0]
#sum_arr_up = start_up["radiance"] * weights_up[0]*wphi[0]
#
#
#
#for imu, mu in enumerate(UMUS_up[1:]):
#    for iphi, phi in enumerate(PHIS[1:]):
#        fpath = "job_" + str(mu) + "_" + str(phi) + "/mc.rad.spc.nc"
#        arr = netCDF4.Dataset(fpath, "r")
#        sum_arr_up += arr["radiance"] * weights_up[imu+1] * wphi[iphi+1]
#        arr.close()
#
#for imu, mu in enumerate(UMUS_down[1:]):
#    for iphi, phi in enumerate(PHIS[1:]):
#        fpath = "job_" + str(mu) + "_" + str(phi) + "/mc.rad.spc.nc"
#        arr = netCDF4.Dataset(fpath, "r")
#        sum_arr_down += arr["radiance"] * weights_down[imu+1] * wphi[iphi+1]
#        arr.close()
#



