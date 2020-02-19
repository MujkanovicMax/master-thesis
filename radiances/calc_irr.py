import numpy as np
from contextlib import closing
import netCDF4

def int_netcdf(x,y,wx,wy,fname,var):

    path = "job_" + x[0] + "_" + y[0] + "/"+ fname
    with closing(netCDF4.Dataset(path)) as D:
        sum_arr = np.zeros_like(D[var][:])

    #sum_arr = netCDF4.Dataset("job_" + x[0] + "_" + y[0] + "/"+ fname, "a")
    #sum_arr[var][:] = sum_arr[var][:] * 0
    I=0

    for ix,m in enumerate(x):
        for iy,p in enumerate(y):
            fpath = "job_" + str(m) + "_" + str(p) + "/"+ fname
            arr = netCDF4.Dataset(fpath, "a")
            sum_arr += abs(float(m)) * arr[var][:] * wx[ix] * wy[iy]        #abs erstmal nur damit Edown nicht negativ wird
            print(ix, iy, p, fpath, sum_arr, m, wx[ix], wy[iy], arr[var][:])
            arr.close()
            I+=1
            #print(str(I),end="\r",flush=True)

    return { var: sum_arr }


def calc_Es(UMUS,PHIS,wumu,wphi,filename,var):


    UMUS_up = UMUS[np.where(UMUS.astype(float) >= 0)[0]]
    UMUS_down = UMUS[np.where(UMUS.astype(float) <= 0)[0]]
    wumu_down = wumu[0:len(UMUS_down)]
    wumu_up = wumu[len(UMUS_down):]


    Eup = int_netcdf(UMUS_up,PHIS,wumu_up,wphi,filename,var)

    Edown = int_netcdf(UMUS_down,PHIS,wumu_down,wphi,filename,var)

    Eavg = calc_avg_E(UMUS,PHIS)

    return Eup, Edown, Eavg


def calc_avg_E(x,y):

    sum_arr = netCDF4.Dataset("job_" + x[0] + "_" + y[0] + "/"+ "mc.flx.spc.nc", "a")
    edn = np.zeros_like(sum_arr["Edown"][:])
    eup = np.zeros_like(sum_arr["Eup"][:])
    #sum_arr["Edown"][:] = sum_arr["Edown"][:] * 0
    #sum_arr["Eup"][:] = sum_arr["Eup"][:] * 0
    I=0

    for ix,m in enumerate(x):
        for iy,p in enumerate(y):
            fpath = "job_" + str(m) + "_" + str(p) + "/"+ "mc.flx.spc.nc"
            arr = netCDF4.Dataset(fpath, "a")
            edn +=  arr["Edown"][:]
            eup +=  arr["Eup"][:]
            arr.close()
            I+=1
            print(str(I),end="\r",flush=True)

    edn /= (len(x)*len(y))
    eup /= (len(x)*len(y))
    
    return {"Edown": edn, "Eup": eup}


UMUS = np.loadtxt("input_params.txt",dtype=str, max_rows=1)
PHIS = np.loadtxt("input_params.txt",dtype=str,skiprows=1, max_rows=1)
wumu = np.loadtxt("numus.txt", skiprows=1, max_rows=1)
wphi = np.loadtxt("nphis.txt", skiprows=1, max_rows=1)

Eup, Edown, Eavg = calc_Es(UMUS,PHIS,wumu,wphi,"mc.rad.spc.nc","radiance")

