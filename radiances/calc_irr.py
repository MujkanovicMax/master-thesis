import numpy as np
import netCDF4

UMUS = np.loadtxt("input_params.txt", max_rows=1)
PHIS = np.loadtxt("input_params.txt", skiprows=1, max_rows=1).astype(int)
nlyr = int(np.loadtxt("input_params.txt", skiprows=3, max_rows=1))
weights = np.loadtxt("gauss_nodes.txt", skiprows=1, max_rows=1)
dphi = (PHIS[1] - PHIS[0])/180. * np.pi


UMUS_up = UMUS[np.where(UMUS >= 0)[0]]
UMUS_down = UMUS[np.where(UMUS <= 0)[0]]
weights_down = weights[0:len(UMUS_down)]
weights_up = weights[len(UMUS_down)-1:]


start_down = netCDF4.Dataset("job_" + str(UMUS_down[0]) + "_" + str(PHIS[0]) + "/mc.rad.spc.nc", "r")
start_up = netCDF4.Dataset("job_" + str(UMUS_up[0]) + "_" + str(PHIS[0]) + "/mc.rad.spc.nc", "r")

sum_arr_down = start_down["radiance"] * weights_down[0]*dphi
sum_arr_up = start_up["radiance"] * weights_up[0]*dphi



for imu, mu in enumerate(UMUS_up[1:]):
    for iphi, phi in enumerate(PHIS[1:]):
        fpath = "job_" + str(mu) + "_" + str(phi) + "/mc.rad.spc.nc"
        arr = netCDF4.Dataset(fpath, "r")
        sum_arr_up += arr["radiance"] * weights_up[imu+1] * dphi
        arr.close()

for imu, mu in enumerate(UMUS_down[1:]):
    for iphi, phi in enumerate(PHIS[1:]):
        fpath = "job_" + str(mu) + "_" + str(phi) + "/mc.rad.spc.nc"
        arr = netCDF4.Dataset(fpath, "r")
        sum_arr_down += arr["radiance"] * weights_down[imu+1] * dphi
        arr.close()


