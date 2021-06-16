import numpy as np
import xarray as xr 



def convertEdiff(fname, outputname):
    
    x = 70
    y = 3
    z = 99
    streams = 10

    ediff_in = np.loadtxt(fname, skiprows=2)
    ediff = ediff_in.reshape((y,x,z,streams))
    Eup = ediff[:,:,:,0]
    Edown = ediff[:,:,:,1]
    E_bl_out = ediff[:,:,:,2]
    E_bl_in = ediff[:,:,:,3]
    E_tl_out = ediff[:,:,:,4]
    E_tl_in = ediff[:,:,:,5]
    E_bba_out = ediff[:,:,:,6]
    E_bba_in = ediff[:,:,:,7]
    E_tba_out = ediff[:,:,:,8]
    E_tba_in = ediff[:,:,:,9]
    
    
    xdim = np.arange(x)
    ydim = np.arange(y)
    zdim = np.arange(z)
    



    E = xr.Dataset({
        "x": xr.DataArray(xdim, dims=("x",)),
        "y": xr.DataArray(ydim, dims=("y",)),
        "z": xr.DataArray(zdim, dims=("z",)),
        "Eup": xr.DataArray(Eup, dims=("y","x","z")),
        "Edown": xr.DataArray(Edown, dims=("y","x","z")),
        "E_bl_out": xr.DataArray(E_bl_out, dims=("y","x","z")),
        "E_bl_in": xr.DataArray(E_bl_in, dims=("y","x","z")),
        "E_tl_out": xr.DataArray(E_tl_out, dims=("y","x","z")),
        "E_tl_in": xr.DataArray(E_tl_in, dims=("y","x","z")),
        "E_bba_out": xr.DataArray(E_bba_out, dims=("y","x","z")),
        "E_bba_in": xr.DataArray(E_bba_in, dims=("y","x","z")),
        "E_tba_out": xr.DataArray(E_tba_out, dims=("y","x","z")),
        "E_tba_in": xr.DataArray(E_tba_in, dims=("y","x","z")),
        })


    E.to_netcdf(outputname)

    return E

def convertEdir(fname, outputname):
    E = 0
    return E





Ediff = convertEdiff("ediff.h5","Ediff_tenstream.nc")








