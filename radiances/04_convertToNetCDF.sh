#!/bin/bash

function convert_librad2nc {
    if [ ! -e mc.rad.spc ]; then
        echo "Error, there is no data to convert!"
        exit
    fi
    cat > convert_data.py << EOF
import xarray as xr
import numpy as np

mystic_output = './mc.rad.spc'

flag = 0

try:
    with open("uvspec.inp") as D:
        lines = D.readlines()
        for line in lines:
            if line.startswith('umu '):
                mu = np.round(float(line.split()[1]),3)
            if line.startswith('phi '):
                phi = np.round(float(line.split()[1]),3)

except:
    flag = 1
    with open("uvspec_panorama.inp") as D:
        lines = D.readlines()
        for line in lines:
            if line.startswith('umu '):
                mu = np.round(float(line.split()[1]),3)
            if line.startswith('phi '):
                phi = np.round(float(line.split()[1]),3)


wvl, i, j, k, l = np.loadtxt(mystic_output, unpack=True)
xdim, x_target_idx = np.unique(i,return_inverse=True)
ydim, y_target_idx = np.unique(j,return_inverse=True)
zdim, z_target_idx = np.unique(k,return_inverse=True)
wdim, w_target_idx = np.unique(wvl,return_inverse=True)

radiance = np.zeros((len(xdim), len(ydim), len(zdim), len(wdim)))
radiance[x_target_idx, y_target_idx, z_target_idx, w_target_idx] = l
    


if flag == 0:
    radiance = radiance[...,np.newaxis,np.newaxis]
    D = xr.Dataset({
        'x': xr.DataArray(xdim, dims=('x',)),
        'y': xr.DataArray(ydim, dims=('y',)),
        'z': xr.DataArray(zdim, dims=('z',)),
        'wvl': xr.DataArray(wdim, dims=('wvl',)),
        'mu': xr.DataArray([mu,], dims=('mu',)),
        'phi': xr.DataArray([phi,], dims=('phi',)), 
        'radiance': xr.DataArray(radiance, dims=('x','y','z','wvl','mu','phi')),
    })

elif flag == 1:
    D = xr.Dataset({
        'x': xr.DataArray(xdim, dims=('x',)),
        'y': xr.DataArray(ydim, dims=('y',)),
        'z': xr.DataArray(zdim, dims=('z',)),
        'wvl': xr.DataArray(wdim, dims=('wvl',)),
        'radiance': xr.DataArray(radiance, dims=('x','y','z','wvl')),
    })

else:
    print("HOW")


D.to_netcdf('{}.nc'.format(mystic_output))
EOF
python3 convert_data.py
rm convert_data.py

}

if [ ! -e mc.rad.spc.nc ]; then convert_librad2nc; fi

