function convert_librad2nc {
  if [ ! -e mc.rad.spc ]; then
    echo "Error, there is no data to convert!"
    exit
  fi
  cat > convert_data.py << EOF
import xarray as xr
import numpy as np

mystic_output = './mc.rad.spc'

with open("uvspec.inp") as D:
	lines = D.readlines()
	for line in lines:
		if line.startswith('umu '):
			mu = float(line.split()[1])
		if line.startswith('phi '):
			phi = float(line.split()[1])
	
wvl, j, i, k, l = np.loadtxt(mystic_output, unpack=True)
xdim, x_target_idx = np.unique(i,return_inverse=True)
ydim, y_target_idx = np.unique(j,return_inverse=True)
zdim, z_target_idx = np.unique(k,return_inverse=True)
wdim, w_target_idx = np.unique(wvl,return_inverse=True)

radiance = np.zeros((len(xdim), len(ydim), len(zdim), len(wdim)))
radiance[x_target_idx, y_target_idx, z_target_idx, w_target_idx] = l
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

D.to_netcdf('{}.nc'.format(mystic_output))
EOF
python3 convert_data.py
rm convert_data.py

}

convert_librad2nc

