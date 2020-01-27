function convert_librad2nc {
  if [ ! -e mc.rad.spc ]; then
    echo "Error, there is no data to convert!"
    exit
  fi
  cat > convert_data.py << EOF
import xarray as xr
import numpy as np

mystic_output = './mc.rad.spc'
wvl, j, i, _, l = np.loadtxt(mystic_output, unpack=True)

xdim, x_target_idx = np.unique(i,return_inverse=True)
ydim, y_target_idx = np.unique(j,return_inverse=True)
wdim, w_target_idx = np.unique(wvl,return_inverse=True)

radiance = np.zeros((len(xdim), len(ydim), len(wdim)))
radiance[x_target_idx, y_target_idx, w_target_idx] = l

D = xr.Dataset({
    'x': xr.DataArray(xdim, dims=('x',)),
    'y': xr.DataArray(ydim, dims=('y',)),
    'wvl': xr.DataArray(wdim, dims=('wvl',)),
    'radiance': xr.DataArray(radiance, dims=('x','y','wvl')),
    })

D.to_netcdf('{}.nc'.format(mystic_output))
EOF
python convert_data.py
rm convert_data.py
mv mc.rad.spc.nc $1
}
