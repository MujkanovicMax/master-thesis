#$linesza=$( grep "sza " uvspec.inp )
#$sza=${linesza:4:2}

function convert_flx2nc {
  if [ ! -e mc.flx.spc ]; then
    echo "Error, there is no data to convert!"
    exit
  fi
  cat > convert_flx.py << EOF
import xarray as xr
import numpy as np

mystic_output = './mc.flx.spc'
wvl, i, j, k, edir, edown, eup, adir, adown, aup  = np.loadtxt(mystic_output, unpack=True)

xdim, x_target_idx = np.unique(i,return_inverse=True)
ydim, y_target_idx = np.unique(j,return_inverse=True)
zdim, z_target_idx = np.unique(k,return_inverse=True)
wdim, w_target_idx = np.unique(wvl,return_inverse=True)

Edir = np.zeros((len(xdim), len(ydim), len(zdim), len(wdim)))
Edown = np.zeros((len(xdim), len(ydim), len(zdim), len(wdim)))
Eup = np.zeros((len(xdim), len(ydim), len(zdim), len(wdim)))

Edir[x_target_idx, y_target_idx, z_target_idx, w_target_idx] = edir
Edown[x_target_idx, y_target_idx, z_target_idx, w_target_idx] = edown
Eup[x_target_idx, y_target_idx, z_target_idx, w_target_idx] = eup



D = xr.Dataset({
    'x': xr.DataArray(xdim, dims=('x',)),
    'y': xr.DataArray(ydim, dims=('y',)),
    'z': xr.DataArray(zdim, dims=('z',)),
    'wvl': xr.DataArray(wdim, dims=('wvl',)),
    'Edir': xr.DataArray(Edir, dims=('x','y','z','wvl')),
    'Edown': xr.DataArray(Edown, dims=('x','y','z','wvl')),    
    'Eup': xr.DataArray(Eup, dims=('x','y','z','wvl')),
})
sza = 45  #45. #np.loadtxt("input_params.txt",dtype=float, skiprows=2, max_rows=1)
mu0 = np.cos(np.deg2rad(sza))
D.attrs["mu0"] = -mu0
D.to_netcdf('{}.nc'.format(mystic_output))
EOF
python3 convert_flx.py
rm convert_flx.py

}

convert_flx2nc

