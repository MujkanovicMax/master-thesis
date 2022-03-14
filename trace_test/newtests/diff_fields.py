import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


obs = xr.open_dataset("output_checkerboard_32x32_baseline.nc")
t = xr.open_dataset("panorama_checkerboard_mystic.nc")

plt1 = obs.image
plt2 = t.radiance[:,:,0,0]
plt3 = plt1 - plt2

pltlist = [plt1,plt2,plt3]
namelist = ["le_panorama_checkerboard","mystic_panorama_checkerboard","diff_le_myst_checkerboard"]
tmax = np.max(np.array([plt2.max(),plt1.max()]))
tmin = np.min(np.array([plt2.min(),plt1.min()]))
for i in range(3):
    fig,ax = plt.subplots()
    if i == 2:
        im = ax.imshow(pltlist[i])
    else:
        im = ax.imshow(pltlist[i],vmin = tmin, vmax=tmax)
    
    plt.xticks([0,20,40,60,80,100,119],[1,20,40,60,80,100,120])

    plt.yticks([0,20,40,60,80,100,119],[1,20,40,60,80,100,120])
    cb = plt.colorbar(im)
    cb.set_label(r'$W m^{-2} sr^{-1}$')
    plt.savefig(namelist[i] + ".pdf")

