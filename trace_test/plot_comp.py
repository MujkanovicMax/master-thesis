import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

mystic_pano = xr.open_dataset("../radiances/job_panorama2/mc.rad.spc.nc")
x = np.arange(1,91,1)



fig, ax = plt.subplots()

#ax.plot(x, mystic_pano["radiance"][:,0,0,0,0,0], ".-r", label = "mystic panorama")



rmse_list = []
for nmu in [2,4,8,16,32,64]:
    for nphi in [2,4,8,16,32,64]:

        #nphi=nmu*1   
        rayli_im = xr.open_dataset("output_mu" + str(nmu) + "_phi" + str(nphi) +".nc")
        rmse = np.sqrt(np.mean((rayli_im["image"][0,0:46] - mystic_pano["radiance"][0:46,0,0,0,0,0])*(rayli_im["image"][0,:] - mystic_pano["radiance"][:,0,0,0,0,0])))
        #ax.plot(x, rayli_im["image"][0,:], label="rayli/mystic_panorama, nmu=" +str(nmu)+", nphi=" + str(nphi))
        rmse_list.append(rmse)


z = np.array(rmse_list).reshape((6,6))
x, y = np.meshgrid([2,4,8,16,32,64],[2,4,8,16,32,64])
cs=ax.contourf(x,y,z, levels=[0,2,3,5,6,9,15,30,60])
ax.set_xlabel("nmu")
ax.set_ylabel("nphi")
cbar=fig.colorbar(cs)
cbar.ax.set_ylabel("RMSE")
#ax.set_xlabel("Pixel")
#ax.set_ylabel("Radiance in W/sr/m²")
#ax.set_xlim(0,90)
#ax.legend()
#plt.grid(True,linestyle=":")


plt.tight_layout()
plt.show()
