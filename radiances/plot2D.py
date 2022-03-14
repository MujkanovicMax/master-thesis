import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

def get_rmse(val,obs):
    diff = val - obs
    return np.sqrt(np.mean(diff*diff))/np.mean(obs)

angles_list = np.array([2,4,8,16,32,64])
rmse_mat = np.zeros((angles_list.shape[0],angles_list.shape[0]))
rmse_g_mat = np.zeros_like(rmse_mat)

irr = xr.open_dataset("./job_flx/flx_base_backup.nc")


for n,i in enumerate(angles_list):
    for m,j in enumerate(angles_list):
        rad = xr.open_dataset("radiances_mu{}_phi{}.nc".format(i,j))
        Edown = (rad.radiance * rad.wmu * rad.wphi  * np.abs(rad.mu)).where(rad.mu<0).sum(dim=["mu","phi"])
        rmse_mat[n,m] = get_rmse(Edown, irr["Edown"])
        rmse_g_mat[n,m] = get_rmse(Edown[:,0,0,0], irr["Edown"][:,0,0,0])

x, y = np.meshgrid(angles_list, angles_list)
z=rmse_g_mat*100

plt.figure(1)
plt.clf()
plt.xlabel("Number of polar angle samples")
plt.ylabel("Number of azimuth angle samples")
plt.title("RMSE of integrated irradiance to mystic irradiance")
lvls=[np.min(z),1.8,2,5,7,10,50,100,np.max(z)]
cs = plt.contour(x,y,z,levels=lvls,colors="k")
plt.clabel(cs, inline=1, fontsize=5)
csf = plt.contourf(x,y,z,levels=lvls)
cb = plt.colorbar(csf)
cb.set_label("RMSE in %")
plt.ylim(np.min(y)-2,np.max(y)+2)
plt.xlim(np.min(x)-2,np.max(x)+2)

#fig, ax = plt.subplots()
#cf = ax.contourf(x,y,rmse_mat)
#ax.set_xlabel("Number of Mu")
#ax.set_ylabel("Number of Phi")
#ax.set_title("RMSE of integrated irradiance to mystic irradiance")
#cbar1 = fig.colorbar(cf)
#cbar1.set_label("RMSE")

plt.tight_layout()
plt.show()
