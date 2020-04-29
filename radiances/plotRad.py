import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import mergerads as mrads

rad = xr.open_dataset("radiances.nc")
irr = xr.open_dataset("./job_flx/mc.flx.spc.nc")

Edown = (rad.radiance * rad.wmu * rad.wphi  * np.abs(rad.mu)).where(rad.mu<0).sum(dim=["mu","phi"]).isel(z=0)

fig, ax = plt.subplots()
ax.plot(rad["x"],Edown[:,0,0],".-b",label="Integrated downward irradiance ("+ str(rad["phi"].shape[0]) +  "azimuth)")
ax.plot(irr["x"],irr["Edown"][:,0,0],".-r" , label="Downward irradiance from file")
for i in [4,8,16]:
    mrads.mergerads(16,-1,1,i,"radiances_{}.nc".format(i))
    radh = xr.open_dataset("radiances_{}.nc".format(i))
    Edownh = (radh.radiance * radh.wmu * radh.wphi  * np.abs(radh.mu)).where(radh.mu<0).sum(dim=["mu","phi"]).isel(z=0)
    ax.plot(radh["x"],Edownh[:,0,0],label="Integrated downward irradiance (" + str(radh["phi"].shape[0]) + " azimuth)")
ax.set_xlabel("Boxindex x")
ax.set_ylabel("Irradiance in W/m²")
ax.set_xlim(0,69)
ax.set_ylim(0)
ax.grid(True, linestyle=":")
ax.legend()
plt.tight_layout()
plt.show()


