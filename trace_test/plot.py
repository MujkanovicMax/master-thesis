import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

a = xr.open_dataset("output.nc")
b = xr.open_dataset("../radiances/job_panorama/mc.rad.spc.nc")
c = xr.open_dataset("../radiances/job_flx/mc.flx.spc.nc")
#d = xr.open_dataset("rayli_panorama.nc")
alt = xr.open_dataset("output2.nc")
groundbox = a["groundbox"][0,:].load().data.astype(int)

fig, ax = plt.subplots()

ax.plot(np.arange(1,91,1),a["image"][0,:],".-b",label="Local estimates/Rayli")
#ax.plot(np.arange(1,91,1),a["Ldown"][0,:], label="ldown")
#ax.plot(np.arange(1,91,1),a["Ldown"][0,:]+a["Lup"][0,:], label="ldown+lup")
#ax.plot(np.arange(1,91,1),a["Ldown"][0,:]+a["Lup"][0,:]+a["Ldir"][0,:], label="ldown+lup+ldir")
ax.plot(np.arange(1,91,1),b["radiance"][:,0,0,0,0,0], ".-r", label="Mystic panorama")
ax.plot(np.arange(1,91,1),alt["image"][0,:],".-k",label="Local estimates/Rayli (no sub sampling)")
#ax.plot(np.arange(1,91,1),a["Ldown"][0,:]+a["Lup"][0,:]+a["Ldir"][0,:]+a["gRdir"][0,:], label="ldown+lup+ldir+gRdir")
#ax.plot(np.arange(1,91,1),a["Ldown"][0,:]+a["Lup"][0,:]+a["Ldir"][0,:]++a["gRdiff"][0,:], label="ldown+lup+ldir+gRdir+grDiff")
#ax.plot(np.arange(1,91,1),groundbox, label="groundbox index")
#ax.plot(np.arange(1,91,1),(c["Edir"][groundbox,0,0]+c["Edown"][groundbox,0,0])*0.2/np.pi , label="ground reflection")
#ax.plot(np.arange(1,91,1),(c["Edir"][groundbox,0,0])*0.2/np.pi , label="ground reflection Edir")
#ax.plot(np.arange(1,91,1),(c["Edown"][groundbox,0,0])*0.2/np.pi , label="ground reflection Edown")
#ax.plot(np.arange(1,91,1), d["transmissivity"][:,0]*1347.745/0.707, label="rayli panorama")
ax.set_xlabel("Pixel")
ax.set_ylabel("Radiance in W/sr/m²")
ax.set_xlim(0,90)
ax.legend()
plt.grid(True,linestyle=":")
plt.tight_layout()
plt.show()
