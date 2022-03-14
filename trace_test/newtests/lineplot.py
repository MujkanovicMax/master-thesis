import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

#x=:, y=59


obs = xr.open_dataset("output2d_homogeneous.nc")
t = xr.open_dataset("panorama_homogeneous.nc")

x = np.arange(1,obs.image.shape[0]+1)


fig,ax = plt.subplots()

ax.plot(x,obs.image[:,59],color="b",label="Full radiance")
ax.fill_between(x,obs.image[:,59],obs.Ldiff[:,59],color="b")
ax.plot(x,obs.Ldiff[:,59],color="r",label="Diffuse radiance")
ax.fill_between(x,obs.Ldiff[:,59],obs.Lup[:,59],color="r")
ax.plot(x,obs.Lup[:,59],color="y",label="Upwards diffuse radiance")
ax.fill_between(x,obs.Lup[:,59],obs.Ldir[:,59],color="y")
ax.plot(x,obs.Ldir[:,59],color="g",label="Direct radiance")
ax.fill_between(x,obs.Ldir[:,59],obs.Ldown[:,59],color="g")
ax.plot(x,obs.Ldown[:,59],color="c",label="Downward diffuse radiance")
ax.fill_between(x,obs.Ldown[:,59],np.zeros(x.shape),color="c")
ax.plot(x,t.radiance[:,59,0,0],color="k",label="Full radiance (MYSTIC)")
ax.set_xlabel("Pixel")
ax.set_ylabel(r'Radiance in $W m^{-2} sr^{-1}$')
ax.set_xlim(1,120)
plt.xticks([1,20,40,60,80,100,120])
plt.grid(True)
plt.legend(loc="lower center",bbox_to_anchor=(0.5,-0.63))
#fig.tight_layout()
fig.subplots_adjust(bottom=0.25)
plt.show()


