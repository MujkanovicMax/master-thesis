import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

def calc_irr(rad):
    
    Edown = (rad.radiance * rad.wmu * rad.wphi  * np.abs(rad.mu)).where(rad.mu<=0).sum(dim=["mu","phi"])
    Eup = (rad.radiance * rad.wmu * rad.wphi  * np.abs(rad.mu)).where(rad.mu>=0).sum(dim=["mu","phi"])

    return Eup, Edown



def irr_to_netcdf(inputname, outputname):

    rad = xr.open_dataset(inputname)
    Eup, Edown = calc_irr(rad)
    E = rad.drop_dims("mu").drop_dims("phi")
    E["Eup"] = Eup
    E["Edown"] = Edown
    E.to_netcdf(outputname)
    return E

fname = "rad_mu_32.nc"

irr_to_netcdf(fname, "irr_from32x32.nc")


#rad = xr.open_dataset(fname)
#rad4 = xr.open_dataset("radiances_h.nc")
#irr = xr.open_dataset("./job_flx/mc.flx.spc.nc")
#irrm = xr.open_dataset("Edir.nc")

#Eup, Edown = calc_irr(rad)
#E = rad.drop_dims("mu").drop_dims("phi")
#E["Eup"] = Eup
#E["Edown"] = Edown
#E.to_netcdf("irr_from32x32.nc")

#Edown4 = (rad4.radiance * rad4.wmu * rad4.wphi  * np.abs(rad4.mu)).where(rad4.mu<0).sum(dim=["mu","phi"]).isel(z=0)


#fig,ax = plt.subplots()
#ax.plot(rad["x"],Edown[:,0,0],".-b",label="Integrated downward irradiance ("+ str(rad["phi"].shape[0]) +  "azimuth)")
#ax.plot(irr["x"],irr["Edown"][:,0,0],".-r" , label="Downward irradiance from file")
##ax.plot(irrm["x"],irrm["Edown"][:,0,0], ".-g", label="Downward irrandiance mean from files")
##ax.plot(rad4["x"],Edown4[:,0,0],".-k",label="Integrated downward irradiance (" + str(rad4["phi"].shape[0]) + " azimuth)")
#ax.set_xlabel("Boxindex x")
#ax.set_ylabel("Irradiance in W/m²")
#ax.set_xlim(0,69)
#ax.set_ylim(0)
#ax.grid(True, linestyle=":")
#ax.legend()
#plt.tight_layout()
#plt.show()
