import xarray as xr
import numpy as np

def calc_irr(ds, var):

    rads = ds[var]
    mus = ds["mu"]
    phis = ds["phi"]
    wmus = ds["wmu"]
    wphis = ds["wphis"]
    x = ds["x"]
    y = ds["y"]
    z = ds["z"]
    
    E = np.zeros(rads.shape)
    Eup = np.zeros((ds["x"].shape[0],ds["y"].shape[0],ds["z"].shape[0]))
    Edown = np.zeros(Eup.shape)
    ind_up = np.where(mus >= 0)
    ind_down = np.where(mus < 0)

    for i,x_i in enumerate(x):
        for j,y_j in enumerate(y):
            for k,z_k in enumerate(z):
                E[i,j,k,:,:] = rads[i,j,k,:,:] * mus * wmus * np.reshape(wphis,(wphis.shape[0],1))
                Edown[i,j,k] = np.sum(E[i,j,k,0:mus.shape[0]//2,:])
                Eup[i,j,k] = np.sum(E[i,j,k,mus.shape[0]//2:,:])
                #for m,mu in enumerate(mus):
                    #for p,phi in enumerate(phis):

                     #   if mu > 0:
                      #      Eup[i,j,k] += mu * rads[i,j,k,m,p] * wmus[m] * wphis[p]
    return Eup, Edown


ds = xr.open_dataset

