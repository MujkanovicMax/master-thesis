import numpy as np
import netCDF4 as nd
import matplotlib.pyplot as plt
from calc_irr import *

def get_sample_points(origin, dir , ns=10):

    gps = np.zeros((origin.shape[0],dir.shape[0]))
    gps[0] = np.tan(dir/180.*np.pi) * origin[2] + origin[0]
    gps[1,:] = origin[1]*1
    gps[2,:] = 0

    ds = np.linspace(0,1,ns)
    v = np.zeros((ds.shape[0],origin.shape[0],dir.shape[0]))

    v = np.reshape(origin,(1,origin.shape[0],1)) + (gps - np.reshape(origin,(origin.shape[0],1))) * np.reshape(ds,(ds.shape[0],1,1))
    
    return v


def get_boxindex(p,grid):
   
    dxi = np.array([grid[3],grid[4],grid[5]])
    bounds = np.array([grid[0],grid[1],grid[2]])
    index = np.zeros_like(p)
   
    try:
        for i in range(len(index)):
            index[i] = p[i]//grid[i+3]
            if index[i] < 0:
                index[i]+=grid[i]
            if index[i] > grid[i]-1:
                index[i] -= grid[i]

    except:

        for i in range(p.shape[0]):
            for xi in range(3):
                for k in range(p.shape[2]):
            
                    index[i,xi,k] = p[i,xi,k]//dxi[xi]
                    if index[i,xi,k] < 0:
                        index[i,xi,k] += bounds[xi]
                    if index[i,xi,k] > bounds[xi]-1:
                        index[i,xi,k] -= bounds[xi]


    return index.astype(int)



def get_e(p,grid,fpath="job_0.183435_36.600028/mc.flx.spc.nc"):

    
    edirs = nd.Dataset(fpath,'r')
    Edir = np.zeros(p.shape[1])
    Edown = np.zeros(p.shape[1])
    Eup = np.zeros(p.shape[1])
    
    for I in range(p.shape[1]):
        
        
        i,j,k = get_boxindex(p[:,I],grid)

        Edir[I] = edirs["Edir"][j,i,k,:]
        Edown[I] = edirs["Edown"][j,i,k,:]
        Eup[I] = edirs["Eup"][j,i,k,:]




        
    edirs.close()

    return Edir,Edown,Eup

def get_rad(p,grid):

    Eup = np.zeros(p.shape[1])
    Edown = np.zeros_like(Eup)
    Eu,Ed = calc_Es(UMUS,PHIS,wumu,wphi,"mc.rad.spc.nc","radiance")

    for I in range(p.shape[1]):

        i,j,k = get_boxindex(p[:,I],grid)
        Eup[I] = Eu["radiance"][j,i,k,:]
        Edown[I] = Ed["radiance"][j,i,k,:]

    return Eup,Edown




UMUS = np.loadtxt("input_params.txt",dtype=str, max_rows=1)
PHIS = np.loadtxt("input_params.txt",dtype=str,skiprows=1, max_rows=1)
wumu = np.loadtxt("numus.txt", skiprows=1, max_rows=1)
wphi = np.loadtxt("nphis.txt", skiprows=1, max_rows=1)



Nx,Ny,dx,dy = np.loadtxt("input_params.txt", skiprows = 6, max_rows=1) 
Zlev = np.loadtxt("input_params.txt", skiprows = 4, max_rows=1)
Nz = 2 #Zlev.shape[0]
dz = 1
sza = np.loadtxt("input_params.txt", skiprows = 2, max_rows=1)
mu = np.cos(sza/180.*np.pi)
albedo = 0.2

grid = np.array([Nx,Ny,Nz,dx,dy,dz])

cloudx = np.array([3,4])
cloudy = np.array([0])
cloudz = np.array([0,1])


camerapos = np.array([1.5,0.01,2])
camerafov = 90.
camerapixels = 90

pixeledges = np.linspace(-camerafov/2.,camerafov/2,camerapixels+1)
pixelangles = (pixeledges[0:-1] + pixeledges[1:])/2.
pixelvalues = np.zeros(pixelangles.shape)

pixelground = get_sample_points(camerapos,pixelangles)[-1]
Edir,Edown,Eup = get_e(pixelground,grid)
Eu,Ed = get_rad(pixelground,grid)

pixelvalues = Edir * albedo / np.pi + Ed*albedo/np.pi 
palt = Edir *albedo/np.pi + Edown*albedo/np.pi

truth = nd.Dataset("job_panorama/mc.rad.spc.nc" , "r")

fig,ax = plt.subplots(3,1)

ax[0].plot(np.arange(pixelvalues.shape[0]),pixelvalues,label="with Edown from radiances")
ax[1].plot(np.arange(palt.shape[0]),palt,label="with Edown from flx file")
ax[2].plot(np.arange(truth["radiance"][0,:,0,0].shape[0]),truth["radiance"][0,:,0,0], label="from mystic panorama")
ax[0].legend()
ax[1].legend()
ax[2].legend()

plt.tight_layout()
plt.show()




