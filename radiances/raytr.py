import numpy as np
import netCDF4 as nd
import matplotlib.pyplot as plt

def get_sample_points(origin, dir, ns=10):

    gps = np.zeros((origin.shape[0],dir.shape[0]))
    gps[0] = np.tan(dir/180.*np.pi) * origin[2] + origin[0]
    gps[1,:] = origin[1]*1
    gps[2,:] = 0

    ds = np.linspace(0,1,ns)
    v = np.zeros((ds.shape[0],origin.shape[0],dir.shape[0]))

    v = np.reshape(origin,(1,origin.shape[0],1)) + (gps - np.reshape(origin,(origin.shape[0],1))) * np.reshape(ds,(ds.shape[0],1,1))
  
    return v


def get_boxindex(p,grid):
    
    index = np.zeros_like(p)
    for i in range(len(index)):
        index[i] = p[i]//grid[i+3]
        print(index[i])
        if index[i] < 0:
            index[i]+=grid[i]
        if index[i] > grid[i]-1:
            index[i] -= grid[i]
            print(index[i])


    return index.astype(int)



def get_e(p,grid,fpath="job_-0.577350_283.923047/mc.flx.spc.nc"):

    
    edirs = nd.Dataset(fpath,'r')
    Edir = np.zeros(p.shape[1])
    Edown = np.zeros(p.shape[1])
    Eup = np.zeros(p.shape[1])
    
    for I in range(p.shape[1]):
        
        
        i,j,k = get_boxindex(p[:,I],grid)
        print(j,i,k)

        Edir[I] = edirs["Edir"][j,i,k,:]
        Edown[I] = edirs["Edown"][j,i,k,:]
        Eup[I] = edirs["Eup"][j,i,k,:]




        
    edirs.close()

    return Edir,Edown,Eup

def get_rad(p,grid):

    data = np.loadtxt("input_params.txt",dtype=str,max_rows=2)
    umus = data[0]
    phis = data[1]
    L = np.zeros((umus.shape[0]*phis.shape[0],p.shape[1]))

    for I in range(p.shape[1]):

        i,j,k = get_boxindex(p[:,I],grid)
        
        N = 0

        for l,mu in enumerate(umus):
            for m,phi in enumerate(phis):
                
                rads = nd.Dataset("job_"+mu+"_"+phi+"/mc.rad.spc.nc","r")
                L[N,I] = rads["radiance"][j,i,k,:]
                rads.close()
                N=N+1

    return L,umus.astype(float),phis.astype(float)


#input = np.loadtxt("input_params.txt")
#sza = input[2]
#phi0 = input[3]
#layers = input[4]
#samplegrid = input[5]

Nx = 30
Ny = 1
Nz = 3
dx = 0.1
dy = 1
dz = 1
sza=45
mu = np.cos(sza/180.*np.pi)
albedo = 0.2

grid = np.array([Nx,Ny,Nz,dx,dy,dz])


camerapos = np.array([1.5,0.01,2])
camerafov = 90.
camerapixels = 90

pixeledges = np.linspace(-camerafov/2.,camerafov/2,camerapixels+1)
pixelangles = (pixeledges[0:-1] + pixeledges[1:])/2.
pixelvalues = np.zeros(pixelangles.shape)

pixelground = get_sample_points(camerapos,pixelangles)[-1]
Edir,Edown,Eup = get_e(pixelground,grid)
L,umus,phis = get_rad(pixelground,grid)
test=np.ravel(umus.reshape((2,1))*phis)

Ldown = L[np.where(test<=0)[0],:]


pixelvalues = Edir * albedo / np.pi + np.sum(Ldown*albedo/np.pi*np.array([umus[0],umus[0]]).reshape((2,1)),axis=1)

fig,ax = plt.subplots()

ax.plot(np.arange(pixelvalues.shape[0]),pixelvalues)
plt.tight_layout()
plt.show(plt.show())


