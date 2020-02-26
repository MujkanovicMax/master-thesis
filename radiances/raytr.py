import numpy as np
import netCDF4 as nd

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



def get_edir(p,grid,fpath="job_-0.577350_283.923047/mc.flx.spc.nc"):

    
    edirs = nd.Dataset(fpath,'r')
    Edir = np.zeros(p.shape[1])


    for I in range(p.shape[1]):
        
        
        i,j,k = get_boxindex(p[:,I],grid)
        print(j,i,k)

        Edir[I] = edirs["Edir"][j,i,k,:]




        
    edirs.close()

    return Edir


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
mu = np.cos(45/180.*np.pi)
albedo = 0.2

grid = np.array([Nx,Ny,Nz,dx,dy,dz])


camerapos = np.array([1.5,0.01,2])
camerafov = 90.
camerapixels = 90

pixeledges = np.linspace(-camerafov/2.,camerafov/2,camerapixels+1)
pixelangles = (pixeledges[0:-1] + pixeledges[1:])/2.
pixelvalues = np.zeros(pixelangles.shape)

pixelground = get_sample_points(camerapos,pixelangles)[-1]
Edir = get_edir(pixelground,grid)

###Edir/mu * albedo *





