import numpy as np
import matplotlib.pyplot as plt
import xarray as xr


def HG(g, mu):
    f = 1+g*g-2*g*mu
    return 1./(4*np.pi) * (1-g*g)/np.sqrt(f*f*f)

def HG_a(g, theta):
    f = 1+g*g-2*g*np.cos(theta)
    return 1./(4*np.pi) * (1-g*g)/np.sqrt(f*f*f)

def HG2D(g,x, y):
    
    val = np.zeros((len(x),len(y)))
    phi = np.linspace(0,360,32, endpoint=False)
    for i in range(32):
            x = np.cos(phi/180*np.pi)
            y = np.sin(phi/180*np.pi)

            init = np.array([x,y,0])/np.sqrt(x*x + y*y)
            end = np.array([x,y,-5])/np.sqrt(x*x + y*y + 5*5)

            mu = np.dot(-end, init)

            val[i,j] =  HG(g, mu)
    
    return val


def HG_2D(g,x,y):

    val = np.zeros((len(x),len(y)))

    for i in range(len(x)):
        for j in range(len(y)):

            init = np.array([0,0,-1])#/np.sqrt(x[i]*x[i] + y[j]*y[j])
            end = np.array([x[i],y[j],-5])/np.sqrt(x[i]*x[i] + y[j]*y[j] + 5*5)

            mu = np.dot(-end, init)

            val[i,j] =  HG(g, mu)
    
    ds = xr.Dataset({
            'x': xr.DataArray(x, dims=('x',)),
            'y': xr.DataArray(y, dims=('y',)),
            'HG': xr.DataArray(val*10000, dims=('x','y')),
            })
    ds.to_netcdf("nsubproblem/HG_radial.nc")
    return val



g = 0.86
mu = np.arange(-1,1,0.01)
theta = np.arange(-180,180,0.1)/180*np.pi

x = np.arange(-5,5,0.01)
y = np.arange(-5,5,0.01)

bg = np.ones((len(x),len(y))) * 0.005
HG_2D(g,x,y)

#fig,ax = plt.subplots()
#
##ax.plot(theta, HG_a(g,theta))
#im=ax.contourf(x,y, HG_2D(g,x,y))
#plt.colorbar(im)
#ax.set_xlabel("x")
#ax.set_ylabel("y")
#plt.tight_layout()
#plt.show()
