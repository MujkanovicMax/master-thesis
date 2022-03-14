import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def getphi(vec):

    angl = np.arctan2(vec[1],vec[0])
    #angl = angl/np.pi*180
    #angl = (angl + 360) % 360
    
    #return angl/180*np.pi
    return (angl + 2*np.pi) % (2*np.pi)



def rot_atob(v1,v2):
    
    v1 = v1/np.linalg.norm(v1)
    v2 = v2/np.linalg.norm(v2)

    ab = np.dot(v1,v2)
    axb = np.cross(v1,v2)
    naxb = np.linalg.norm(axb)
    G = np.array([[ab,-naxb,0],[naxb,ab,0],[0,0,1]])
    u = v1*1
    v = v2-(ab*v1)
    v = v/np.linalg.norm(v)
    w = axb
    w /= np.linalg.norm(w)
    Finv = np.array([[u[0],v[0],w[0]],[u[1],v[1],w[1]],[u[2],v[2],w[2]]])
    print(Finv)
    F = np.linalg.inv(Finv)
    print(F) 
    R = np.matmul(Finv,np.matmul(G,F))
    return R




streams = np.loadtxt("streamdirs.txt")
substreams = np.loadtxt("substreamdirs.txt")


s = streams.reshape((1,3))
ss = substreams.reshape((100,3))

s = s[0,:]

ss = ss[0:100,:]

U=ss[:,0]
V=ss[:,1]
W=ss[:,2]
x=0
y=0
z=0
u=s[0]
v=s[1]
w=s[2]
lim = 1.5

firstss = np.array([U[0],V[0],W[0]])
ray = np.array([-7,44,10])
ray = ray/np.linalg.norm(ray)
#rayphi = getphi(ray)
#ssphi = getphi(firstss)
#
#rotangl = np.pi + rayphi-ssphi
#
#a = np.cos(rotangl)
#b = np.sin(rotangl)
#R = np.array([[a,-b,0],[b,a,0],[0,0,1]])
R = rot_atob(firstss, -ray)

vecs = np.array([U,V,W])
vecs = np.matmul(R,vecs)
Ud = vecs[0]
Vd = vecs[1]
Wd = vecs[2]

firstssd = np.array([Ud[0],Vd[0],Wd[0]])
angl = np.arccos(np.dot(firstssd[0:2]/np.linalg.norm(firstssd[0:2]), -ray[0:2]/np.linalg.norm(ray[0:2])))


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim([-1.5, 1.5])
ax.set_ylim([-1.5, 1.5])
ax.set_zlim([-1.5, 1.5])
ax.quiver(1,1,z,ray[0],ray[1],ray[2], color="r")
#ax.quiver(1+ray[0],1+ray[1],ray[2],Ud[1:],Vd[1:],Wd[1:])
ax.quiver(1+ray[0],1+ray[1],ray[2],U[0],V[0],W[0],color="g")
ax.quiver(1+ray[0],1+ray[1],ray[2],Ud[0],Vd[0],Wd[0],color="k")
plt.show()

