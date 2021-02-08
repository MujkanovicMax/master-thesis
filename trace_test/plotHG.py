import numpy as np
import matplotlib.pyplot as plt



def HG(g, mu):
    f = 1+g*g-2*g*mu
    return 1./(4*np.pi) * (1-g*g)/np.sqrt(f*f*f)


g = 0.86
mu = np.arange(-1,1,0.01)


fig,ax = plt.subplots()

ax.plot(mu, HG(g,mu))
ax.set_xlabel("mu")
ax.set_ylabel("HG")
plt.tight_layout()
plt.show()
