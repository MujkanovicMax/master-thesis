import numpy as np
import matplotlib.pyplot as plt


mie = np.loadtxt("phase_function_500nm_reff10.dat")
def hg(g,theta):
    #return 1/(4*np.pi) * (1-g*g)/np.power(1+g*g-2*g*np.cos(np.deg2rad(theta)),3./2.)

    return  (1-g*g)/np.power(1+g*g-2*g*np.cos(np.deg2rad(theta)),3./2.)


fig, ax = plt.subplots()

mieplot = ax.plot(mie[:,0], mie[:,1],color="grey",label="Mie phase function")
hgplot = ax.plot(mie[:,0], hg(0.8,mie[:,0]),color="k", label="Henyey-Greenstein phase function")
ax.set_xlabel("Theta in degrees")
ax.set_ylabel("Phase function")
plt.yscale("log")
plt.grid(True)
plt.legend()
plt.savefig("hg_mie_plot_log.pdf")


