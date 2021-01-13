import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


t = xr.open_dataset("mystic_cam_10000photons.nc")["radiance"]

divs = [1,2,5,10,25,50]
method = [1,2,3,4,5]
rmse = np.zeros((len(divs),len(method)))
for i,div in enumerate(divs):
    for j,met in enumerate(method):
        obs = xr.open_dataset("output_levels_div" + str(div) + "method" + str(met) + ".nc")
        rmse[i,j] = np.sqrt(np.mean((obs["image"] - t)*(obs["image"] - t)))


fig,ax = plt.subplots()
ax.axis("tight")
ax.axis("off")

cols = ["Default","Mean of Up and Down", "Only Top", "Only Bottom", "Adaptive Weights"]
rows = ["50 Cloud Levels", "25 Cloud Levels", "10 Cloud Levels", "5 Cloud Levels", "2 Cloud Levels", "1 Cloud Level"]

ax.table(cellText=rmse,rowLabels=rows, colLabels=cols, loc="center")
plt.show()
