import numpy as np
import xarray as xr
import matplotlib.pyplot as plt



t = xr.open_dataset("../radiances/rad_checkerboard_32x32/job_panorama/panorama_checkerboard.nc")


rmse = np.zeros((7,6))
bias = np.zeros((7,6))

for i,mu in enumerate([32,16,8,4,2]):
    for j,div in enumerate([1,2,5,10,25,50]):
        obs = xr.open_dataset("checkerboard_perf/output2d_rad_checkerboard_32x32_meaned_to_{}x{}_zdiv{}nc".format(str(mu),str(mu),str(div)))
        rmse[i,j] = np.sqrt(((obs.image - t.radiance)*(obs.image - t.radiance)).mean())/t.radiance.mean()
        bias[i,j] = (obs.image.mean()/t.radiance.mean() - 1) * 100

for j,div in enumerate([1,2,5,10,25,50]):
        obs = xr.open_dataset("checkerboard_perf/output2d_rad_checkerboard_32x32_meaned_to_2x1_zdiv{}nc".format(str(div)))
        rmse[5,j] = np.sqrt(((obs.image - t.radiance)*(obs.image - t.radiance)).mean())/t.radiance.mean()
        bias[5,j] = (obs.image.mean()/t.radiance.mean() - 1) * 100

for j,div in enumerate([1,2,5,10,25,50]):
        obs = xr.open_dataset("checkerboard_perf/output2d_irr_checkerboard_32x32_meaned_to_2x1_zdiv{}.nc".format(str(div)))
        rmse[6,j] = np.sqrt(((obs.image - t.radiance)*(obs.image - t.radiance)).mean())/t.radiance.mean()
        bias[6,j] = (obs.image.mean()/t.radiance.mean() - 1) * 100




fig,(ax1, ax2) = plt.subplots(2,1,figsize=(20,5))
ax1.axis("tight")
ax1.axis("off")
ax1.set_title("Normalized RMSE vs mystic panorama")
ax2.axis("tight")
ax2.axis("off")
ax2.set_title("Bias in %")
#cols = ["Default","Mean of Up and Down", "Only Top", "Only Bottom", "Adaptive Weights"]
rows = ["Nmu = 32, Nphi = 32","Nmu = 16, Nphi = 16","Nmu = 8, Nphi = 8", "Nmu = 4, Nphi = 4", "Nmu = 2, Nphi = 2", "Nmu = 2, Nphi = 1", "Irradiances"]
cols = ["50 Cloud Levels","25 Cloud Levels","10 Cloud Levels","5 Cloud Levels","2 Cloud Levels","1 Cloud Levels"]
table = ax1.table(cellText=np.around(rmse,4),rowLabels=rows, colLabels=cols, loc="center")
table2 = ax2.table(cellText=np.around(bias,4),rowLabels=rows, colLabels=cols, loc="center")
table.auto_set_font_size(False)
table.set_fontsize(8)
table.auto_set_column_width(col=list(range(len(cols))))
#table.scale(0.5, 0.5)

table2.auto_set_font_size(False)
table2.set_fontsize(8)
table2.auto_set_column_width(col=list(range(len(cols))))
#table2.scale(0.5, 0.5)
#plt.tight_layout()
#plt.show()
#plt.savefig("checkerboard_32x32_rmse_bias_v_mystic_table.pdf")





#fig,(ax1, ax2) = plt.subplots(2,1,figsize=(20,5))
#ax1.axis("tight")
#ax1.axis("off")
#ax1.set_title("Normalized RMSE vs 32x32 50 cloud layers")
#ax2.axis("tight")
#ax2.axis("off")
#ax2.set_title("Bias in %")
##cols = ["Default","Mean of Up and Down", "Only Top", "Only Bottom", "Adaptive Weights"]
#cols = ["Nmu = 2, Nphi = 1","Nmu = 2, Nphi = 2","Nmu = 4, Nphi = 4","Nmu = 8, Nphi = 8","Nmu = 16, Nphi = 16","Nmu = 32, Nphi = 32", "Ten-Stream Irradiances two stream only", "Mystic Irradiances mc ipa", "Mystic Irradiances normal", "Tenstream Irradiances normal"]
#rows = ["50 Cloud Levels", "25 Cloud Levels", "10 Cloud Levels", "5 Cloud Levels", "2 Cloud Levels", "1 Cloud Level"]
#table = ax1.table(cellText=np.around(rmse,4),rowLabels=rows, colLabels=cols, loc="center")
#table2 = ax2.table(cellText=np.around(bias,4),rowLabels=rows, colLabels=cols, loc="center")
#table.auto_set_font_size(False)
#table.set_fontsize(8)
#table.auto_set_column_width(col=list(range(len(cols))))
##table.scale(0.5, 0.5)
#
#table2.auto_set_font_size(False)
#table2.set_fontsize(8)
#table2.auto_set_column_width(col=list(range(len(cols))))
##table2.scale(0.5, 0.5)
##plt.tight_layout()
##plt.show()
#plt.savefig("rmse_bias_v_mystic_table.pdf")
