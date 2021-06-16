import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

t = xr.open_dataset("mystic_cam_10000photons.nc")["radiance"]
#t = xr.open_dataset("output_zdiv1_mu_32_phi_32.nc")["image"]

divs = [1,2,5,10,25,50]
#method = [1,2,3,4,5]
na = [1,2,4,8,16,32]
#nsub = [1,2,4,10,20,35,50,100,200]
#rmse_sub = np.zeros(len(nsub))
rmse = np.zeros((len(divs),len(na)+4))
bias = np.zeros((len(divs),len(na)+4))

#for i in range(len(nsub)):
#    obs = xr.open_dataset("output_subsamp_plot_4x4_nsub_" + str(nsub[i]) + ".nc")
#    rmse_sub[i] = (np.sqrt(np.mean((obs["image"] - t)*(obs["image"] - t))))/np.mean(t)
#
#
#fig, ax = plt.subplots()
#
#ax.plot(nsub, rmse_sub)
#ax.set_xlabel("Number of subdivisions")
#ax.set_ylabel("Normalized RMSE")
#ax.grid(True)
#ax.set_ylim(0,0.7)
#ax.set_title("Normalized RMSE of different subdivisions of phase function against mystic panorama (2 mus, 2 phis)") 
#fig.tight_layout()
#plt.savefig("subsamp_RMSE.png")
for i,div in enumerate(divs):
    for j,n in enumerate(na):
        if n == 0:
            obs = xr.open_dataset("output_zdiv" + str(div) + "_mu_" + "2" + "_phi_" + "1" +  ".nc")
        else:
            obs = xr.open_dataset("output_zdiv" + str(div) + "_mu_" + str(n) + "_phi_" + str(n) +  ".nc")
        rmse[i,j] = (np.sqrt(np.mean((obs["image"] - t)*(obs["image"] - t))))/np.mean(t)#/(np.max(t) - np.min(t)) # oder np mean (t)
        bias[i,j] = ((np.mean(obs["image"])/np.mean(t))-1) *100
        print(div,n,"x",n,rmse[i,j])



for j,name in enumerate(["irr_from_10s_twostr_only","irr_from_mystic_mcipa","irr_from_mystic","irr_from_10s_base"]):
    for i,div in enumerate(divs):
        
        obs = xr.open_dataset("output_" + name  + "_div" + str(div) + ".nc")
        rmse[i,6+j] = (np.sqrt(np.mean((obs["image"] - t)*(obs["image"] - t))))/np.mean(t)#/(np.max(t) - np.min(t)) # oder np mean (t)
        bias[i,6+j] = (np.mean(obs["image"])/np.mean(t) - 1)*100

#for i,div in enumerate(divs):
#
#    obs = xr.open_dataset("output_irr_from32x32_myst_zdiv" + str(div) + ".nc")
#    rmse[i,7] = (np.sqrt(np.mean((obs["image"] - t)*(obs["image"] - t))))/np.mean(t)#/(np.max(t) - np.min(t)) # oder np mean (t)
#    bias[i,7] = (np.mean(obs["image"])/np.mean(t)-1)*100


fig,(ax1, ax2) = plt.subplots(2,1,figsize=(20,5))
ax1.axis("tight")
ax1.axis("off")
ax1.set_title("Normalized RMSE vs 32x32 50 cloud layers")
ax2.axis("tight")
ax2.axis("off")
ax2.set_title("Bias in %")
#cols = ["Default","Mean of Up and Down", "Only Top", "Only Bottom", "Adaptive Weights"]
cols = ["Nmu = 2, Nphi = 1","Nmu = 2, Nphi = 2","Nmu = 4, Nphi = 4","Nmu = 8, Nphi = 8","Nmu = 16, Nphi = 16","Nmu = 32, Nphi = 32", "Ten-Stream Irradiances two stream only", "Mystic Irradiances mc ipa", "Mystic Irradiances normal", "Tenstream Irradiances normal"]
rows = ["50 Cloud Levels", "25 Cloud Levels", "10 Cloud Levels", "5 Cloud Levels", "2 Cloud Levels", "1 Cloud Level"]
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
plt.savefig("rmse_bias_v_mystic_table.pdf")
