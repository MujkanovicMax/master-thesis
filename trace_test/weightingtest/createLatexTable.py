import numpy as np
import xarray as xr
import pandas as pd
from calc_rmse import *

zdivs = [10,25,50]
modes = ["Exponent 6","Exponent 4","Exponent 8","Only top face radiances","Only bottom face radiances","Averaged top and bottom radiances"]
t = "panorama_checkerboard_mystic.nc"
rmse = np.zeros(len(modes))
bias = np.zeros(len(modes))

for zdiv in zdivs:
    for i,mode in enumerate(modes):
        obs = "output/output_rad_checkerboard_32x32_meaned_to_32x32_zdiv"+str(zdiv)+"_mode"+str(i)+".nc"
        rmse[i],bias[i] = calc_rmse(obs,t)
    
    print(zdiv)
    print("\n")
    df = pd.DataFrame(dict(MODE=modes,
                        RMSE=rmse,
                        BIAS=bias))
    print(df.to_latex(index=False))
                        
