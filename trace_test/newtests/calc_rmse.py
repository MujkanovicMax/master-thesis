import numpy as np
import xarray as xr
import argparse as ap

def calc_rmse(fobs, ft):
    obs = xr.open_dataset(fobs)
    t = xr.open_dataset(ft)

    try:
        rmse = np.sqrt(((obs.image - t.radiance)*(obs.image - t.radiance)).mean())/t.radiance.mean()
        bias = (obs.image.mean()/t.radiance.mean() - 1) * 100
    except:
        try:
            rmse = np.sqrt(((obs.image - t.image)*(obs.image - t.image)).mean())/t.image.mean()
            bias = (obs.image.mean()/t.image.mean() - 1) * 100
        except:
            rmse = np.sqrt(((obs.radiance - t.radiance)*(obs.radiance - t.radiance)).mean())/t.radiance.mean()
            bias = (obs.radiance.mean()/t.radiance.mean() - 1) * 100


    return float(rmse.values), float(bias.values)

def _main():

    parser = ap.ArgumentParser(description="calc rmse and bias of two nc files")
    parser.add_argument('-o', type=str)
    parser.add_argument('-t', type=str)
    args = parser.parse_args()
    
    rmse , bias = calc_rmse(args.o, args.t)
    print("RMSE = " + str(rmse) + ", bias = " + str(bias) + "%")

if __name__ == "__main__":
    _main()
