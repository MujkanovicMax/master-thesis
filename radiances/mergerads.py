#!/usr/bin/env python

import xarray as xr
import numpy as np
import angles as a
import argparse
import os



def mergerads(nmus,m1,m2, nphis, fname="radiances.nc"):
    #xarray.set_options(file_cache_maxsize
    if os.path.exists(fname):
        print("File " + fname + " already exists. Radiances not merged")
        return 

    mus,wmus = a.gen_mus(nmus,m1,m2)
    phis,wphis = a.gen_phis(nphis)
    Dlist = [xr.open_mfdataset(["job_mu" + str(int(nmus))  + "/job_{}_{}/mc.rad.spc.nc".format(i,j) for j in phis],concat_dim="phi") for i in mus]
    D = xr.concat(Dlist, dim="mu")
    D["wphi"] = wphis
    D["wmu"] = wmus
    D["wmu"] = xr.DataArray(D.wmu.data, dims=("mu",))
    D["wphi"] = xr.DataArray(D.wphi.data, dims=("phi",))
    print("Merging radiances into netcdf. Info: ") 
    print(D)
    D.to_netcdf(fname)
    D.close()


def meanE(nmus,m1,m2, nphis, fname="Edir.nc"):
    mus,wmus = a.gen_mus(nmus,m1,m2)
    phis,wphis = a.gen_phis(nphis)
    E = xr.open_mfdataset(["job_mu" +str(int(nmus)) + "/job_{}_{}/mc.flx.spc.nc".format(i,j) for i in mus for j in phis],concat_dim="files")
    Em = E.mean(dim = "files")
    print("Mean irradiance from flx files. Info: ")
    print(Em)
    Em.to_netcdf(fname)
    E.close()



def _main():
    parser = argparse.ArgumentParser(description="merge radiances/ mean irradiances")
    parser.add_argument('-Nphi', type=int)
    parser.add_argument('-Nmu', nargs="+", type=float)
    parser.add_argument('-radfn', type=str)
    parser.add_argument('-meanfn', type=str)

    args = parser.parse_args()
    
    if args.Nphi is not None and args.Nmu is not None and args.radfn is not None:
        mergerads(args.Nmu[0], args.Nmu[1], args.Nmu[2], args.Nphi, args.radfn)
    if args.Nphi is not None and args.Nmu is not None and args.meanfn is not None:    
        meanE(args.Nmu[0], args.Nmu[1], args.Nmu[2], args.Nphi, args.meanfn)


if __name__ == "__main__":
    _main()




