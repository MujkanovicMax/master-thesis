import xarray as xr
import numpy as np








def cutoff_x(ds, radkey, xkey, op, flx, outpath=None):
    rads = ds[radkey]
    l = np.floor(np.log2(rads.shape[0]))
    l = int(pow(2,l))
    rads = rads[:l]
    x = ds[xkey]
    x = x[:l]
    out = ds.drop(radkey)
    out[xkey] = x
    out[radkey] = rads
   
    Edir = flx["Edir"][:l,:,:,:]
    Edown = flx["Edown"][:l,:,:,:]
    Eup = flx["Eup"][:l,:,:,:]
    out_flx = flx.drop_dims("x")
    out_flx["Edown"] = Edown
    out_flx["Eup"] = Eup
    out_flx["Edir"] = Edir

    
    kext = op["caoth3d_0_wc_ext"][:,:l,:]
    w0 = op["caoth3d_0_wc_ssa"][:,:l,:]
    g1 = op["caoth3d_0_wc_g1"][:,:l,:]
    out_op = op.drop_dims("caoth3d_0_wc_Nx")
    out_op["caoth3d_0_wc_ext"] = kext
    out_op["caoth3d_0_wc_ssa"] = w0
    out_op["caoth3d_0_wc_g1"] = g1

    
    if outpath != None:
        out.to_netcdf(outpath)
    return out, out_op, out_flx

def half_levels(ds, radkey, hkey, op, flx, div, outpath=None):

    lwc = op["caoth3d_0_wc_lwc"]
    ind = np.where(lwc > 2.e-30)
    cz_start = ind[0][0]
    cz_end = ind[0][-1] + 1


    rads = xr.concat([ds[radkey][:,:,cz_start:cz_end:div,:,:,:],ds[radkey][:,:,cz_end:,:,:,:]],dim=hkey)
    out = ds.drop_dims(hkey)
    out[radkey] = rads
    
    Edir = xr.concat([flx["Edir"][:,:,cz_start:cz_end:div],flx["Edir"][:,:,cz_end:]],dim="z")
    Edown = xr.concat([flx["Edown"][:,:,cz_start:cz_end:div],flx["Edown"][:,:,cz_end:]],dim="z")
    Eup = xr.concat([flx["Eup"][:,:,cz_start:cz_end:div],flx["Eup"][:,:,cz_end:]],dim="z")
    out_flx = flx.drop_dims("z")
    out_flx["Edown"] = Edown
    out_flx["Eup"] = Eup
    out_flx["Edir"] = Edir
    
    kext = xr.concat([op["caoth3d_0_wc_ext"][0:50//div,:,:],op["caoth3d_0_wc_ext"][cz_end:,:,:]],dim="caoth3d_0_wc_nlyr")
    w0 = xr.concat([op["caoth3d_0_wc_ssa"][0:50//div,:,:],op["caoth3d_0_wc_ssa"][cz_end:,:,:]],dim="caoth3d_0_wc_nlyr")
    g1 = xr.concat([op["caoth3d_0_wc_g1"][0:50//div,:,:],op["caoth3d_0_wc_g1"][cz_end:,:,:]],dim="caoth3d_0_wc_nlyr")
    opz = xr.concat([op["output_mc.z"][cz_start:cz_end:div],op["output_mc.z"][cz_end:]],dim="nlev")


    out_op = op.drop_dims("caoth3d_0_wc_nlyr").drop_dims("nlev")
    out_op["caoth3d_0_wc_ext"] = kext
    out_op["caoth3d_0_wc_ssa"] = w0
    out_op["caoth3d_0_wc_g1"] = g1
    out_op["output_mc.z"] = opz

    if outpath != None:
        out.to_netcdf(outpath)
    return out, out_op, out_flx


#def half_cloud_levels(rad, op, flx, div):
#
#    lwc = op["caoth3d_0_wc_lwc"]
#    ind = np.where(lwc > 2.e-30)
#    cz_start = ind[0][0]
#    cz_end = ind[0][-1] + 1
#    
#    rads = xr.concat([ds["radiance"][:,:,:cz_start,:,:,:],ds["radiance"][:,:,cz_start:cz_end+1:div,:,:,:],ds["radiance"][:,:,cz_end+1:,:,:,:]],dim="z")
#    out = ds.drop_dim("z")
#    out[radkey] = rads
#    
#    Edir = flx["Edir"][:,:,cz_start:cz_end+1:div,:]
#    Edown = flx["Edown"][:,:,cz_start:cz_end+1:div,:]
#    Eup = flx["Eup"][:,:,cz_start:cz_end+1:div,:]
#    out_flx = flx.drop_dims("z")
#    out_flx["Edown"] = Edown
#    out_flx["Eup"] = Eup
#    out_flx["Edir"] = Edir
#
#    
#
#    try:
#        kext = (op["caoth3d_0_wc_ext"][cz_start:cz_end:div,:,:] + op["caoth3d_0_wc_ext"][[cz_start+1:cz_end:div,:,:]) / 2
#        w0 = (op["caoth3d_0_wc_ssa"][[cz_start:cz_end:div,:,:] + op["caoth3d_0_wc_ssa"][[cz_start+1:cz_end:div,:,:]) / 2
#        g1 = (op["caoth3d_0_wc_g1"][[cz_start:cz_end:div,:,:] + op["caoth3d_0_wc_g1"][[cz_start+1:cz_end:div,:,:]) / 2
#        opz = op["output_mc.z"][[cz_start:cz_end:div]
#    except:
#        kext = (op["caoth3d_0_wc_ext"][::2,:,:][0:-1,:,:] + op["caoth3d_0_wc_ext"][1::2,:,:]) / 2
#        w0 = (op["caoth3d_0_wc_ssa"][::2,:,:][0:-1,:,:] + op["caoth3d_0_wc_ssa"][1::2,:,:]) / 2
#        g1 = (op["caoth3d_0_wc_g1"][::2,:,:][0:-1,:,:] + op["caoth3d_0_wc_g1"][1::2,:,:]) / 2
#        opz = op["output_mc.z"][::2][0:-1]
#
#
#    out_op = op.drop_dims("caoth3d_0_wc_nlyr").drop_dims("nlev")
#    out_op["caoth3d_0_wc_ext"] = kext
#    out_op["caoth3d_0_wc_ssa"] = w0
#    out_op["caoth3d_0_wc_g1"] = g1
#    out_op["output_mc.z"] = opz
#
#    if outpath != None:
#        out.to_netcdf(outpath)
#    return out, out_op, out_flx







def half_grid(ds, radkey, xkey, op, flx, outpath=None):

    rads = ds[radkey]
    rads_1 = rads[::2,:,:,:,:,:]
    rads_2 = rads[1::2,:,:,:,:,:]
    rads = rads[::2,:,:,:,:,:]
    rads.values = (rads_1.values + rads_2.values) / 2
    out = ds.drop_dims(xkey)
    out[radkey] = rads

    Edir = flx["Edir"][::2,:,:,:]
    Edown = flx["Edown"][::2,:,:,:]
    Eup = flx["Eup"][::2,:,:,:]
    Edir.values = flx["Edir"][::2,:,:,:].values*0.5 + flx["Edir"].values[1::2,:,:,:]
    Edown.values = flx["Edown"][::2,:,:,:].values*0.5 + flx["Edown"].values[1::2,:,:,:]
    Eup.values = flx["Eup"][::2,:,:,:].values*0.5 + flx["Eup"].values[1::2,:,:,:]
    out_flx = flx.drop_dims("x")
    out_flx["Edown"] = Edown
    out_flx["Eup"] = Eup
    out_flx["Edir"] = Edir


    kext = (op["caoth3d_0_wc_ext"][:,::2,:] + op["caoth3d_0_wc_ext"][:,1::2,:]) / 2
    w0 = (op["caoth3d_0_wc_ssa"][:,::2,:] + op["caoth3d_0_wc_ssa"][:,1::2,:]) / 2
    g1 = (op["caoth3d_0_wc_g1"][:,::2,:] + op["caoth3d_0_wc_g1"][:,1::2,:]) / 2
    
    out_op = op.drop_dims("caoth3d_0_wc_Nx")
    out_op["caoth3d_0_wc_ext"] = kext
    out_op["caoth3d_0_wc_ssa"] = w0
    out_op["caoth3d_0_wc_g1"] = g1

    if outpath != None:
        out.to_netcdf(outpath)
    return out, out_op, out_flx


def half_mu(ds, radkey, mukey, outpath=None):

    rads = ds[radkey]
    rads_1 = rads[:,:,:,:,::2,:]
    rads_2 = rads[:,:,:,:,1::2,:]
    rads = rads[:,:,:,:,::2,:]
    rads.values = (rads_1.values + rads_2.values)/2
    rads[mukey].values = (rads_1[mukey].values + rads_2[mukey].values)/2
    rads["wmu"].values = (rads_1["wmu"].values + rads_2["wmu"].values)

    out = ds.drop(radkey)
    out[mukey] = rads[mukey]
    out["wmu"] = rads["wmu"]
    out[radkey] = rads
    if outpath != None:
        out.to_netcdf(outpath)
    return out

def half_phi(ds, radkey, mukey, outpath=None):
    
    rads = ds[radkey]
    rads_1 = rads[:,:,:,:,:,::2]
    rads_2 = rads[:,:,:,:,:,1::2]
    rads = rads[:,:,:,:,:,::2]
    rads.values = (rads_1.values + rads_2.values)/2
    rads[mukey].values = (rads_1[mukey].values + rads_2[mukey].values)/2
    rads["wphi"].values = (rads_1["wphi"].values + rads_2["wphi"].values)

    out = ds.drop(radkey)
    out[mukey] = rads[mukey]
    out["wphi"] = rads["wphi"]
    out[radkey] = rads
    
    
    if outpath != None:
        out.to_netcdf(outpath)
    return out


    




rad = xr.open_dataset("../radiances/radiances_mu32_phi32.nc")
op = xr.open_dataset("test.optical_properties.nc")
flx = xr.open_dataset("../radiances/job_flx/mc.flx.spc.nc")
rad = rad.load()

rkey = "radiance"

for div in [1,2,5,10,25,50]:

    out_rad, out_op, out_flx = half_levels(rad,rkey,"z",op,flx,div)
    out_rad.to_netcdf("rad_levels_div" + str(int(div)) + ".nc")
    out_op.to_netcdf("op_levels_div" + str(int(div)) + ".nc")
    out_flx.to_netcdf("flx_levels_div" + str(int(div)) + ".nc")



