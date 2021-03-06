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
    
    if cz_start == 0:
        rads = xr.concat([ds[radkey][:,:,cz_start:cz_end:div,:,:,:],ds[radkey][:,:,cz_end:,:,:,:]],dim=hkey)
    else:
        rads = xr.concat([ds[radkey][:,:,0:cz_start,:,:,:],ds[radkey][:,:,cz_start:cz_end:div,:,:,:],ds[radkey][:,:,cz_end:,:,:,:]],dim=hkey)
    out = ds.drop_dims(hkey)
    out[radkey] = rads
    
    if cz_start == 0:
        Edir = xr.concat([flx["Edir"][:,:,cz_start:cz_end:div],flx["Edir"][:,:,cz_end:]],dim="z")
        Edown = xr.concat([flx["Edown"][:,:,cz_start:cz_end:div],flx["Edown"][:,:,cz_end:]],dim="z")
        Eup = xr.concat([flx["Eup"][:,:,cz_start:cz_end:div],flx["Eup"][:,:,cz_end:]],dim="z")
    else:
        Edir = xr.concat([flx["Edir"][:,:,0:cz_start],flx["Edir"][:,:,cz_start:cz_end:div],flx["Edir"][:,:,cz_end:]],dim="z")
        Edown = xr.concat([flx["Edown"][:,:,0:cz_start],flx["Edown"][:,:,cz_start:cz_end:div],flx["Edown"][:,:,cz_end:]],dim="z")
        Eup = xr.concat([flx["Eup"][:,:,0:cz_start],flx["Eup"][:,:,cz_start:cz_end:div],flx["Eup"][:,:,cz_end:]],dim="z")

    out_flx = flx.drop_dims("z")
    out_flx["Edown"] = Edown
    out_flx["Eup"] = Eup
    out_flx["Edir"] = Edir
   
    if cz_start == 0:
        kext = xr.concat([op["caoth3d_0_wc_ext"][0:50//div,:,:],op["caoth3d_0_wc_ext"][cz_end:,:,:]],dim="caoth3d_0_wc_nlyr")
        w0 = xr.concat([op["caoth3d_0_wc_ssa"][0:50//div,:,:],op["caoth3d_0_wc_ssa"][cz_end:,:,:]],dim="caoth3d_0_wc_nlyr")
        g1 = xr.concat([op["caoth3d_0_wc_g1"][0:50//div,:,:],op["caoth3d_0_wc_g1"][cz_end:,:,:]],dim="caoth3d_0_wc_nlyr")
        opz = xr.concat([op["output_mc.z"][cz_start:cz_end:div],op["output_mc.z"][cz_end:]],dim="nlev")
    else:
        kext = xr.concat([op["caoth3d_0_wc_ext"][0:cz_start,:,:],op["caoth3d_0_wc_ext"][cz_start:cz_end:div,:,:],op["caoth3d_0_wc_ext"][cz_end:,:,:]],dim="caoth3d_0_wc_nlyr")
        w0 = xr.concat([op["caoth3d_0_wc_ssa"][0:cz_start,:,:],op["caoth3d_0_wc_ssa"][cz_start:cz_end:div,:,:],op["caoth3d_0_wc_ssa"][cz_end:,:,:]],dim="caoth3d_0_wc_nlyr")
        g1 = xr.concat([op["caoth3d_0_wc_g1"][0:cz_start,:,:],op["caoth3d_0_wc_g1"][cz_start:cz_end:div,:,:],op["caoth3d_0_wc_g1"][cz_end:,:,:]],dim="caoth3d_0_wc_nlyr")
        opz = xr.concat([op["output_mc.z"][0:cz_start],op["output_mc.z"][cz_start:cz_end:div],op["output_mc.z"][cz_end:]],dim="nlev")


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
    mu = ds[mukey]
    wmu = ds["wmu"]
    rads_1 = rads[:,:,:,:,::2,:]
    rads_2 = rads[:,:,:,:,1::2,:]
    rads = rads[:,:,:,:,::2,:]
    tmp_rads = (rads_1.values + rads_2.values)/2
    #print(rads.shape)
    tmp = (mu[::2].values + mu[1::2].values)/2
    mu = tmp*1
    tmp = wmu[::2].values + wmu[1::2].values
    wmu = tmp*1

    out = ds.drop_dims(mukey)
    out[mukey] = mu
    out["wmu"] = wmu
    #print(rads.shape)
    #print(out[mukey].shape, out["wmu"].shape)
    out[radkey] = xr.DataArray(tmp_rads, dims=("x", "y", "z", "wvl", "mu", "phi"))
    if outpath != None:
        out.to_netcdf(outpath)
    #return out
    return out


def half_phi(ds, radkey, phikey, outpath=None):

    rads = ds[radkey]
    phi = ds[phikey]
    wphi = ds["wphi"]
    rads_1 = rads[:,:,:,:,:,::2]
    rads_2 = rads[:,:,:,:,:,1::2]
    rads = rads[:,:,:,:,:,::2]
    tmp_rads = (rads_1.values + rads_2.values)/2
    #print(rads.shape)
    tmp = (phi[::2].values + phi[1::2].values)/2
    phi = tmp*1
    tmp = wphi[::2].values + wphi[1::2].values
    wphi = tmp*1

    out = ds.drop_dims(phikey)
    out[phikey] = phi
    out["wphi"] = wphi
    #print(rads.shape)
    #print(out[mukey].shape, out["wmu"].shape)
    out[radkey] = xr.DataArray(tmp_rads, dims=("x", "y", "z", "wvl", "mu", "phi"))
    if outpath != None:
        out.to_netcdf(outpath)
    #return out
    return out



#for div in [2,4,8]:
#    rad = half_mu(half_phi(rad,"radiance", "phi"), "radiance", "mu")
#    rad.to_netcdf("rad_philipp_16x16_meaned_to_" + str(16//div) + "x" + str(16//div) + ".nc")
#
#rad = half_phi(rad,"radiance","phi")
#rad.to_netcdf("rad_philipp_16x16_meaned_to_meaned_to_2x1.nc")

rad = xr.open_dataset("rad_checkerboard_32x32_meaned_to_32x32.nc")
op = xr.open_dataset("op_checkerboard_32x32.nc")
flx = xr.open_dataset("flx_checkerboard_32x32.nc")

for i in [32,16,8,4,2]:
    op = xr.open_dataset("op_checkerboard_32x32.nc")
    flx = xr.open_dataset("flx_checkerboard_32x32.nc")
    rad = xr.open_dataset("rad_checkerboard_32x32_meaned_to_{}x{}.nc".format(str(i),str(i)))
    for div in [1,2,5,10,25,50]:
        r,o,f = half_levels(rad, "radiance", "z", op, flx, div)
        r.to_netcdf("checkerboard_perf/rad_checkerboard_32x32_meaned_to_{}x{}_zdiv{}.nc".format(str(i),str(i),str(div)))
        o.to_netcdf("checkerboard_perf/op_checkerboard_32x32_zdiv{}.nc".format(str(div)))
        f.to_netcdf("checkerboard_perf/flx_checkerboard_32x32_zdiv{}.nc".format(str(div)))

rad = xr.open_dataset("rad_checkerboard_32x32_meaned_to_2x1.nc")
for div in [1,2,5,10,25,50]:
    r,o,f = half_levels(rad, "radiance", "z", op, flx, div)
    r.to_netcdf("checkerboard_perf/rad_checkerboard_32x32_meaned_to_2x1_zdiv{}.nc".format(str(div)))

rad = xr.open_dataset("irr_checkerboard_32x32.nc")
for div in [1,2,5,10,25,50]:
    r,o,f = half_levels(rad, "radiance", "z", op, flx, div)
    r.to_netcdf("checkerboard_perf/irr_checkerboard_32x32_zdiv{}.nc".format(str(div)))



#names = ["irr_from_10s_twostr_only","irr_from_mystic_mcipa","irr_from_mystic","irr_from_10s_base"]
#   
#flxn = ["../radiances/flx_10s_twostr_only.nc","../radiances/job_flx/flx_mcipa.nc","../radiances/job_flx/flx_base_backup.nc","../radiances/flx_10s_base.nc"]
#
#for iflx,fname in enumerate(names):
#
#
#    rad =xr.open_dataset(fname + ".nc")
#    op = xr.open_dataset("test.optical_properties.nc")
#    flx = xr.open_dataset(flxn[iflx])
#    #rad = rad.load()
#
#    rkey = "radiance"
#
#    for div in [1,2,5,10,25,50]:
#
#        out_rad, out_op, out_flx = half_levels(rad,rkey,"z",op,flx,div)
#        out_rad.to_netcdf(fname + "_div" +str(int(div)) + ".nc")
#        #out_op.to_netcdf("op_levels_div" + str(int(div)) + ".nc")
#        out_flx.to_netcdf("flx_levels_"+ fname  +"_div" + str(int(div)) + ".nc")
#    
#    rad.close()
#    op.close()
#    flx.close()
#

#nlist = [50,25,10,5,2,1]
#for n in nlist:
#
#    rad_org = xr.open_dataset("rad_levels_div"+ str(n)  +".nc")
#
#    rad_mu_phi = rad_org
#    rad_mu_phi.to_netcdf("rad_zdiv" + str(n) + "_mu_32_phi_32.nc")
#    #rad_phi = rad_org
#    #rad_mu.to_netcdf("rad_mu_32.nc")
#    #rad_phi.to_netcdf("rad_phi_32.nc")
#    for i in range(5):
#        rad_mu_phi = half_phi(half_mu(rad_mu_phi,"radiance", "mu"),"radiance", "phi")
#
#        rad_mu_phi.to_netcdf("rad_zdiv"+ str(n)  +"_mu_" + str(32//int(pow(2,i+1))) +"_phi_" + str(32//int(pow(2,i+1))) + ".nc")

#rad = xr.open_dataset("rad_mu_2_phi_2.nc")
#rad_n = half_phi(rad,"radiance","phi")
#rad_n.to_netcdf("rad_mu_2_phi_1.nc")



