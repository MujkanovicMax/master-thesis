#include <vgrid/MysticCloud.hpp>
#include <netcdf>
#include <iostream>
#include <array>
#include <camera.hpp>
#include <vgrid/Cyclic.hpp>
#include <cmath>
#include <numeric>
#include <Eigen/Dense>
#include "functions.hpp"


int main(int argc, char** argv) {
    std::cout << "Reading netcdf data" << "\n";
    //Reading in netcdf data
    std::vector<double> radiances;
    std::vector<double> mus;
    std::vector<double> phis;
    std::vector<double> wmus;
    std::vector<double> wphis;
    size_t nmu,nphi;

    {
        using namespace netCDF;

        NcFile file("/home/m/Mujkanovic.Max/ma/radiances/radiances_mu16_phi16.nc", NcFile::FileMode::read);
        size_t nx = file.getDim("x").getSize(); 
        size_t ny = file.getDim("y").getSize(); 
        size_t nz = file.getDim("z").getSize(); 
        nphi = file.getDim("phi").getSize(); 
        nmu = file.getDim("mu").getSize(); 
        size_t nwvl  = file.getDim("wvl").getSize();
        std::cout << nx*ny*nz*nmu*nphi*nwvl << "\n";
        radiances.resize(nx*ny*nz*nmu*nphi*nwvl);
        mus.resize(nmu);
        phis.resize(nphi);
        wmus.resize(nmu);
        wphis.resize(nphi);

        file.getVar("radiance").getVar(radiances.data());
        file.getVar("mu").getVar(mus.data());
        file.getVar("phi").getVar(phis.data());
        file.getVar("wmu").getVar(wmus.data());
        file.getVar("wphi").getVar(wphis.data());
    }
    
    std::vector<double> Edir;
    std::vector<double> Edown;
    Eigen::Vector3d sza_dir;
    double muEdir;

    {
        using namespace netCDF;

        NcFile file("/home/m/Mujkanovic.Max/ma/radiances/job_flx/mc.flx.spc.nc", NcFile::FileMode::read);
        size_t Nx = file.getDim("x").getSize(); 
        size_t Ny = file.getDim("y").getSize(); 
        size_t Nz = file.getDim("z").getSize(); 
        size_t Nwvl = file.getDim("wvl").getSize();

        file.getAtt("mu0").getValues(&muEdir);
        sza_dir = angleToVec(muEdir,270);
//        std::cout << "sza_dir = " << sza_dir.transpose() << "\n";
        Edir.resize(Nx*Ny*Nz*Nwvl);
        Edown.resize(Nx*Ny*Nz*Nwvl);
        file.getVar("Edir").getVar(Edir.data());
        file.getVar("Edown").getVar(Edown.data());


    }




    std::vector<double> kext;
    std::vector<double> zlev;
    std::vector<double> w0,g1;
    size_t nlev, nlyr, nx, ny;

    { 
        using namespace netCDF;

        NcFile file("test.optical_properties.nc", NcFile::FileMode::read);
        nlev = file.getDim("nlev").getSize();
        nlyr = file.getDim("caoth3d_0_wc_nlyr").getSize();
        nx = file.getDim("caoth3d_0_wc_Nx").getSize();
        ny = file.getDim("caoth3d_0_wc_Ny").getSize();

        zlev.resize(nlev);
        kext.resize(nlyr*nx*ny);
        w0.resize(nlyr*nx*ny);
        g1.resize(nlyr*nx*ny);
        file.getVar("caoth3d_0_wc_ext").getVar(kext.data());
        file.getVar("output_mc.z").getVar(zlev.data());
        file.getVar("caoth3d_0_wc_ssa").getVar(w0.data());
        file.getVar("caoth3d_0_wc_g1").getVar(g1.data());

    }
    for(auto& ke: kext) {
        ke *= 1000;
    }
    std::cout << "done" << "\n";


    // Grid Parameters and grid declaration
    double dx = 0.1;
    double dy = 1;
    double albedo = 0.2;
    auto mgrid = rayli::vgrid::MysticCloud(dx, dy, nx, ny, zlev);
    auto grid = rayli::vgrid::Cyclic(mgrid, {{0, 0, -std::numeric_limits<double>::infinity()},{nx*dx, ny*dy, std::numeric_limits<double>::infinity()}});
    std::cout <<  dx << " " << dy << " " << " " << nx << " " << ny << "\n";

    // Camera parameters and camera declaration
    size_t Nxpixel = 90;
    size_t Nypixel = 1;
    double fov = 2;
    auto loc = Eigen::Vector3d{4,0.01,2};
    double fovx = fov;
    double fovy = fov * Nxpixel / Nypixel;
    size_t rays = 90;
    auto cam = MysticPanoramaCamera(loc, 0, 0, -45, 45, 90, 90, 90);

    // Stream and substream calculations
    size_t nsub = 10;
    std::vector<double> streams = calcStreamDirs(mus,phis,nmu,nphi);
    std::vector<double> substreams = calcSubstreamDirs(mus,phis,wmus,wphis,nmu,nphi,nsub);
    
    //Image Declarations
    std::vector<double> image(Nxpixel * Nypixel);
    std::vector<double> opthick_image(Nxpixel * Nypixel);
    std::vector<double> Ldiff_i(Nxpixel * Nypixel);
    std::vector<double> Lup_i(Nxpixel * Nypixel);
    std::vector<double> Ldown_i(Nxpixel * Nypixel);
    std::vector<double> Ldir_i(Nxpixel * Nypixel);
    std::vector<double> groundbox(Nxpixel * Nypixel);
    std::vector<double> gRdir_i(Nxpixel * Nypixel);
    std::vector<double> gRdiff_i(Nxpixel * Nypixel);
    
    //Loop variables

    //Main loop
    std::cout << "Starting ray tracing..." << "\n";
    for(size_t i = 0; i < Nypixel; ++i) {
        double ypx = (i + 0.5) / Nypixel;
        //for(size_t j = 0; j < 1; ++j) {
        for(size_t j = 0; j < Nxpixel; ++j) {
            double xpx = (j + 0.5) / Nxpixel;
            auto ray = cam.compute_ray(Eigen::Vector2d{xpx, ypx});
            //auto ray = Ray{loc, Eigen::Vector3d{0.5, 0, -1}.normalized()};
            double optical_thickness = 0;
            double radiance = 0;
            double Lups = 0;
            double Ldowns = 0;
            double Ldirs = 0;
            double Ldiffs = 0;
            size_t groundidx;
            for(auto slice: grid.walk_along(ray, 0., std::numeric_limits<double>::infinity())) {
                if(auto pvol = std::get_if<VolumeSlice>(&slice)) {
                    auto [x,y,z] = indexDecompose<3>(pvol->idx, {nx,ny,nlyr});
                    size_t optprop_index = indexRecompose(std::array{z,x,y},std::array{nlyr,nx,ny});
                    size_t rad_index = indexRecompose(std::array{x,y,z+1},std::array{nx,ny,nlyr+1});
                    groundidx = indexRecompose(std::array{x,y,z},std::array{nx,ny,nlyr+1});
                    double L = Edir[rad_index]/fabs(muEdir);
                    double transmission = exp(-optical_thickness);
                    double dtau = (pvol->tfar - pvol->tnear) * kext[optprop_index];
                    optical_thickness += dtau;
                    double phase_function = phase_HG(g1[optprop_index], (-ray.d).dot(sza_dir.normalized()));
                    auto [tmp, Lup, Ldown] = calc_Ldiff(ray,pvol->tfar, pvol->tnear, pvol->idx, g1[optprop_index], nx, ny, nlyr, nmu, nphi, mus, phis, wmus, wphis, radiances, streams, nsub, substreams);
                    double scatter_prob = 1 - exp(-dtau*w0[optprop_index]);
                    std::cout << " Lup = " << Lup << " Lup_sum = " << Lups <<" transm = "<< transmission << " w0 = " << w0[optprop_index] << " dtau = " << dtau  << " scatter_prob = " << scatter_prob << "\n";
                    Lup *= transmission * scatter_prob;
                    Ldown *= transmission * scatter_prob;
                    double Ldir = transmission * L * scatter_prob * phase_function;
                    Ldiffs += Lup + Ldown;
                    Lups += Lup;
                    Ldowns += Ldown;
                    Ldirs += Ldir;
                    radiance += Ldir + Lup + Ldown;

                }

            }

            auto [gR, gRdir, gRdiff] = groundReflection_lambert(ray,groundidx,albedo, nx, ny, nlyr, nmu, nphi, mus, wmus, wphis, Edir, radiances); 
            double ground_E = gR *  exp(-optical_thickness);            
            image[j + i * Nxpixel] = radiance + ground_E;
            opthick_image[j+i*Nxpixel] = optical_thickness;
            Ldiff_i[j+i*Nxpixel] = Ldiffs;
            Lup_i[j+i*Nxpixel] = Lups;
            Ldown_i[j+i*Nxpixel] = Ldowns;
            Ldir_i[j+i*Nxpixel] = Ldirs;
            groundbox[j+i*Nxpixel] = indexDecompose<3>(groundidx,std::array<size_t,3>{nx,ny,nlyr})[0];
            gRdir_i[j+i*Nxpixel] = gRdir * exp(-optical_thickness);
            gRdiff_i[j+i*Nxpixel] = gRdiff * exp(-optical_thickness);

           
            std::cout << "Pixel " << j+1 << " done" << "  ival = " << image[j+i*Nxpixel] << " " <<  Ldiff_i[j+i*Nxpixel] << " "  << Lup_i[j+i*Nxpixel] << " " << Ldown_i[j+i*Nxpixel]<< " " << Ldir_i[j+i*Nxpixel] <<"\n";
        }
    }

    {
        using namespace netCDF;    
        NcFile file("output2.nc", NcFile::FileMode::replace);
        auto xdim = file.addDim("x", Nxpixel);
        auto ydim = file.addDim("y", Nypixel);
        file.addVar("image", NcType::nc_DOUBLE, {ydim, xdim}).putVar(image.data());
        file.addVar("optical_thickness_image", NcType::nc_DOUBLE, {ydim, xdim}).putVar(opthick_image.data());
        file.addVar("Ldiff", NcType::nc_DOUBLE, {ydim, xdim}).putVar(Ldiff_i.data());
        file.addVar("Lup", NcType::nc_DOUBLE, {ydim, xdim}).putVar(Lup_i.data());
        file.addVar("Ldown", NcType::nc_DOUBLE, {ydim, xdim}).putVar(Ldown_i.data());
        file.addVar("Ldir", NcType::nc_DOUBLE, {ydim, xdim}).putVar(Ldir_i.data());
        file.addVar("groundbox", NcType::nc_UINT64, {ydim, xdim}).putVar(groundbox.data());
        file.addVar("gRdir", NcType::nc_DOUBLE, {ydim, xdim}).putVar(gRdir_i.data());
        file.addVar("gRdiff", NcType::nc_DOUBLE, {ydim, xdim}).putVar(gRdiff_i.data());

    }
    


}
