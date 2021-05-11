#include <vgrid/MysticCloud.hpp>
#include <netcdf>
#include <iostream>
#include <array>
#include <camera.hpp>
#include <vgrid/Cyclic.hpp>
#include <cmath>
#include <numeric>
#include <Eigen/Dense>
#include <string>
#include "functions.hpp"


int main(int argc, char** argv) {
    //for(int xg = 2; xg <= 64; xg = 2*xg){
      //  std::cout << "nmu = " << MU << "\n";
        for(int zdiv = 0; zdiv <= 0; ++zdiv){
        for(int adiv = 0; adiv <= 0; ++adiv){
            //std::cout << "nphi = " << PHI << "\n\n";
    //int zdiv = 1;
    int zdivs[9] {1,2,4,10,20,35,50,100,200}; // for level benchmark
    int adivs[6] {32,16,8,4,2,1}; // for radiance samples benchmark
    //int divs[1] {2};
    std::cout << "Reading netcdf data" << "\n";
    //Reading in netcdf data
    //std::string radfpath = "/home/m/Mujkanovic.Max/ma/radiances/radiances_mu" + std::to_string(MU) + "_phi" + std::to_string(PHI) + ".nc";
    //std::string radfpath = "rad_zdiv" + std::to_string(zdivs[zdiv])  + "_mu_" + std::to_string(adivs[adiv])  + "_phi_" + std::to_string(adivs[adiv]) + ".nc";
    std::string radfpath = "rad_mu_2_phi_2.nc";  //"irr_from32x32_myst_zdiv" + std::to_string(zdivs[zdiv])  + ".nc";
    //std::string radfpath = "rad_zdiv1_mu_2_phi_2.nc";
    //std::string radfpath = "../radiances/rad_philipp_16x16.nc";

    std::cout << "radfpath = " << radfpath << "\n\n";
    std::vector<double> radiances;
    std::vector<double> mus;
    std::vector<double> phis;
    std::vector<double> wmus;
    std::vector<double> wphis;
    size_t nmu,nphi;
    
    read_radiances(radfpath, radiances, mus, phis, wmus, wphis,  nmu, nphi);
   

    std::string flxfpath = "/home/m/Mujkanovic.Max/ma/radiances/job_flx/mc.flx.spc.nc";
    //std::string flxfpath = "flx_levels_div" + std::to_string(zdivs[zdiv]) + ".nc";
    //std::string flxfpath = "../radiances/flx_philipp_16x16.nc";
    std::cout << "flxfpath = " << flxfpath << "\n\n";
    std::vector<double> Edir;
    std::vector<double> Edown;
    Eigen::Vector3d sza_dir;
    double muEdir;
    
    read_flx( flxfpath, Edir, Edown, sza_dir, muEdir );
   

    //std::string opfpath = "op_levels_div" + std::to_string(zdivs[zdiv]) + ".nc";
    std::string opfpath = "test.optical_properties.nc";
    //std::string opfpath = "../radiances/opprop_philipp_16x16.nc";

    std::cout << "opfpath = " << opfpath << "\n\n";
    std::vector<double> kext;
    std::vector<double> zlev;
    std::vector<double> w0,g1;
    size_t nlev, nlyr, nx, ny;
    
    read_opprop( opfpath, kext, zlev, w0, g1, nlev, nlyr, nx, ny );
    
    for(auto& ke: kext) {
        ke *= 1000;
    }
    std::cout << "done" << "\n";
    

    // Filename for output
    std::string outputfname = "refactest_n";



    // Grid Parameters and grid declaration
    double dx = 0.1;
    double dy = 1;
    double albedo = 0.2;
    auto mgrid = rayli::vgrid::MysticCloud(dx, dy, nx, ny, zlev);
    auto grid = rayli::vgrid::Cyclic(mgrid, {{0, 0, -std::numeric_limits<double>::infinity()},{nx*dx, ny*dy, std::numeric_limits<double>::infinity()}});
    //auto grid = mgrid;
    std::cout << "\nGrid parameters:\n" << "dx = "  << dx << " dy = " << dy << " nx = " << nx << " ny = " << ny << " nlyr = " << nlyr << " nlev = " << nlev << " nmu = " << nmu << " nphi = " << nphi <<"\n\n";

    // Camera parameters and camera declaration
    size_t Nxpixel = 90;
    size_t Nypixel = 1;
    double fov = 2;
    double fov_phi1 = -45;
    double fov_phi2 = 45;
    double fov_theta1 = 90-45;
    double fov_theta2 = 90+45;
    double xloc = 4; //in km
    double yloc = 0.01;
    double zloc = 2;
    auto loc = Eigen::Vector3d{xloc,yloc,zloc};     
    double fovx = fov;
    double fovy = fov * Nxpixel / Nypixel;
    size_t rays = Nxpixel*Nypixel;
    auto cam = MysticPanoramaCamera(loc, 0, 0, fov_phi1, fov_phi2, fov_theta1, fov_theta2, rays); 

    // Stream and substream calculations
    size_t nsub = 30; //50;
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
    
    calc_image(grid, cam, Nxpixel, Nypixel, dx, dy, zlev, kext,  g1, w0, albedo, muEdir, nx, ny, nlyr, nmu, nphi, mus, phis, wmus,  wphis,
                sza_dir, Edir, radiances, streams, nsub, substreams, image, opthick_image, Ldiff_i, Lup_i, Ldown_i, Ldir_i, groundbox, gRdir_i, gRdiff_i);

    
    std::cout << "Stop Raytracing...\n";    
    
    print_parameters(radfpath, flxfpath, opfpath, outputfname, dx, dy, albedo, xloc, yloc, zloc, Nxpixel, Nypixel, fov, fov_phi1, fov_phi2, fov_theta1,
                        fov_theta2, rays, nsub);


    {
        using namespace netCDF;    
        NcFile file(outputfname + ".nc", NcFile::FileMode::replace);
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
       }


}
