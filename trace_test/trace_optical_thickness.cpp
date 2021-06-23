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

    int zdivs[9] {1,2,5,10,25,50}; // for level benchmark
    std::string names[4] {"irr_from_10s_twostr_only","irr_from_mystic_mcipa","irr_from_mystic","irr_from_10s_base"};
    int adivs[5] {2,4,8,16,32} ;
    for(int zdiv = 0; zdiv <= 0; ++zdiv){
        for(int adiv = 0; adiv <= 0; ++adiv){

            // Filenames for Input

            //std::string radfpath = "rad_zdiv" + std::to_string(zdivs[zdiv])  + "_mu_" + "2"  + "_phi_" + "1" + ".nc";
            //std::string radfpath = names[in] + "_div" + std::to_string(zdivs[zdiv]) +  ".nc";  //"irr_from32x32_myst_zdiv" + std::to_string(zdivs[zdiv])  + ".nc";
            //std::string radfpath = "../radiances/rad_testcases/checkerboard_full/rad_checkerboard.nc";
            std::string radfpath = "/project/meteo-scratch/Mujkanovic.Max/rad_philipp_new/rad_philipp_new.nc";

            //std::string flxfpath = "../radiances/rad_testcases/checkerboard_full/flx_checkerboard.nc";
            //std::string flxfpath = "flx_levels_div" + std::to_string(zdivs[zdiv]) + ".nc";
            std::string flxfpath = "/project/meteo-scratch/Mujkanovic.Max/rad_philipp_new/flx_philipp_new.nc";

            //std::string opfpath = "op_levels_div" + std::to_string(zdivs[zdiv]) + ".nc";
            //std::string opfpath = "../radiances/rad_testcases/checkerboard_full/op_checkerboard.nc";
            std::string opfpath = "/project/meteo-scratch/Mujkanovic.Max/rad_philipp_new/op_philipp_new_expanded_nobs.nc";

            // Filename for output
            std::string outputfname = "output2d_philipp_test";

            // Grid Parameters 
            double dx;
            double dy;
            double albedo;

            // Camera parameters and camera declaration
            int Nxpixel;
            int Nypixel;
            double fov = 2;
            double fov_phi1;
            double fov_phi2;
            double fov_theta1;
            double fov_theta2;
            double xloc; //in km
            double yloc;
            double zloc;
            

            //Subdivisions for Substreams
            int nsub; //50;
            

            //Parsing config file
            std::string configname = "config.txt";
            configparser(configname, radfpath, flxfpath, opfpath, outputfname, dx, dy, albedo, Nxpixel, Nypixel, fov_phi1, fov_phi2, fov_theta1, fov_theta2, xloc, yloc, zloc, nsub);

            auto loc = Eigen::Vector3d{xloc,yloc,zloc};     
            double fovx = fov;
            double fovy = fov * Nxpixel / Nypixel;
            size_t rays = Nxpixel*Nypixel;
            auto cam = MysticPanoramaCamera(loc, 0, 0, fov_phi1, fov_phi2, fov_theta1, fov_theta2, rays); 
            std::cout << "Reading netcdf data" << "\n\n";
            std::cout << "radfpath = " << radfpath << "\n";
            std::cout << "flxfpath = " << flxfpath << "\n";
            std::cout << "opfpath = " << opfpath << "\n\n";

            //NetCDF declarations
            std::vector<double> radiances;
            std::vector<double> mus;
            std::vector<double> phis;
            std::vector<double> wmus;
            std::vector<double> wphis;
            std::vector<double> Edir;
            std::vector<double> Edown;
            std::vector<double> kext;
            std::vector<double> zlev;
            std::vector<double> w0,g1;
            Eigen::Vector3d sza_dir;
            size_t nlev, nlyr, nx, ny;
            size_t nmu,nphi;
            double muEdir;

            //read NetCDF data
            read_radiances(radfpath, radiances, mus, phis, wmus, wphis,  nmu, nphi);
            read_flx(flxfpath, Edir, Edown, sza_dir, muEdir);
            read_opprop(opfpath, kext, zlev, w0, g1, nlev, nlyr, nx, ny);
            for(auto& ke: kext) {
                ke *= 1000;
            }

            std::cout << "done" << "\n";



            // Stream and substream calculations
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


            //Grid declaration
            auto mgrid = rayli::vgrid::MysticCloud(dx, dy, nx, ny, zlev);
            auto grid = rayli::vgrid::Cyclic(mgrid, {{0, 0, -std::numeric_limits<double>::infinity()},{nx*dx, ny*dy, std::numeric_limits<double>::infinity()}});
            std::cout << "wah\n";
            auto ray = cam.compute_ray(Eigen::Vector2d{0, 0});
            std::cout << "wah2\n";
            grid.walk_along(ray, 0., std::numeric_limits<double>::infinity());

            //outputting parameters 
            std::cout << "\nCamera parameters:\n\n" << "X Pixel: " << Nxpixel << "    Y Pixel: " << Nypixel << "\n";
            std::cout << "Camera at x = " << xloc << "km  y = " << yloc << "km  z = " << zloc << "km\n";
            std::cout << "Camera opening angles: Phi1 = " + std::to_string(fov_phi1) + "    Phi2 = " + std::to_string(fov_phi2) + "     Theta1 = " + std::to_string(fov_theta1) + "     Theta2 = " + \
                        std::to_string(fov_theta2) + "\n";
            std::cout << "\nGrid parameters:\n\n" << "dx = "  << dx << " dy = " << dy << " nx = " << nx << " ny = " << ny << " nlyr = " << nlyr 
                << " nlev = " << nlev << " nmu = " << nmu << " nphi = " << nphi <<"\n";
            std::cout << "Calculating phase functions with " + std::to_string(nsub) + " subdivisions.\n\n";

            //main image calculation
            calc_image(grid, cam, Nxpixel, Nypixel, dx, dy, zlev, kext,  g1, w0, albedo, muEdir, nx, ny, nlyr, nmu, nphi, mus, phis, wmus,  wphis,
                    sza_dir, Edir, radiances, streams, nsub, substreams, image, opthick_image, Ldiff_i, Lup_i, Ldown_i, Ldir_i, groundbox, gRdir_i, gRdiff_i);

            //printing calculation parameters 
            //print_parameters(radfpath, flxfpath, opfpath, outputfname, dx, dy, albedo, xloc, yloc, zloc, Nxpixel, Nypixel, fov, fov_phi1, fov_phi2, fov_theta1,
              //      fov_theta2, rays, nsub);



            //saving image data to NetCDF
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
