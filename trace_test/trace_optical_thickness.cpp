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
    std::string radfpath = "irr_from_10s_mcipa_pp.nc";  //"irr_from32x32_myst_zdiv" + std::to_string(zdivs[zdiv])  + ".nc";
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
   

    //std::string flxfpath = "/home/m/Mujkanovic.Max/ma/radiances/job_flx/mc.flx.spc.nc";
    std::string flxfpath = "flx_levels_div" + std::to_string(zdivs[zdiv]) + ".nc";
    //std::string flxfpath = "../radiances/flx_philipp_16x16.nc";
    std::cout << "flxfpath = " << flxfpath << "\n\n";
    std::vector<double> Edir;
    std::vector<double> Edown;
    Eigen::Vector3d sza_dir;
    double muEdir;
    
    read_flx( flxfpath, Edir, Edown, sza_dir, muEdir );
   

    std::string opfpath = "op_levels_div" + std::to_string(zdivs[zdiv]) + ".nc";
    //std::string opfpath = "test.optical_properties.nc";
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
    std::string outputfname = "2dtest_10s_ipa";



    // Grid Parameters and grid declaration
    double dx = 0.1;
    double dy = 1;
    double albedo = 0.2;
    auto mgrid = rayli::vgrid::MysticCloud(dx, dy, nx, ny, zlev);
    auto grid = rayli::vgrid::Cyclic(mgrid, {{0, 0, -std::numeric_limits<double>::infinity()},{nx*dx, ny*dy, std::numeric_limits<double>::infinity()}});
    //auto grid = mgrid;
    std::cout << "\nGrid parameters:\n" << "dx = "  << dx << " dy = " << dy << " nx = " << nx << " ny = " << ny << " nlyr = " << nlyr << " nlev = " << nlev << " nmu = " << nmu << " nphi = " << nphi <<"\n\n";

    // Camera parameters and camera declaration
    size_t Nxpixel = 1000;
    size_t Nypixel = 1000;
    double fov = 2;
    double fov_phi1 = -45;
    double fov_phi2 = 45;
    double fov_theta1 = 90-45;
    double fov_theta2 = 90+45;
    double xloc = 4; //in km
    double yloc = 0.01;
    double zloc = 2;
    auto loc = Eigen::Vector3d{xloc,yloc,zloc};     //def  Eigen::Vector3d{4,0.01,2}
    double fovx = fov;
    double fovy = fov * Nxpixel / Nypixel;
    size_t rays = 1000*1000;
    auto cam = MysticPanoramaCamera(loc, 0, 0, fov_phi1, fov_phi2, fov_theta1, fov_theta2, rays); //def: 0,0,-45,45,90,90,90

    // Stream and substream calculations
    size_t nsub = 10; //50;
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
    size_t a,b,c;
    //Main loop
    std::cout << "Starting ray tracing..." << "\n";
    
    for(size_t i = 0; i < Nypixel; ++i) {
        
        double ypx = (i + 0.5) / Nypixel;
        //for(size_t j = 0; j < 1; ++j) {
        
        for(size_t j = 0; j < Nxpixel; ++j) {
            //std::cout << "ray " << j << "\n";
            double xpx = (j + 0.5) / Nxpixel;
            
            auto ray = cam.compute_ray(Eigen::Vector2d{xpx, ypx});
            //std::cout << ray << "\n";
            //auto ray = Ray{loc, Eigen::Vector3d{-0.700909264299851,1.82977702939428e-16,-0.713250449154182}.normalized()};
            //auto ray = Ray{loc, Eigen::Vector3d{-0.0610485395348568,1.6787439687354e-16,-0.998134798421867}.normalized()};
            //sum zeroing
            double optical_thickness = 0;
            double radiance = 0;
            double Lups = 0;
            double Ldowns = 0;
            double Ldirs = 0;
            double Ldiffs = 0;
            size_t groundidx;

            for(auto slice: grid.walk_along(ray, 0., std::numeric_limits<double>::infinity())) {
                
                if(auto pvol = std::get_if<VolumeSlice>(&slice)) {
                    
                    //index calculations ( 3D -> 1D, 1D -> 3D )
                    auto [x,y,z] = indexDecompose<3>(pvol->idx, {nx,ny,nlyr});
                    a = x;
                    b = y;
                    c = z;
                    //std::cout << "x = " << x << "   y = " << y << "     z = " << z << "\n";
                    size_t optprop_index = indexRecompose(std::array{z,x,y},std::array{nlyr,nx,ny});
                    size_t rad_index = indexRecompose(std::array{x,y,z+1},std::array{nx,ny,nlyr+1});
                    groundidx = indexRecompose(std::array{x,y,z},std::array{nx,ny,nlyr+1});
                    
                    //main pixel calculation ( radiance summation, etc.)
                    double L = Edir[rad_index]/fabs(muEdir);
                    //std::cout << Edir[rad_index] << "  " << Edir[rad_index-1] << "  " << Edir[rad_index+1] << "     " << x << "     " << y << "     " << z+1 << "\n"; 
                    double transmission = exp(-optical_thickness);
                    double dtau = (pvol->tfar - pvol->tnear) * kext[optprop_index];
                    optical_thickness += dtau;
                    double phase_function = phase_HG(g1[optprop_index], (-ray.d).dot(sza_dir.normalized()));
                    //std::cout << "dtau = " << dtau << "\n";
                    auto [Lup_Plus_Ldown, Lup, Ldown] = calc_Ldiff(ray, dx, dy, zlev, pvol->tfar, pvol->tnear, pvol->idx, kext[optprop_index], dtau, g1[optprop_index],\
                                       nx, ny, nlyr, nmu, nphi, mus, phis, wmus, wphis, radiances, streams, nsub, substreams);
                    double scatter_prob = 1 - exp(-dtau*w0[optprop_index]);
                    //std::cout << " Ldiff = " << Lup + Ldown << "   Lup = " << Lup << "   Ldown = " << Ldown << "\n";
                    //std::cout << " Lup = " << Lup << " Lup_sum = " << Lups <<" transm = "<< transmission << " w0 = "\
                            << w0[optprop_index] << " dtau = " << dtau  << " scatter_prob = " << scatter_prob << "\n";
                    
                    // summation for pixelvalues        
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
            //bottom-most index ground reflection calculation
            auto [gR, gRdir, gRdiff] = groundReflection_lambert(ray,groundidx,albedo, nx, ny, nlyr, nmu, nphi, mus, wmus, wphis, Edir, radiances); 
            double ground_E = gR *  exp(-optical_thickness);            
            
            //std::cout << "Search Index x= " << a << " y= " << b << " z= " << "  RAY = " << ray <<"\n";
            
            //writing pixelvalues to image
            image[j + i * Nxpixel]      = radiance + ground_E;
            opthick_image[j+i*Nxpixel]  = optical_thickness;
            Ldiff_i[j+i*Nxpixel]        = Ldiffs;
            Lup_i[j+i*Nxpixel]          = Lups;
            Ldown_i[j+i*Nxpixel]        = Ldowns;
            Ldir_i[j+i*Nxpixel]         = Ldirs;
            groundbox[j+i*Nxpixel]      = indexDecompose<3>(groundidx,std::array<size_t,3>{nx,ny,nlyr})[0];
            gRdir_i[j+i*Nxpixel]        = gRdir * exp(-optical_thickness);
            gRdiff_i[j+i*Nxpixel]       = gRdiff * exp(-optical_thickness);
            //std::cout << "image val = " << image[j + i * Nxpixel] << "  Ldiff = " << Ldiff_i[j+i*Nxpixel] << "  Lup = " << Lup_i[j+i*Nxpixel] << "  Ldown = " << Ldown_i[j+i*Nxpixel] << "  Ldir = " << Ldir_i[j+i*Nxpixel] << "\n";
           
//            std::cout << "Pixel " << j+1 << " done" << "  ival = " << image[j+i*Nxpixel] << " "\
                <<  Ldiff_i[j+i*Nxpixel] << " "  << Lup_i[j+i*Nxpixel] << " " << Ldown_i[j+i*Nxpixel]<< " " << Ldir_i[j+i*Nxpixel] <<"\n";


        }
    }
    std::cout << "Stop Raytracing...\n";    
    
    print_parameters(radfpath, flxfpath, opfpath, outputfname, dx, dy, albedo, xloc, yloc, zloc, Nxpixel, Nypixel, fov, fov_phi1, fov_phi2, fov_theta1,
                        fov_theta2, rays, nsub);


    {
        using namespace netCDF;    
        //NcFile file("output_mu" + std::to_string(MU) + "_phi" + std::to_string(PHI) +  ".nc", NcFile::FileMode::replace);
        //NcFile file("output_zdiv" + std::to_string(zdivs[zdiv])  + "_mu_" + std::to_string(adivs[adiv])  + "_phi_" + std::to_string(adivs[adiv]) + ".nc", NcFile::FileMode::replace);
        //NcFile file("output_subsamp_plot_4x4_nsub_" + std::to_string(zdivs[zdiv])  + ".nc", NcFile::FileMode::replace);
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
