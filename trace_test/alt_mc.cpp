#include <vgrid/MysticCloud.hpp>
#include <netcdf>
#include <iostream>
#include <array>
#include <camera.hpp>
#include <vgrid/Cyclic.hpp>
#include <cmath>
#include <numeric>
#include <Eigen/Dense>
#include <progress.hpp>
#include <recorder.hpp>
#include <FarSink.hpp>
#include <PathTracer.hpp>
#include <atm.hpp>
#include <material.hpp>
#include <texture/texture.hpp>
#include <vgrid/AASurfVGrid.hpp>
#include <vgrid/AddVGrid.hpp>
#include "functions.hpp"

struct BSurfaceIndexer {
    using a_type = size_t;
    using b_type = size_t;
    using value_type = size_t;

    size_t surface_a(size_t ia) const {
        return 0;
    }
    size_t surface_b(size_t ib) const {
        return ib;
    }
};

struct AVolumeIndexer {
    using a_type = size_t;
    using b_type = size_t;
    using value_type = size_t;

    size_t volume_a(size_t ia) const {
        return ia;
    }
    size_t volume_b(size_t ib) const {
        return 0;
    }
    size_t volume_ab(size_t ia, size_t ib) const {
        return 0;
    }
};


int main(int argc, char** argv) {

    std::vector<double> kext;
    std::vector<double> zlev;
    std::vector<double> w0,g1;
    size_t nlev, nlyr, nx, ny;
    double dx = 0.1;
    double dy = 1.;

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
   
    Eigen::Vector3d sza_dir;
    double muEdir;

    {
        using namespace netCDF;

        NcFile file("Edir_t.nc", NcFile::FileMode::read);
        file.getAtt("mu0").getValues(&muEdir);

        sza_dir[0] = -sin(acos(muEdir));
        sza_dir[1] = 0;
        sza_dir[2] = muEdir;
    }


    auto grid2optprop = [&](size_t i){
        auto [x,y,z] = indexDecompose<3>(i,{nx,ny,nlyr});
        return indexRecompose(std::array{z,x,y},std::array{nlyr,nx,ny});
    };

    auto k_ext_tex =
        function_texture<size_t>([&](size_t i) { 
                return kext[grid2optprop(i)]; });


    auto hg_tex =
        function_texture<size_t>([&](size_t i) { 
                return HGScat<std::mt19937_64>(g1[grid2optprop(i)]); });

    auto w_tex =
        function_texture<size_t>([&](size_t i) { 
                return w0[grid2optprop(i)]; });
    
    auto scat = hg_tex * k_ext_tex * w_tex;
    
    auto mgrid = rayli::vgrid::MysticCloud(dx, dy, nx, ny, zlev);
    auto grid = rayli::vgrid::Cyclic(mgrid, {{0, 0, -std::numeric_limits<double>::infinity()},{nx*dx, ny*dy, std::numeric_limits<double>::infinity()}});
    auto ground = rayli::vgrid::AASurfVGrid{0, 2};
    auto ground_cloud_grid = rayli::vgrid::AddVGrid(grid,ground, BSurfaceIndexer{}, AVolumeIndexer{}); 
    
    
    size_t Nxpixel = 90;
    size_t Nypixel = 1;
    double fov = 2;

    auto loc = Eigen::Vector3d{3,0.01,2};
    double fovx = fov;
    double fovy = fov * Nxpixel / Nypixel;
    size_t rays = 100000;

    double albedo = 0.2;
    auto surf = constant_texture<size_t>(LambertianMaterial{albedo});
    auto sink = FarSink<std::mt19937_64>{-sza_dir.normalized()};
    auto progress = ProgressBar(rays);

    auto cam = MysticPanoramaCamera(loc, 0, 0, -45, 45, 90, 90, rays, {}, progress.reporter());
    auto scn = Scn{ground_cloud_grid, k_ext_tex, scat, surf, sink};

    auto tracer = PathTracer{scn};
    auto recorder = ImageRecorder{Nxpixel,Nypixel};


    progress.display_while([&]() { tracer.solve_par(cam, recorder); });
    
    {

        using namespace netCDF;

        NcFile file("rayli_panorama.nc", NcFile::FileMode::replace);
        auto xdim = file.addDim("x", recorder.sx);
        auto ydim = file.addDim("y", recorder.sy);
        auto datavar = file.addVar("data", NcType::nc_DOUBLE, {xdim, ydim});
        auto countvar = file.addVar("count", NcType::nc_UINT64, {xdim, ydim});
        auto radvar =
          file.addVar("transmissivity", NcType::nc_DOUBLE, {xdim, ydim});
        auto stdvar =
          file.addVar("transmissivity_std", NcType::nc_DOUBLE, {xdim, ydim});
        datavar.putVar(recorder.data.data());
        countvar.putVar(recorder.count.data());
        auto rad = recorder.radiance();
        radvar.putVar(rad.data());
        auto rad_std = recorder.radiance_std();
        stdvar.putVar(rad_std.data());
    }

} 
