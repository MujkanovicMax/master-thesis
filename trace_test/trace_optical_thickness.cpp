#include <vgrid/MysticCloud.hpp>
#include <netcdf>
#include <iostream>
#include <array>
#include <camera.hpp>
#include <vgrid/Cyclic.hpp>
#include <cmath>

template <size_t N>
constexpr std::array<size_t,N> indexDecompose(size_t index, const std::array<size_t,N>& shape) {
	std::array<size_t,N> out = {};
	for(size_t i = N; i>0; --i) {
		out[i-1] = index % shape[i-1];
		index /= shape[i-1];	
	}
	return out;
}	

template <size_t N>
constexpr size_t indexRecompose(std::array<size_t,N> index, const std::array<size_t,N>& shape) {
	size_t out = 0;//index[0] * shape[1]
	for(size_t i = 0; i<N-1; ++i) {
		out += index[i];
		out *= shape[i+1];
	}
	return out + index[N-1];
}

int main(int argc, char** argv) {
    std::cout << "trace optical thickness\n";

    std::vector<double> radiances;


    {
        using namespace netCDF;

        NcFile file("radiances.nc", NcFile::FileMode::read);
	size_t Nx = file.getDim("x").getSize(); 
	size_t Ny = file.getDim("y").getSize(); 
	size_t Nz = file.getDim("z").getSize(); 
	size_t Nphi = file.getDim("phi").getSize(); 
	size_t Nmu = file.getDim("mu").getSize(); 
        size_t Nwvl  = file.getDim("wvl").getSize();

	radiances.resize(Nx*Ny*Nz*Nphi*Nmu*Nwvl);

        file.getVar("radiance").getVar(radiances.data());
    }

    std::vector<double> Edir;
   // double muEdir;

    {
	   using namespace netCDF;

	   NcFile file("Edir.nc", NcFile::FileMode::read);
   	   size_t Nx = file.getDim("x").getSize(); 
	   size_t Ny = file.getDim("y").getSize(); 
	   size_t Nz = file.getDim("z").getSize(); 
	   size_t Nwvl = file.getDim("wvl").getSize();

//	   muEdir = file.getAtt("mu").getValues();

	   Edir.resize(Nx*Ny*Nz*Nwvl);
	   file.getVar("Edir").getVar(Edir.data());
	   
	
    }




    std::vector<double> kext;
    std::vector<double> zlev;
    std::vector<double> w0;
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
	    file.getVar("caoth3d_0_wc_ext").getVar(kext.data());
            file.getVar("output_mc.z").getVar(zlev.data());
	    file.getVar("caoth3d_0_wc_ssa").getVar(w0.data());

    }



    auto mgrid = rayli::vgrid::MysticCloud(dx, dy, nx, ny, zlev);
    auto grid = rayli::vgrid::Cyclic(mgrid, {{0, 0, -std::numeric_limits<double>::infinity()},{nx*dx, ny*dy, std::numeric_limits<double>::infinity()}});
    
    std::cout <<  dx << " " << dy << " " << " " << nx << " " << ny << "\n";

    size_t Nxpixel = 400;
    size_t Nypixel = 400;
    double fov = 2;

    auto loc = Eigen::Vector3d{4,0.5,10};
    double fovx = fov;
    double fovy = fov * Nxpixel / Nypixel;
    size_t rays = 90;
    

    auto cam = SimpleCamera(loc, lookat({0,0,-1},{0,1,0}), fovx, fovy, rays);
    double transmission;
    

    std::vector<double> image(Nxpixel * Nypixel);
    double optical_thickness;
    double radiance; 
    for(size_t i = 0; i < Nypixel; ++i) {
        double ypx = (i + 0.5) / Nypixel;
        for(size_t j = 0; j < Nxpixel; ++j) {
            double xpx = (j + 0.5) / Nxpixel;
            auto ray = cam.compute_ray(Eigen::Vector2d{xpx, ypx});
	    //auto ray = Ray{loc, {0, 0, -1}};
	    optical_thickness = 0;
	    radiance = 0;
	    
            for(auto slice: grid.walk_along(ray, 0., std::numeric_limits<double>::infinity())) {
		    if(auto pvol = std::get_if<VolumeSlice>(&slice)) {
                        auto [x,y,z] = indexDecompose<3>(pvol->idx, {nx,ny,nlyr});
			size_t optprop_index = indexRecompose(std::array{z,x,y},std::array{nlyr,nx,ny});//((n_index[2] * nx) + n_index[0]) * ny + n_index[1];
			size_t rad_index = indexRecompose(std::array{x,y,z+1},std::array{nx,ny,nlyr+1});//((n_index[0] * ny) + n_index[1]) * (nlyr+1) + n_index[2];
		       	double L = Edir[rad_index];
			transmission = exp(-optical_thickness);
			double dtau = (pvol->tfar - pvol->tnear) * kext[optprop_index];
			optical_thickness += dtau;
                        double phase_function = 1./(4*M_PI);
			radiance += transmission * L * w0[optprop_index] * dtau * phase_function;  
//			radiance += transmission * sum_over_radiances_around_box(phase_function * incoming_radiance);
//		        std::cout << "VolIdx = " << pvol->idx << "  L= " << L << "  x= " << x << "  y= " << y << "  z= " << z << " radidx = " << rad_index << " optpropidx = " << optprop_index << " w0=" << w0[optprop_index] << "  dtau=" << dtau << "  radiance=" << radiance << " kext=" << kext[optprop_index] << "\n";   
		    }
			    
	    }

            image[j + i * Nxpixel] = radiance;
	}
    }

    {
	using namespace netCDF;    
        NcFile file("output.nc", NcFile::FileMode::replace);
        auto xdim = file.addDim("x", Nxpixel);
        auto ydim = file.addDim("y", Nypixel);
        file.addVar("image", NcType::nc_DOUBLE, {ydim, xdim}).putVar(image.data());

    }



//    for(auto slice: grid.walk_along(ray, 0., std::numeric_limits<double>::infinity())) {
//        if(auto pvol = std::get_if<VolumeSlice>(&slice)) {
//            std::cout << "passing volume " << pvol->idx << " from " << pvol->tnear << " to " << pvol->tfar << "\n";
//        }
//    }
}
