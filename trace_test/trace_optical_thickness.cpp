#include <vgrid/MysticCloud.hpp>
#include <netcdf>
#include <iostream>
#include <camera.hpp>

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

    std::vector<double> oprop;
    std::vector<double> zlev;
    size_t nlev, nlyr, nx, ny;
    double dx = 1.;
    double dy = 1.;

    { 
	    using namespace netCDF;

	    // key: caoth3d_0_wc_ext
	    NcFile file("test.optical_properties.nc", NcFile::FileMode::read);
	    nlev = file.getDim("nlev").getSize();
	    nlyr = file.getDim("caoth3d_0_wc_nlyr").getSize();
	    nx = file.getDim("caoth3d_0_wc_Nx").getSize();
	    ny = file.getDim("caoth3d_0_wc_Ny").getSize();

	    zlev.resize(nlev);
	    oprop.resize(nlyr*nx*ny);
	    file.getVar("caoth3d_0_wc_ext").getVar(oprop.data());
            file.getVar("output_mc.z").getVar(zlev.data());
    }



    

    auto grid = rayli::vgrid::MysticCloud(dx, dy, nx, ny, zlev);

    
    

    size_t Nxpixel = 400;
    size_t Nypixel = 400;
    double fov = 2;

    auto loc = Eigen::Vector3d{3.5,0.5,10};
    double fovx = fov;
    double fovy = fov * Nxpixel / Nypixel;
    size_t rays = 90;
    

    auto cam = SimpleCamera(loc, lookat({0,0,-1},{0,1,0}), fovx, fovy, rays);


    std::vector<double> image(Nxpixel * Nypixel);
    double optical_thickness; 
    for(size_t i = 0; i < Nypixel; ++i) {
        double ypx = (i + 0.5) / Nypixel;
        for(size_t j = 0; j < Nxpixel; ++j) {
            double xpx = (j + 0.5) / Nxpixel;
            auto ray = cam.compute_ray(Eigen::Vector2d{xpx, ypx});
	    optical_thickness = 0;
            for(auto slice: grid.walk_along(ray, 0., std::numeric_limits<double>::infinity())) {
		    if(auto pvol = std::get_if<VolumeSlice>(&slice)) {
                        size_t z_index = pvol->idx % nlyr;
			size_t y_index = (pvol->idx / nlyr) % ny;
			size_t x_index = (pvol->idx / (nlyr * ny)) % nx;

			size_t optprop_index = ((z_index * nx) + x_index) * ny + y_index;
			    optical_thickness += (pvol->tfar - pvol->tnear) * oprop[optprop_index];
			radiance += exp(-optical_thickness) * sum_over_radiances_around_box(phase_function * incoming_radiance);
		    }
			    
	    }

            image[j + i * Nxpixel] = optical_thickness;
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
