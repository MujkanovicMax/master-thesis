#include <vgrid/MysticCloud.hpp>
#include <netcdf>
#include <iostream>
#include <array>
#include <camera.hpp>
#include <vgrid/Cyclic.hpp>
#include <cmath>
#include <numeric>
#include <Eigen/Dense>

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

Eigen::Vector3d angleToVec(double mu, double phi) {
	double a=sin(acos(mu));
	double x = a*cos(phi);
	double y = a*sin(phi);
	double z = mu;

	return Eigen::Vector3d{x,y,z};
}


double groundReflection_lambert(auto ray, size_t idx, double albedo, size_t nx, size_t ny, size_t nlyr, size_t nmu, size_t nphi, std::vector<double> mu, std::vector<double> wmu, std::vector<double> wphi, const std::vector<double>& Edir, const std::vector<double>& rad) {
	double E = Edir[idx] * albedo / M_PI;
	auto [x,y,z] = indexDecompose(idx, std::array<size_t,3>{nx,ny,nlyr+1});
	double irradiance = 0;
	for(size_t i = 0; i<nmu; ++i) {
		for(size_t j = 0; j<nphi; ++j) {
			size_t index = indexRecompose(std::array<size_t,5>{x,y,z,i,j},std::array<size_t,5>{nx,ny,nlyr+1,nmu,nphi});
		        irradiance += rad[index] * abs(mu[i])  * wmu[i] * wphi[j];	
		}
	}
	E += irradiance * albedo / M_PI;
	return E;
}

double phase_HG(double g, double mu) {
	return 1./(4*M_PI) * (1-g*g)/pow(1+g*g-2*g*mu,3./2.);
}



double calc_Ldiff(auto ray, double tfar, double tnear, size_t idx, double g1, size_t nx, size_t ny, size_t nlyr, size_t nmu, size_t nphi,
	       	  std::vector<double> mu, std::vector<double> phi, std::vector<double> wmu, std::vector<double> wphi, std::vector<double> rad) {
	double L = 0;
	auto [x,y,z] = indexDecompose(idx, std::array<size_t,3>{nx,ny,nlyr+1});
	//std::cout << "nmu = " << nmu << "  nphi = " << nphi << "\n";
	for(size_t i = 0; i < nmu; ++i) {
		for(size_t j = 0; j < nphi; ++j) {
			//std::cout << "Part 1\n";
			Eigen::Vector3d lvec = angleToVec(mu[i],phi[i]);
			size_t index = indexRecompose(std::array<size_t,5>{x,y,z+1,i,j},std::array<size_t,5>{nx,ny,nlyr+1,nmu,nphi});
			if( lvec.dot(ray((tfar+tnear)/2)-ray(tnear)) >= 0) { 
				L +=  abs(mu[i]) * rad[index] * wmu[i] * wphi[i] * phase_HG(g1, ray.d.dot(lvec));
			}
			index = indexRecompose(std::array<size_t,5>{x,y,z,i,j},std::array<size_t,5>{nx,ny,nlyr+1,nmu,nphi});
			if( lvec.dot(ray((tfar+tnear)/2)-ray(tfar)) >= 0) { 
				L +=  abs(mu[i]) * rad[index] * wmu[i] * wphi[i] * phase_HG(g1, ray.d.dot(lvec));
			}
			//std::cout << "PArt 2\n";
		}
	}
	
	return L;
}

int main(int argc, char** argv) {
    std::cout << "trace optical thickness\n";

    std::vector<double> radiances;
    std::vector<double> mus;
    std::vector<double> phis;
    std::vector<double> wmus;
    std::vector<double> wphis;
    size_t nmu,nphi;

    {
        using namespace netCDF;

        NcFile file("radiances.nc", NcFile::FileMode::read);
	size_t nx = file.getDim("x").getSize(); 
	size_t ny = file.getDim("y").getSize(); 
	size_t nz = file.getDim("z").getSize(); 
	nphi = file.getDim("phi").getSize(); 
	nmu = file.getDim("mu").getSize(); 
        size_t nwvl  = file.getDim("wvl").getSize();

	radiances.resize(nx*ny*nz*nphi*nmu*nwvl);
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
    Eigen::Vector3d sza_dir;
    double muEdir;

    {
	   using namespace netCDF;

	   NcFile file("Edir.nc", NcFile::FileMode::read);
   	   size_t Nx = file.getDim("x").getSize(); 
	   size_t Ny = file.getDim("y").getSize(); 
	   size_t Nz = file.getDim("z").getSize(); 
	   size_t Nwvl = file.getDim("wvl").getSize();

	   file.getAtt("mu0").getValues(&muEdir);
	   sza_dir[0] = sin(acos(muEdir));
	   sza_dir[1] = 0;
	   sza_dir[2] = muEdir;

	   Edir.resize(Nx*Ny*Nz*Nwvl);
	   file.getVar("Edir").getVar(Edir.data());
	   
	
    }




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

    double albedo = 0.2;
    

    auto cam = SimpleCamera(loc, lookat({0,0,-1},{0,1,0}), fovx, fovy, rays);
    double transmission;
    

    std::vector<double> image(Nxpixel * Nypixel);
    std::vector<double> opthick_image(Nxpixel * Nypixel);
    double optical_thickness;
    double radiance;
    size_t I=0;
    for(size_t i = 0; i < Nypixel; ++i) {
        double ypx = (i + 0.5) / Nypixel;
        for(size_t j = 0; j < Nxpixel; ++j) {
            double xpx = (j + 0.5) / Nxpixel;
            auto ray = cam.compute_ray(Eigen::Vector2d{xpx, ypx});
	    //auto ray = Ray{loc, {0, 0, -1}};
	    optical_thickness = 0;
	    radiance = 0;
	    size_t groundidx;
	    std::cout << I << "\n";
	    I++;
            for(auto slice: grid.walk_along(ray, 0., std::numeric_limits<double>::infinity())) {
		    if(auto pvol = std::get_if<VolumeSlice>(&slice)) {
                        auto [x,y,z] = indexDecompose<3>(pvol->idx, {nx,ny,nlyr});
			size_t optprop_index = indexRecompose(std::array{z,x,y},std::array{nlyr,nx,ny});//((n_index[2] * nx) + n_index[0]) * ny + n_index[1];
			size_t rad_index = indexRecompose(std::array{x,y,z+1},std::array{nx,ny,nlyr+1});//((n_index[0] * ny) + n_index[1]) * (nlyr+1) + n_index[2];
		       	groundidx = indexRecompose(std::array{x,y,z},std::array{nx,ny,nlyr+1});
			double L = Edir[rad_index];
			transmission = exp(-optical_thickness);
			double dtau = (pvol->tfar - pvol->tnear) * kext[optprop_index];
			optical_thickness += dtau;
                        double phase_function = phase_HG(g1[optprop_index], ray.d.dot(sza_dir));
			double Ldiff = calc_Ldiff(ray,pvol->tfar, pvol->tnear, optprop_index, g1[optprop_index], nx, ny, nlyr, nmu, nphi, mus, phis, wmus, wphis, radiances); 
			radiance += transmission * L * w0[optprop_index] * dtau * phase_function + transmission * Ldiff * w0[optprop_index] * dtau;  
//			radiance += transmission * sum_over_radiances_around_box(phase_function * incoming_radiance);
//		        std::cout << "VolIdx = " << pvol->idx << "  L= " << L << "  x= " << x << "  y= " << y << "  z= " << z << " radidx = " << rad_index << " optpropidx = " << optprop_index << " w0=" << w0[optprop_index] << "  dtau=" << dtau << "  radiance=" << radiance << " kext=" << kext[optprop_index] << "\n";   
		    }
	            
	    }

            double ground_E = groundReflection_lambert(ray,groundidx,albedo, nx, ny, nlyr, nmu, nphi, mus, wmus, wphis, Edir, radiances) * exp(-optical_thickness);            
            image[j + i * Nxpixel] = radiance + ground_E;
	    opthick_image[j+i*Nxpixel] = optical_thickness;
	}
    }

    {
	using namespace netCDF;    
        NcFile file("output.nc", NcFile::FileMode::replace);
        auto xdim = file.addDim("x", Nxpixel);
        auto ydim = file.addDim("y", Nypixel);
        file.addVar("image", NcType::nc_DOUBLE, {ydim, xdim}).putVar(image.data());
        file.addVar("optical_thickness_image", NcType::nc_DOUBLE, {ydim, xdim}).putVar(opthick_image.data());

    }



//    for(auto slice: grid.walk_along(ray, 0., std::numeric_limits<double>::infinity())) {
//        if(auto pvol = std::get_if<VolumeSlice>(&slice)) {
//            std::cout << "passing volume " << pvol->idx << " from " << pvol->tnear << " to " << pvol->tfar << "\n";
//        }
//    }
}
