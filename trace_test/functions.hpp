#include <vgrid/MysticCloud.hpp>
#include <netcdf>
#include <iostream>
#include <array>
#include <camera.hpp>
#include <vgrid/Cyclic.hpp>
#include <cmath>
#include <numeric>
#include <Eigen/Dense>
#include <stdexcept>

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
    
    size_t out = 0;

    for(size_t i = 0; i<N; ++i) {
        if(index[i] >= shape[i]){
            throw std::out_of_range(std::string("given index is out of bounds: ") +
                        "Indexvalue " + std::to_string(index[i]) + 
                        " exceeds shape " + std::to_string(shape[i]) +
                        " in dimension " + std::to_string(i)); 
        }
        out *= shape[i];
        out += index[i];
    }
    return out;
}

Eigen::Vector3d angleToVec(double mu, double phi) {
    if(mu>1.01){std::cout << "mu to big:  " << mu << "\n"; std::abort();}
    if(mu>1){mu=1;}
    double a = sin(acos(mu));
    double x = a*cos(phi);
    double y = a*sin(phi);
    double z = mu;
    
    if(std::isnan(a) ==1){ std::cout << mu << "  " << "\n";}

    return Eigen::Vector3d{x,y,z};
}


std::array<double,2> groundReflection_lambert(auto ray, size_t idx, double albedo, size_t nx, size_t ny, size_t nlyr, size_t nmu, size_t nphi, std::vector<double> mu, std::vector<double> wmu, std::vector<double> wphi, const std::vector<double>& Edir, const std::vector<double>& rad) {    
    double E = Edir[idx] * albedo / M_PI;
    auto [x,y,z] = indexDecompose(idx, std::array<size_t,3>{nx,ny,nlyr+1});  
    //std::cout << x << "  " << y << "  " << z << "  " << Edir[idx]  << "  " << albedo << "  " << M_PI << "  " << E  << "\n"; 
    double irradiance = 0;
    for(size_t i = 0; i<nmu; ++i) {
        for(size_t j = 0; j<nphi; ++j) {
            size_t index = indexRecompose(std::array<size_t,5>{x,y,z,i,j},std::array<size_t,5>{nx,ny,nlyr+1,nmu,nphi});
            if(mu[i] < 0) {
                irradiance += rad[index] * fabs(mu[i])  * wmu[i] * wphi[j];	
            }
        }
    }
    E += irradiance * albedo / M_PI;
    return std::array<double,2>{E,irradiance * albedo / M_PI};
}

double phase_HG(double g, double mu) {
    return 1./(4*M_PI) * (1-g*g)/pow(1+g*g-2*g*mu,3./2.);
}


double calc_pHG(double wmu_s, double wmu_e, double phi, double wphi, auto ray, double g, double n=50){
    double delmu = (wmu_e-wmu_s)/n;
    double delphi = wphi/n;
    double pf = 0;
    double weightsum = 0;
    double wphi_s = phi - wphi/2;
    double wphi_e = phi + wphi/2;

    

    for(size_t i = 0; i<n; ++i){
        double m = wmu_s + (i+0.5) * delmu;
        for(size_t j = 0; j<n; ++j){
            double p = wphi_s + (j+0.5) * delphi;
            double sa = (-ray.d).dot(angleToVec(m,p));
            pf += phase_HG(g, sa);
            weightsum += 1;
            if(pf < 0){std::cout << " negative \n";}
            //std::cout << wmu_s << "  " << wmu_e << "  " << m << "  " << p << "  "  << g << "  " << sa << "  " << delmu << "  " << delphi << " " << pf << "\n";
        }
    }
    //std::cout << weightsum << " " << pf << "\n";
    return pf/(n*n);
}


std::array<double,3> calc_Ldiff(auto ray, double tfar, double tnear, size_t idx, double g1, size_t nx, size_t ny, size_t nlyr, size_t nmu, size_t nphi,
        std::vector<double> mu, std::vector<double> phi, std::vector<double> wmu, std::vector<double> wphi, std::vector<double> rad) {
    double Lup = 0;
    double Ldown = 0;
    auto [x,y,z] = indexDecompose<3>(idx, {nx,ny,nlyr});
    //std::cout << "in func: " << x << " " << y << " " << z << "\n";
    double wmu_s = -1;
    double wmu_e,wphi_e;
    double pf,pf_n,pf_i;

    for(size_t i = 0; i < nmu; ++i) {
        wmu_e = wmu_s + wmu[i];
        for(size_t j = 0; j < nphi; ++j) {
            Eigen::Vector3d lvec = angleToVec(mu[i],phi[i]);
            auto muscatter = (-ray.d).dot(lvec); 
            double weight = wmu[i] * wphi[j];
            //std::cout << "Part 1\n";
            size_t index_t = indexRecompose(std::array<size_t,5>{x,y,z+1,i,j},std::array<size_t,5>{nx,ny,nlyr+1,nmu,nphi});
            size_t index_b = indexRecompose(std::array<size_t,5>{x,y,z,i,j},std::array<size_t,5>{nx,ny,nlyr+1,nmu,nphi});
            if(0==0){
                pf = weight * phase_HG(g1, muscatter);
            }
            else{
                pf = weight * calc_pHG(wmu_s, wmu_e, phi[j], wphi[j], ray, g1, 1);
            }
            //pf_n = weight * phase_HG(g1, muscatter);
            //pf_i = calc_pHG(wmu_s, wmu_e, wphi_s, wphi_e, ray, g1);
            if(mu[i] < 0){
                Ldown +=  rad[index_t] * pf;
                //if(z < 25){
                //    auto [a,b,c,d,e] = indexDecompose<5>(index_t, {nx,ny,nlyr+1,nmu,nphi});
                //    std::cout << "down: "  << a << " " << b << " " << c << " " << d << " " << e <<  "\n";
                //    std::cout << rad[index_t] << "     " << pf << "\n";
                //}
            }
            if(mu[i] > 0){
                Lup += rad[index_b] * pf;
                //if(z < 25){
                //    auto [a,b,c,d,e] = indexDecompose<5>(index_b, {nx,ny,nlyr+1,nmu,nphi});
                //    std::cout << "up: " << a << " " << b << " " << c << " " << d << " " << e <<  "\n";
                //    std::cout << rad[index_b] << "     " << pf << "\n";
                //}
            }

            //std::cout << phase_HG(g1, muscatter) << " " << muscatter  << "  mu = " << mu[i] << "  phi = " << phi[j] << "\n";
            //if(g1 != 0){
            //    std::cout << "pf_n = " << pf_n << "  pf_i = " << pf_i << "  wmu = " << wmu[i] << "  wphi = " << wphi[i] << "  pf_i/pf_n = " << pf_i/pf_n << "\n" ;
            //}
        }
        wmu_s = wmu_e;
    }
    
    return std::array<double,3>{Lup + Ldown,Lup,Ldown}; 
}


