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
    double p = phi/180. * M_PI;
    double a = sin(acos(mu));
    double y = a*cos(p);
    double x = a*sin(p);
    double z = mu;
    
    if(std::isnan(a) ==1){ std::cout << mu << "  " << "\n";}

    return Eigen::Vector3d{x,y,z};
}

std::vector<double> calcStreamDirs(const std::vector<double>& mus, const std::vector<double>& phis, size_t nmu, size_t nphi) {
    std::vector<double> streams;
    streams.resize(nmu*nphi*3);
    for(size_t i = 0; i < nmu; ++i) {
        for(size_t j = 0; j < nphi; ++j) {
            Eigen::Vector3d dir = angleToVec(mus[i],phis[j]);
            for(size_t x = 0; x < 3; ++x) {
                size_t idx = indexRecompose(std::array<size_t,3>{i,j,x}, std::array<size_t,3>{nmu,nphi,3});
                streams[idx] = dir[x];
            }
        }
    }
    return streams;
}

std::vector<double> calcSubstreamDirs(const std::vector<double>& mus, const std::vector<double>& phis, const std::vector<double>& wmus, const std::vector<double>& wphis, size_t nmu, size_t nphi, size_t nsub) {
    double wmu_s = -1;
    double wmu_e;
    std::vector<double> substreams = {};
    substreams.resize(nmu*nphi*nsub*nsub*3);
    for(size_t i = 0; i < nmu; ++i) {
        wmu_e = wmu_s + wmus[i];
        for(size_t j = 0; j < nphi; ++j) {
            double delmu = (wmu_e-wmu_s)/nsub;
            double delphi = wphis[j]/nsub;
            double wphi_s = phis[j] - wphis[j]/2;
            double wphi_e = phis[j] + wphis[j]/2;
            for(size_t is = 0; is<nsub; ++is){
                double m = wmu_s + (is+0.5) * delmu;
                for(size_t js = 0; js<nsub; ++js){
                    double p = wphi_s + (js+0.5) * delphi;
                    Eigen::Vector3d dir = angleToVec(m,p);
                    for(size_t x = 0; x < 3; ++x) {
                        size_t idx = indexRecompose(std::array<size_t,5>{i,j,is,js,x}, std::array<size_t,5>{nmu,nphi,nsub,nsub,3});
                        substreams[idx] = dir[x];
                    }
                }
            }
        }
        wmu_s = wmu_e;
    }
    return substreams;
}


std::array<double,3> groundReflection_lambert(auto ray, size_t idx, double albedo, size_t nx, size_t ny, size_t nlyr, size_t nmu, size_t nphi, std::vector<double> mu, std::vector<double> wmu, std::vector<double> wphi, const std::vector<double>& Edir, const std::vector<double>& rad) {    
    double E_dir_ref = Edir[idx] * albedo / M_PI;
    auto [x,y,z] = indexDecompose(idx, std::array<size_t,3>{nx,ny,nlyr+1});  
    double E_diff = 0;
    for(size_t i = 0; i<nmu; ++i) {
        for(size_t j = 0; j<nphi; ++j) {
            size_t index = indexRecompose(std::array<size_t,5>{x,y,z,i,j},std::array<size_t,5>{nx,ny,nlyr+1,nmu,nphi});
            if(mu[i] < 0) {
                E_diff += rad[index] * fabs(mu[i])  * wmu[i] * wphi[j];	
            }
        }
    }
    double E_diff_ref = E_diff * albedo / M_PI;
    double E_ref = E_dir_ref + E_diff_ref;
    return std::array<double,3>{E_ref,E_dir_ref,E_diff_ref};
}

double phase_HG(double g, double mu) {
    double f = 1+g*g-2*g*mu;
    return 1./(4*M_PI) * (1-g*g)/std::sqrt(f*f*f);
}


double calc_pHG(double wmu_s, double wmu_e, double phi, double wphi, const Ray& ray, double g, size_t n, const std::vector<double>& substreams, size_t mu_i, size_t phi_j, size_t nmu, size_t nphi){
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
            //Eigen::Vector3d dir = angleToVec(m,p);
            Eigen::Vector3d dir;
            for(size_t x = 0; x<3; ++x) {
                size_t li = indexRecompose(std::array<size_t,5>{mu_i,phi_j,i,j,x}, std::array<size_t,5>{nmu,nphi,n,n,3});
                dir[x] = substreams[li];
            }
            double sa = (-ray.d).dot(dir);
            pf += phase_HG(g, sa);
            weightsum += 1;
            if(pf < 0){std::cout << " negative \n";}
        }
    }
    return pf/(n*n);
}

std::array<double,3> calc_Ldiff(const Ray& ray, double tfar, double tnear, size_t idx, double g1, size_t nx, size_t ny, size_t nlyr, size_t nmu, size_t nphi,
        const std::vector<double>& mu, const std::vector<double>& phi, const std::vector<double>& wmu, const std::vector<double>& wphi, const std::vector<double>& rad, 
        const std::vector<double>& streams, size_t nsub, const std::vector<double>& substreams) {
    double Lup = 0;
    double Ldown = 0;
    auto [x,y,z] = indexDecompose<3>(idx, {nx,ny,nlyr});
    double wmu_s = -1;
    double wmu_e,wphi_e;
    double pf;

    for(size_t i = 0; i < nmu; ++i) {
        wmu_e = wmu_s + wmu[i];
        for(size_t j = 0; j < nphi; ++j) {
            //Eigen::Vector3d lvec = angleToVec(mu[i],phi[j]);
            Eigen::Vector3d lvec;
            for(size_t x = 0; x<3; ++x) {
                size_t li = indexRecompose(std::array<size_t,3>{i,j,x}, std::array<size_t,3>{nmu,nphi,3});
                lvec[x] = streams[li];
            }

            auto muscatter = (-ray.d).dot(lvec); 
            double weight = wmu[i] * wphi[j];
            size_t index_t = indexRecompose(std::array<size_t,5>{x,y,z+1,i,j},std::array<size_t,5>{nx,ny,nlyr+1,nmu,nphi});
            size_t index_b = indexRecompose(std::array<size_t,5>{x,y,z,i,j},std::array<size_t,5>{nx,ny,nlyr+1,nmu,nphi});
            if(0==0){
                pf = weight * phase_HG(g1, muscatter);
            }
            else{
                pf = weight * calc_pHG(wmu_s, wmu_e, phi[j], wphi[j], ray, g1, nsub, substreams,i,j,nmu,nphi);
            }
            if(mu[i] < 0){
                Ldown +=  rad[index_t] * pf;
            }
            if(mu[i] > 0){
                Lup += rad[index_b] * pf;
            //    std::cout << "rad = " << rad[index_b] << " muscatter = " << muscatter << " ray = "  << -ray.d.transpose() << " lvec = " << lvec.transpose() << " mu = " << mu[i] << "phi = "<< phi[j] << " pf ="  << pf << "\n";
            }

        }
        wmu_s = wmu_e;
    }
    
    return std::array<double,3>{Lup + Ldown,Lup,Ldown}; 
}


