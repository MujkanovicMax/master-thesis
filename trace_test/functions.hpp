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
#include <string>
#include <fstream>
#include <algorithm>


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
            double delphi = wphis[j]/M_PI*180/nsub;
            double wphi_s = phis[j] - wphis[j]/M_PI*180/2; //def wphi_s = phis[j] - wphis[j]/2 wahrscheinlich falsch weil phi in ° und wphi in rad
            double wphi_e = phis[j] + wphis[j]/M_PI*180/2; //def wphi_e = phis[j] + wphis[j]/2

            //if(wphi_s < 0){ wphi_s = 0;}
            //if(wphi_e > 360){ wphi_e = 360;}
            //std::cout << wphi_s << "    " << phis[j] << "   " << wphis[j] << "\n";
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
    //std::cout << Edir[idx] << "     " << x << "     " << y << "     " << z << "\n"; 
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

double integrated_HG(double g, double mu_s, double mu_e){

    if(g==0){
        return 1/(4*M_PI);//*(mu_e-mu_s);
    }
    else{
        return (1-g*g)/(4*M_PI*g*std::sqrt(1+g*g-2*g*mu_e)) - (1-g*g)/(4*M_PI*g*std::sqrt(1+g*g-2*g*mu_s));
    }
}

double calc_IHG( double wmu_s, double wmu_e, double phi, double wphi, const Ray& ray, double g){

    return 0;

}

double calc_pHG(double wmu_s, double wmu_e, double phi, double wphi, const Ray& ray, double g, size_t n, const std::vector<double>& substreams, size_t mu_i, size_t phi_j, size_t nmu, size_t nphi){
    //double delmu = (wmu_e-wmu_s)/n;
    //double delphi = wphi/n;
    double pf = 0;
    //double weightsum = 0;
    //double wphi_s = phi - wphi/2;
    //double wphi_e = phi + wphi/2;

    ///test///
    // Eigen::Vector3d dir_s, dir_e;
    // for(size_t x = 0; x<3; ++x) {
    //         size_t l1 = indexRecompose(std::array<size_t,5>{mu_i,phi_j,0,0,x}, std::array<size_t,5>{nmu,nphi,n,n,3});
    //         size_t l2 = indexRecompose(std::array<size_t,5>{mu_i,phi_j,n-1,n-1,x}, std::array<size_t,5>{nmu,nphi,n,n,3});
    //         dir_s[x] = substreams[l1];
    //         dir_e[x] = substreams[l2];
    //     }
    // double mu_s = (-ray.d).dot(dir_s);
    // double mu_e = (-ray.d).dot(dir_e);
    // 
    // double pf_alt = integrated_HG(g,mu_s,mu_e); 
    // //return pf;
    //////////
    //double sa_min=10;
    double sa;
    //double sa_max=-10;
    for(size_t i = 0; i<n; ++i){
        //double m = wmu_s + (i+0.5) * delmu;
        for(size_t j = 0; j<n; ++j){
            //double p = wphi_s + (j+0.5) * delphi;
            //Eigen::Vector3d dir = angleToVec(m,p);
            Eigen::Vector3d dir;
            for(size_t x = 0; x<3; ++x) {
                size_t li = indexRecompose(std::array<size_t,5>{mu_i,phi_j,i,j,x}, std::array<size_t,5>{nmu,nphi,n,n,3});
                dir[x] = substreams[li];
            }
            sa = (-ray.d).dot(dir);
            //sa = -ray.d[0] * dir[0] - ray.d[1] * dir[1] - ray.d[2] * dir[2];
            //if(sa_min > sa){sa_min=sa;}
            //if(sa_max < sa){sa_max=sa;}
            pf += phase_HG(g, sa);
            //weightsum += 1;
            if(pf < 0){std::cout << " negative \n";}
        }
    }
    //double pf_alt = integrated_HG(g,sa_min,sa_max);
    //std::cout << "pf = " << pf/(n*n) << " pf_alt = " << pf_alt << " sa_s = " << sa_min << " sa_e = " << sa_max  << "\n";
    //return pf_alt;
    return pf/(n*n);
}

//std::array<double, 3> calc_Ldiff_tenstream(const Ray& ray, double dx, double dy, std::vector<double>& zlev,double tfar, double tnear, 
//        size_t idx, double kext, double dtau, double g1, size_t nx, size_t ny, size_t nlyr, size_t nmu, size_t nphi,
//        const std::vector<double>& mu, const std::vector<double>& phi, const std::vector<double>& wmu, const std::vector<double>& wphi, 
//        const std::vector<double>& rad, const std::vector<double>& streams, size_t nsub, const std::vector<double>& substreams) {
//
//    //    Tenstream Streams: For Box: top 2, left 4, back 4 
//    //    0     Eup
//    //    1     Edown
//    //    2     E x bottom left out
//    //    3     E x bottom left in
//    //    4     E x top left out
//    //    5     E x top left in
//    //    6     E y bottom back out
//    //    7     E y bottom back in
//    //    8     E y top back out
//    //    9    E y top back in
//    //    for all inward streams information from neighbouring boxes is needed !!
//    
//    Eigen::MatrixX3 dist(3,nstreams); 
//    Eigen::VectorXd w(nstreams);
//    Eigen::VectorXd w2(nstreams);
//    Eigen::MatrixXd m(3,nstreams);
//
//    Eigen::Vector3d P = ray(tnear + dist);
//    Eigen::Vector3d top((x+0.5)*dx,(y+0.5)*dy,zlev[z+1]);
//    Eigen::Vector3d bottom((x+0.5)*dx,(y+0.5)*dy,zlev[z]);
//    Eigen::Vector3d left_bottom(x*dx,(y+0.5)*dy,(zlev[z+1] - zlev[z])*0.25 + zlev[z]);
//    Eigen::Vector3d left_top(x*dx,(y+0.5)*dy,(zlev[z+1] - zlev[z])*0.75 + zlev[z]);
//    Eigen::Vector3d right_bottom((x+1)*dx,(y+0.5)*dy,(zlev[z+1] - zlev[z])*0.25 + zlev[z]);
//    Eigen::Vector3d right_top((x+1)*dx,(y+0.5)*dy,(zlev[z+1] - zlev[z])*0.75 + zlev[z]);
//    Eigen::Vector3d back_bottom((x+0.5)*dx,(y+1)*dy,(zlev[z+1] - zlev[z])*0.25 + zlev[z]);
//    Eigen::Vector3d back_top((x+0.5)*dx,(y+1)*dy,(zlev[z+1] - zlev[z])*0.75 + zlev[z]);
//    Eigen::Vector3d front_bottom((x+0.5)*dx,(y+0)*dy,(zlev[z+1] - zlev[z])*0.25 + zlev[z]);
//    Eigen::Vector3d front_top((x+0.5)*dx,(y+0)*dy,(zlev[z+1] - zlev[z])*0.75 + zlev[z]);
//    
//    m.col(0) = top;
//    m.col(1) = left_top;
//    m.col(2) = right_top;
//    m.col(3) = back_top;
//    m.col(4) = front_top;
//    m.col(5) = bottom;
//    m.col(6) = left_bottom;
//    m.col(7) = right_bottom;
//    m.col(8) = back_bottom;
//    m.col(9) = front_bottom;
//    
//    double total_dist = 0;
//    for(size_t j=0; j<nstreams; ++j){
//        
//        dist.col(j) = P-m.col(j);
//        for(size_t i; i<3; ++i){
//            
//            double x1 = dist.col(j)[0];
//            double x2 = dist.col(j)[1];
//            double x3 = dist.col(j)[2];
//            
//            w[j] = std::sqrt(x1*x1 + x2*x2 + x3*x3);
//            total_dist += w[j]; 
//
//        }
//
//    }
//    
//    for(size_t j=0; j<nstreams; ++j){
//
//        w[j] *= w[j]*w[j];
//        w[j] = 1-w[j]/total_dist;
//
//    }
//    
//    Eigen::VectorXd iup(nstreams);
//    Eigen::VectorXd idown(nstreams);
//
//    //downward streams
//    idown[0] = indexRecompose(std::array<size_t,3>{y,x,z,1},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//    idown[1] = indexRecompose(std::array<size_t,3>{y,x,z,5},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//    idown[2] = indexRecompose(std::array<size_t,3>{y,x+1,z,4},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//    idown[3] = indexRecompose(std::array<size_t,3>{y,x,z,9},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//    idown[4] = indexRecompose(std::array<size_t,3>{y+1,x,z,8},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//    idown[5] = indexRecompose(std::array<size_t,3>{y,x,z-1,1},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//    idown[6] = indexRecompose(std::array<size_t,3>{y,x,z,2},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//    idown[7] = indexRecompose(std::array<size_t,3>{y,x+1,z,3},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//    idown[8] = indexRecompose(std::array<size_t,3>{y,x,z,6},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//    idown[9] = indexRecompose(std::array<size_t,3>{y+1,x,z,7},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//    
//    //upward streams
//    iup[0] = indexRecompose(std::array<size_t,3>{y,x,z,0},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//    iup[1] = indexRecompose(std::array<size_t,3>{y,x,z,4},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//    iup[2] = indexRecompose(std::array<size_t,3>{y,x+1,z,5},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//    iup[3] = indexRecompose(std::array<size_t,3>{y,x,z,8},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//    iup[4] = indexRecompose(std::array<size_t,3>{y+1,x,z,9},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//    iup[5] = indexRecompose(std::array<size_t,3>{y,x,z-1,0},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//    iup[6] = indexRecompose(std::array<size_t,3>{y,x,z,3},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//    iup[7] = indexRecompose(std::array<size_t,3>{y,x+1,z,2},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//    iup[8] = indexRecompose(std::array<size_t,3>{y,x,z,7},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//    iup[9] = indexRecompose(std::array<size_t,3>{y+1,x,z,6},std::array<size_t,5>{ny,nx,nlyr+1,nstream});
//
//
//
//
//    for(size_t i = 0; i<nstreams; ++i){
//        
//        Ldown += Et[idown[i]]*w[i]/cos[i];
//        Lup += Et[iup[i]]*w[i]/cos[i];
//
//    }
//
//
//
//
//
//
//}


std::array<double,3> calc_Ldiff(const Ray& ray, double dx, double dy, const std::vector<double>& zlev,double tfar, double tnear, size_t idx, double kext, 
        double dtau, double g1, size_t nx, size_t ny, size_t nlyr, size_t nmu, size_t nphi, const std::vector<double>& mu, const std::vector<double>& phi, 
        const std::vector<double>& wmu, const std::vector<double>& wphi, const std::vector<double>& rad, 
        const std::vector<double>& streams, size_t nsub, const std::vector<double>& substreams) {

    double Lup = 0;
    double Ldown = 0;
    auto [x,y,z] = indexDecompose<3>(idx, {nx,ny,nlyr});
    double wmu_s = -1;
    double wmu_e,wphi_e;
    double pf;
    double tau_thres = 1;
    double alpha, beta;

    if(kext == 0){
        alpha = 0.5;
        beta = 0.5;
    }

    else{
        double dist = tau_thres/kext;
        Eigen::Vector3d P = ray(tnear + dist);
        double l1 = (x+0.5)*dx;
        double l2 = (y+0.5)*dy;
        double l3 = zlev[z];
        double u3 = zlev[z+1];
        double pu = (P[0]-l1)*(P[0]-l1) + (P[1]-l2)*(P[1]-l2) + (P[2]-u3)*(P[2]-u3);
        double pl = (P[0]-l1)*(P[0]-l1) + (P[1]-l2)*(P[1]-l2) + (P[2]-l3)*(P[2]-l3);
        pu *= pu*pu;
        pl *= pl*pl;
        alpha =  1-pu/(pl+pu);
        beta = 1 - alpha;
    }

    if(ray.d[2] > 0){

        alpha = 1 - alpha;
        beta = 1 - alpha;

    }

    for(size_t i = 0; i < nmu; ++i) {
        wmu_e = wmu_s + wmu[i];
        for(size_t j = 0; j < nphi; ++j) {
            Eigen::Vector3d lvec;
            double weight = wmu[i] * wphi[j];
            size_t index_t = indexRecompose(std::array<size_t,5>{x,y,z+1,i,j},std::array<size_t,5>{nx,ny,nlyr+1,nmu,nphi});
            size_t index_b = indexRecompose(std::array<size_t,5>{x,y,z,i,j},std::array<size_t,5>{nx,ny,nlyr+1,nmu,nphi});

            if(nsub == 0){

                if(g1==0){  pf = weight * 1./(4*M_PI); }

                else{

                    for(size_t x = 0; x<3; ++x) {
                        size_t li = indexRecompose(std::array<size_t,3>{i,j,x}, std::array<size_t,3>{nmu,nphi,3});
                        lvec[x] = streams[li];
                    }

                    double muscatter = (-ray.d).dot(lvec);
                    pf = weight * phase_HG(g1, muscatter);

                }

            }

            else{

                if(g1 == 0){ pf = weight * 1./(4*M_PI); }

                else{

                    for(size_t x = 0; x<3; ++x) {
                        size_t li = indexRecompose(std::array<size_t,3>{i,j,x}, std::array<size_t,3>{nmu,nphi,3});
                        lvec[x] = streams[li];
                    }
                    double muscatter = (-ray.d).dot(lvec);
                    pf = weight * calc_pHG(wmu_s, wmu_e, phi[j], wphi[j], ray, g1, nsub, substreams,i,j,nmu,nphi);

                }

            }

            if(mu[i] < 0){
                Ldown += alpha * rad[index_t] * pf  + beta * rad[index_b] * pf;
            }
            else{
                Lup += alpha * rad[index_t] * pf + beta * rad[index_b] * pf;
            }


        }
        wmu_s = wmu_e;
    }
    return std::array<double,3>{Lup + Ldown,Lup,Ldown}; 
}
void calc_image(auto grid, auto cam, double Nxpixel, double Nypixel, double dx, double dy, const std::vector<double>& zlev, const std::vector<double>& kext, 
        const std::vector<double>& g1, const std::vector<double>& w0, double albedo, double muEdir, size_t nx, size_t ny, size_t nlyr, 
        size_t nmu, size_t nphi, const std::vector<double>& mus, const std::vector<double>& phis, const std::vector<double>& wmus, 
        const std::vector<double>& wphis, const Eigen::Vector3d& sza_dir, const std::vector<double>& Edir, const std::vector<double>& radiances, 
        const std::vector<double>& streams, size_t nsub, 
        const std::vector<double>& substreams, std::vector<double>& image, std::vector<double>& opthick_image, std::vector<double>& Ldiff_i, 
        std::vector<double>& Lup_i, std::vector<double>& Ldown_i, std::vector<double>& Ldir_i, std::vector<double>& groundbox, 
        std::vector<double>& gRdir_i, std::vector<double>& gRdiff_i)
{

    //Loop variables
    size_t a,b,c;
    //Main loop
    std::cout << "Starting Raytracer..." << "\n";

    for (size_t i = 0; i < Nxpixel; ++i){
        double xpx = (i + 0.5) / Nxpixel;

        for(size_t j = 0; j < Nypixel; ++j) {

            double ypx = (j + 0.5) / Nypixel;
            auto ray = cam.compute_ray(Eigen::Vector2d{xpx, ypx});

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
                    size_t optprop_index = indexRecompose(std::array<size_t,3>{z,x,y},std::array<size_t,3>{nlyr,nx,ny});
                    size_t rad_index = indexRecompose(std::array<size_t,3>{x,y,z+1},std::array<size_t,3>{nx,ny,nlyr+1});
                    groundidx = indexRecompose(std::array<size_t,3>{x,y,z},std::array<size_t,3>{nx,ny,nlyr+1});

                    //main pixel calculation ( radiance summation, etc.)
                    double L = Edir[rad_index]/fabs(muEdir);
                    double transmission = exp(-optical_thickness);
                    double dtau = (pvol->tfar - pvol->tnear) * kext[optprop_index];
                    optical_thickness += dtau;
                    double phase_function = phase_HG(g1[optprop_index], (-ray.d).dot(sza_dir.normalized()));
                    auto [Lup_Plus_Ldown, Lup, Ldown] = calc_Ldiff(ray, dx, dy, zlev, pvol->tfar, pvol->tnear, pvol->idx, kext[optprop_index], 
                            dtau, g1[optprop_index], nx, ny, nlyr, nmu, nphi, mus, phis, wmus, wphis, radiances, streams, nsub, substreams);
                    double scatter_prob = - expm1(-dtau*w0[optprop_index]);

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
            double gR, gRdir, gRdiff = 0;
            std::array<double,3> ground{gR, gRdir, gRdiff};
            auto [x,y,z] = indexDecompose(groundidx, std::array<size_t,3>{nx,ny,nlyr+1}); 
            if( z == 0){
                ground = groundReflection_lambert(ray,groundidx,albedo, nx, ny, nlyr, nmu, nphi, mus, wmus, wphis, Edir, radiances); 
            }
            double ground_E = ground[0] *  exp(-optical_thickness);            

            //writing pixelvalues to image
            image[j + i * Nypixel]      = radiance + ground_E;
            opthick_image[j+i*Nypixel]  = optical_thickness;
            Ldiff_i[j+i*Nypixel]        = Ldiffs;
            Lup_i[j+i*Nypixel]          = Lups;
            Ldown_i[j+i*Nypixel]        = Ldowns;
            Ldir_i[j+i*Nypixel]         = Ldirs;
            groundbox[j+i*Nypixel]      = indexDecompose<3>(groundidx,std::array<size_t,3>{nx,ny,nlyr})[0];
            gRdir_i[j+i*Nypixel]        = ground[1] * exp(-optical_thickness);
            gRdiff_i[j+i*Nypixel]       = ground[2] * exp(-optical_thickness);


            std::cout << "Pixel " << i << " ," << j << " done." << "\n";
        }
    }

    std::cout << "Stopping Raytracer...\n";
}


void read_radiances( std::string fpath, std::vector<double>& radiances, std::vector<double>& mus, std::vector<double>& phis, std::vector<double>& wmus, std::vector<double>& wphis, size_t& nmu, size_t& nphi )
{

    //try{

    using namespace netCDF;

    NcFile file(fpath, NcFile::FileMode::read);

    size_t nx = file.getDim("x").getSize();
    size_t ny = file.getDim("y").getSize();
    size_t nz = file.getDim("z").getSize();
    nphi = file.getDim("phi").getSize();
    nmu = file.getDim("mu").getSize();
    size_t nwvl  = file.getDim("wvl").getSize();
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

    //}

    //catch(const netCDF::exceptions::NcBadId& e){

    //    using namespace netCDF;

    //    NcFile file(fpath, NcFile::FileMode::read);
    //    
    //    size_t nx = file.getDim("x").getSize();
    //    size_t ny = file.getDim("y").getSize();
    //    size_t nz = file.getDim("z").getSize();
    //    size_t nwvl  = file.getDim("wvl").getSize();
    //    radiances.resize(nx*ny*nz*nwvl);

    //    file.getVar("radiance").getVar(radiances.data());
    //
    //}




}

void read_opprop( std::string fpath, std::vector<double>& kext, std::vector<double>& zlev, std::vector<double>& w0, std::vector<double>& g1, size_t& nlev, size_t& nlyr, size_t& nx, size_t& ny )
{

    using namespace netCDF;

    NcFile file(fpath, NcFile::FileMode::read);
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
void read_flx( std::string fpath, std::vector<double>& Edir, std::vector<double>& Edown, Eigen::Vector3d& sza_dir, double& muEdir )
{

    using namespace netCDF;

    NcFile file(fpath, NcFile::FileMode::read);
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

void print_parameters(std::string radpath, std::string flxpath, std::string oppath, std::string fname, double dx, double dy, double albedo, double xloc,
        double yloc, double zloc, double Nxpixel, double Nypixel, double fov, double fov_phi1, double fov_phi2, double fov_theta1, 
        double fov_theta2, double rays, double nsub)
{
    std::ofstream file;
    file.open("parameter_logs/" + fname + "_input_parameters.txt");

    file << "radpath: " << radpath  << "\n";
    file << "flxpath: " << flxpath  << "\n";
    file << "oppath: " << oppath  << "\n";
    file << "dx = " << dx  << "\n";
    file << "dy = " << dy  << "\n";
    file << "albedo = " << albedo  << "\n";
    file << "xloc = " << xloc << "\n";
    file << "yloc = " << yloc << "\n";
    file << "zloc = " << zloc << "\n";
    file << "Nxpixel = " << Nxpixel  << "\n";
    file << "Nypixel = " << Nypixel  << "\n";
    file << "fov = " << fov  << "\n";
    file << "fov_phi1 = " << fov_phi1  << "\n";
    file << "fov_phi2 = " << fov_phi2  << "\n";
    file << "fov_theta1 = " << fov_theta1  << "\n";
    file << "fov_theta2 = " << fov_theta2  << "\n";
    file << "rays = " << rays  << "\n";
    file << "nsub = " << nsub  << "\n";
    file.close();


}

void configparser(std::string fname, std::string& radfn, std::string& flxfn, std::string& opfn, std::string& outfn, double& dx, double& dy, double& albedo, \
        int& xpixel, int& ypixel, double& fov_phi1, double& fov_phi2, double& fov_theta1, double& fov_theta2, double& xloc, double& yloc, \
        double& zloc, int& nsub){

    ///////default values//////
    radfn = "../radiances/rad_testcases/checkerboard_full/rad_checkerboard.nc";
    flxfn = "../radiances/rad_testcases/checkerboard_full/flx_checkerboard.nc";
    opfn = "../radiances/rad_testcases/checkerboard_full/op_checkerboard.nc";
    outfn = "output2d_default_checkerboard.nc";
    dx = 1;
    dy = 1;
    albedo = 0.2;
    xpixel = 100; 
    ypixel = 100; 
    fov_phi1 = -225;
    fov_phi2 = -135;
    fov_theta1 = 45; 
    fov_theta2 = 135;
    xloc = 0; //in km
    yloc = 0;
    zloc = 0;
    nsub = 1;
    //////////////////////////

    std::ifstream file(fname);
    if(file.is_open())
    {
        std::string line;
        while(getline(file,line))
        {

            line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
            if(line[0] == '#' || line.empty()){continue;}

            auto delimPos = line.find("=");
            auto name = line.substr(0, delimPos);
            auto value = line.substr(delimPos + 1);

            if(name == "rad_filename"){ radfn = value; continue;}
            if(name == "flx_filename"){ flxfn = value; continue;}
            if(name == "op_filename"){ opfn = value; continue;}
            if(name == "output_filename"){ outfn = value; continue;}
            if(name == "dx"){ dx = std::stod(value); continue;}
            if(name == "dy"){ dy = std::stod(value); continue;}
            if(name == "albedo"){ albedo = std::stod(value); continue;}
            if(name == "xpixel"){ xpixel = std::stoi(value); continue;}
            if(name == "ypixel"){ ypixel = std::stoi(value); continue;}
            if(name == "fov_phi1"){ fov_phi1 = std::stod(value); continue;}
            if(name == "fov_phi2"){ fov_phi2 = std::stod(value); continue;}
            if(name == "fov_theta1"){ fov_theta1 = std::stod(value); continue;}
            if(name == "fov_theta2"){ fov_theta2 = std::stod(value); continue;}
            if(name == "xloc"){ xloc = std::stod(value); continue;}
            if(name == "yloc"){ yloc = std::stod(value); continue;}
            if(name == "zloc"){ zloc = std::stod(value); continue;}
            if(name == "nsub"){ nsub = std::stoi(value); continue;}
            else{ std::cerr << "Name " << name << " not recognized\n";}
        }


    }

    else{
        std::cerr << "Couldn't open/find config file.\n";
    }

    file.close();

    std::ifstream inFile(fname);
    std::ofstream outFile("parameter_logs/config_" + outfn + ".txt");

    outFile << inFile.rdbuf();

}


