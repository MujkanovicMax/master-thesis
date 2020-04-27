#include "catch.hpp"
#include "functions.hpp"
#include <netcdf>
#include <iostream>
#include <array>
#include <cmath>
#include <numeric>
#include <Eigen/Dense>
#include <stdexcept>

TEST_CASE("out of bounds index composition"){

    std::array<size_t,3> shape = {4,3,2};
    CHECK_THROWS_AS(indexRecompose(std::array<size_t,3>{0,0,2},shape), std::out_of_range);
    CHECK_THROWS_AS(indexRecompose(std::array<size_t,3>{4,0,0},shape), std::out_of_range);



}

TEST_CASE("check index"){
    std::array<size_t,3> shape = {2,3,4};
    REQUIRE(indexRecompose(indexDecompose<3>(17, shape),shape) == 17);
}

TEST_CASE("check values/index"){
    std::vector<double> radiances;
    std::vector<double> mus;
    std::vector<double> phis;
    std::vector<double> wmus;
    std::vector<double> wphis;
    size_t nx,ny,nz,nwvl,nmu,nphi;

    {
        using namespace netCDF;

        NcFile file("radiances_t.nc", NcFile::FileMode::read);
        nx = file.getDim("x").getSize();
        ny = file.getDim("y").getSize();
        nz = file.getDim("z").getSize();
        nphi = file.getDim("phi").getSize();
        nmu = file.getDim("mu").getSize(); 
        nwvl  = file.getDim("wvl").getSize();

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
        
    }

    size_t x = 31;
    size_t y = 0;
    size_t z = 22;
    size_t w = 0;
    size_t m = 5;
    size_t p = 6;

    auto idx = indexRecompose(std::array<size_t,6>{x,y,z,w,m,p},std::array<size_t,6>{ny,nx,nz,nwvl,nmu,nphi}); 
    auto [a,b,c,d,e,f] = indexDecompose<6>(idx,std::array<size_t,6>{ny,nx,nz,nwvl,nmu,nphi});
    idx = indexRecompose(std::array<size_t,6>{b,a,c,d,e,f},std::array<size_t,6>{nx,ny,nz,nwvl,nmu,nphi});

    REQUIRE(radiances[idx] == 159.804);
}





