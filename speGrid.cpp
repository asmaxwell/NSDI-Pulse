/*
 * speGrid.cpp
 *
 *  Created on: 12 Dec 2023
 *      Author: andy
 */

#include "speGrid.h"
#include <catch2/catch_test_macros.hpp>

spePoint::spePoint(vec4 pf_) : pf(pf_){
	ti=tr=k=0.;
}

speGrid::speGrid(const speGridParameters& gridParam) : Nx1(gridParam.Nx1), Nz1(gridParam.Nz1)
, Nx2(gridParam.Nx2), Nz2(gridParam.Nz2), pfStart(gridParam.pfStart), pfEnd(gridParam.pfEnd){
	//Set grid step size
	ddpx1 = (pfEnd[0] - pfStart[0])/(Nx1-1);
	ddpz1 = (pfEnd[1] - pfStart[1])/(Nz1-1);
	ddpx2 = (pfEnd[0] - pfStart[0])/(Nx2-1);
	ddpz2 = (pfEnd[1] - pfStart[1])/(Nz2-1);
}

speGrid::~speGrid() {
	// TODO Auto-generated destructor stub
}

void speGrid::populateGrid(){
	/*
	 * Function to populate the grid (grid) with momentum points
	 */
	for(size_t ix1=0; ix1!=Nx1; ++ix1){
		for(size_t iz1=0; iz1!=Nz1; ++iz1){
			for(size_t ix2=0; ix2!=Nx2; ++ix2){
				for(size_t iz2=0; iz2!=Nz2; ++iz2){
					double px1 = pfStart[0] + ix1*ddpx1, pz1 = pfStart[1] +iz1*ddpz1;
					double px2 = pfStart[0] + ix2*ddpx2, pz2 = pfStart[1] +iz2*ddpz2;
					vec4 pf ={px1, pz1, px2, pz2};
					spePoint point(pf);
					grid.push_back(point);
				}
			}
		}
	}
}

//get and set
const auto& speGrid::getGrid() const { return grid;}
const auto& speGrid::getPfStart() const { return pfStart;}
const auto& speGrid::getPfEnd() const { return pfEnd;}
const auto speGrid::getNx1() const { return Nx1;}
const auto speGrid::getNz1() const { return Nz1;}
const auto speGrid::getNx2() const { return Nx2;}
const auto speGrid::getNz2() const { return Nz2;}

/// --- Testing Suite --- ///
TEST_CASE( "Testing SPEGrid Methods" ,"[SPEGrid]") {
	speGridParameters gridParam;
	gridParam.Nx1=gridParam.Nz1=gridParam.Nx2=gridParam.Nz2=10;
	gridParam.pfStart={-3, -3};
	gridParam.pfEnd={3, 3};
	speGrid saddlePointGrid(gridParam);

	//check populateGrid method
	saddlePointGrid.populateGrid();
	auto grid = saddlePointGrid.getGrid();
	SECTION( "Check the grid size" ) {
		//check grid not empty
		REQUIRE(!grid.empty());
		//validate grid size
		size_t gridLength = saddlePointGrid.getNx1()*saddlePointGrid.getNz1()*saddlePointGrid.getNx2()*saddlePointGrid.getNz2();
		REQUIRE(grid.size()==gridLength);
	}
	SECTION( "Check values") {
		//check first element, vars that should be zero
		REQUIRE(grid[0].ti == std::complex<double>(0.,0.));
		REQUIRE(grid[0].tr == std::complex<double>(0.,0.));
		REQUIRE(grid[0].k == std::complex<double>(0.,0.));
		//check some random values (not too high to account for small grids)
		REQUIRE(grid[3].ti == std::complex<double>(0.,0.));
		REQUIRE(grid[7].tr == std::complex<double>(0.,0.));
		REQUIRE(grid[11].k == std::complex<double>(0.,0.));
		//check momentum values
		auto pfStart = saddlePointGrid.getPfStart();
		auto pfEnd = saddlePointGrid.getPfEnd();
		REQUIRE( ((grid[0].pf[0] == pfStart[0]) && (grid[0].pf[1] == pfStart[1]) ) );
		REQUIRE( ((grid.back().pf[0] == pfEnd[0]) && (grid.back().pf[1] == pfEnd[1]) ) );
	}
}


