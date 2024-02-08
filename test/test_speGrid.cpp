/*
 * test_speGrid.cpp
 *
 *  Created on: 08 Feb 2024
 *      Author: andy
 */
#include <iostream>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>

#include "../src/speGrid.h"

const double pi = 3.1415926535897932384626433832795028841971693993751;
constexpr dcmplx I(0,1.);

/// --- Testing Suite --- ///
TEST_CASE( "Testing SPEGrid Methods" ,"[SPEGrid]") {
	speGridParameters gridParam;
	gridParam.Nx1=gridParam.Nz1=gridParam.Nx2=gridParam.Nz2=11;
	gridParam.pfStart={-3, -3};
	gridParam.pfEnd={3, 3};
	gridParam.E01 = 0.79;
	gridParam.E02 = 1.51;

	double omega=0.057, rtUp=std::sqrt(1.2);
	auto phi = GENERATE(0, 0.33*pi);
	size_t N = GENERATE(1, 2);
	fieldTypes fieldType = GENERATE(fieldTypes::monochromatic, fieldTypes::sin2);
	laserField LF;
	if(fieldType == monochromatic){
		LF = std::make_shared<fields::monochromaticField>(omega, rtUp, phi, N, fieldType);
	}else{
		LF = std::make_shared<fields::sin2>(omega, rtUp, phi, N, fieldType);
	}
	speGrid saddlePointGrid(gridParam, LF);

	//check populateGrid method
	saddlePointGrid.populateGrid();
	auto grid = saddlePointGrid.getGrid();
	SECTION( "Check the grid size" ) {
		//check grid not empty
		REQUIRE(!grid.empty());
		//validate grid size
		size_t gridLength = saddlePointGrid.getNx1()*saddlePointGrid.getNz1()*saddlePointGrid.getNx2()*saddlePointGrid.getNz2();
		REQUIRE(grid.size()==gridParam.Nx1*gridParam.Nz1*gridParam.Nx2*gridParam.Nz2);
	}
	SECTION( "Check values") {
		//check first element, vars that should be zero
		REQUIRE(grid[0].ti.empty());
		REQUIRE(grid[0].tr.empty());
		REQUIRE(grid[0].k.empty());
		REQUIRE(grid[0].pointSolvedQ==false);
		//check some random values (not too high to account for small grids)
		REQUIRE(grid[3].ti.empty());
		REQUIRE(grid[7].tr.empty());
		REQUIRE(grid[11].k.empty());
		REQUIRE(grid[13].pointSolvedQ==false);
		//check momentum values
		std::array<double,2> pfStart = saddlePointGrid.getPfStart();
		std::array<double,2> pfEnd = saddlePointGrid.getPfEnd();
		REQUIRE( ((grid[0].pf[0] == pfStart[0]) && (grid[0].pf[1] == pfStart[1]) ) );
		REQUIRE( ((grid.back().pf[0] == pfEnd[0]) && (grid.back().pf[1] == pfEnd[1]) ) );
	}
	//check solutions
	size_t midPoint = gridParam.Nx1/2;
	saddlePointGrid.solvePointRandom(midPoint, midPoint, midPoint, midPoint);
	spePoint& point = saddlePointGrid.at(midPoint, midPoint, midPoint, midPoint);
	std::cout<<" num sols = "<<point.ti.size()<<"\n";
	REQUIRE(!point.ti.empty());
	REQUIRE(!point.tr.empty());
	REQUIRE(!point.k.empty());
	REQUIRE(point.pointSolvedQ);
	dcmplx ti = point.ti[0], tr = point.tr[0];
	vec4 pf = point.pf;
	dcmplx SPE1val = saddlePointGrid.getSaddlePointEquations().SPE_ti(ti, tr, pf);
	dcmplx SPE2val = saddlePointGrid.getSaddlePointEquations().SPE_tr(ti, tr, pf);
	REQUIRE(std::abs(SPE1val+SPE2val) < 100 * 1e-12);

//	spePoint adjacentPoint = saddlePointGrid.at(4,3,4,4);
//	saddlePointGrid.solveWithAdjacent(point, adjacentPoint);
//	std::cout<<" num sols = "<<adjacentPoint.ti.size()<<"\n";
//	REQUIRE(!adjacentPoint.ti.empty());
//	REQUIRE(!adjacentPoint.tr.empty());
//	REQUIRE(!adjacentPoint.k.empty());
//	REQUIRE(adjacentPoint.pointSolvedQ);
//	ti = adjacentPoint.ti[0];
//	tr = adjacentPoint.tr[0];
//	pf = adjacentPoint.pf;
//	SPE1val = saddlePointGrid.getSaddlePointEquations().SPE_ti(ti, tr, pf);
//	SPE2val = saddlePointGrid.getSaddlePointEquations().SPE_tr(ti, tr, pf);
//	REQUIRE(std::abs(SPE1val+SPE2val) < 100 * 1e-12);

	saddlePointGrid.propagateSolutionOverGrid(midPoint, midPoint, midPoint, midPoint);


}
