#include <iostream>
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include "speGrid.h"
/*
 * Created by Andrew S Maxwell 12/12/2023
 *
 * NSDI code using a sin square pulse envelope
 */

int main(int argc, char **argv) {
	speGridParameters gridParam;
	gridParam.Nx1 = gridParam.Nz1 = gridParam.Nx2 = gridParam.Nz2 = 10;
	gridParam.pfStart={-3, -3};
	gridParam.pfEnd={3, 3};
	gridParam.E01 = 0.79;
	gridParam.E02 = 1.51;

	//laser parameters
	double omega = 0.057, rtUp = std::sqrt(1.2), phi = 0.;
	int N = 1;
	fieldTypes fieldType = fieldTypes::monochromatic;

	laserField LF = std::make_shared<fields::monochromaticField>(omega, rtUp, phi, N, fieldType);
	speGrid saddlePointGrid(gridParam, LF);
//	saddlePointGrid.populateGrid();
//	saddlePointGrid.solveAllRandom();
	int result = Catch::Session().run( argc, argv );
//	int result = 0;


	return result;
}

TEST_CASE( "1: All test cases reside in other .cpp files (empty)", "[multi-file:1]" ) {
}
