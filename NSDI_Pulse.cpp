#include <iostream>
#include <catch2/catch_session.hpp>
#include "speGrid.h"
/*
 * Created by Andrew S Maxwell 12/12/2023
 *
 * NSDI code using a sin square pulse envelope
 */

int main(int argc, char **argv) {
	speGridParameters gridParam;
	gridParam.Nx1=gridParam.Nz1=gridParam.Nx2=gridParam.Nz2=10;
	gridParam.pfStart={-3, -3};
	gridParam.pfEnd={3, 3};

	//laser parameters
	double omega=0.057, rtUp=std::sqrt(1.2), phi=0.;
	int N = 1;
	fieldTypes fieldType1 = fieldTypes::monochromatic, fieldType2 = fieldTypes::sin2;

	fields::monochromaticField LF1(omega, rtUp, phi, N, fieldType1);
	fields::sin2 LF2(omega, rtUp, phi, N, fieldType2);
	speGrid saddlePointGrid(gridParam);
	int result = Catch::Session().run( argc, argv );

	std::cout<<result;

	return 0;
}
