#include <iostream>
#include <sstream>
#include <sys/stat.h>

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>

#include "actionData.h"
#include "speGrid.h"
/*
 * Created by Andrew S Maxwell 12/12/2023
 *
 * NSDI code using a sin square pulse envelope
 */

int main(int argc, char **argv) {

	const std::string PATH="Data/";
	//ensure Data folder exists not will not overwrite existing one
	struct stat st;
	if(stat(PATH.c_str(),&st)==-1){
		mkdir(PATH.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}

	speGridParameters gridParam;
	gridParam.Nz1 = gridParam.Nz2 = 201;
	gridParam.Nx1 = gridParam.Nx2 = 3;
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
	saddlePointGrid.populateGrid();
	size_t midPointx = gridParam.Nx1/2;
	size_t midPointz = gridParam.Nz1/2;
	saddlePointGrid.solvePointRandom(midPointx, midPointz, midPointx, midPointz);
	saddlePointGrid.propagateSolutionOverGrid(midPointx, midPointz, midPointx, midPointz);
	saddlePointGrid.printToFile(PATH+"SaddlePointSoluitons.dat");

	//compute action
	actionData actionGrid(saddlePointGrid);
	actionGrid.calculateAllAction();
	actionGrid.printToFile(PATH+"ActionData.dat");



	return 0;
}

