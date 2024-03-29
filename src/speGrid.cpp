/*
 * speGrid.cpp
 *
 *  Created on: 12 Dec 2023
 *      Author: andy
 */
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ranges>
#include <tuple>

#include "speGrid.h"




const double pi = 3.1415926535897932384626433832795028841971693993751;

spePoint::spePoint(vec4 pf_) : pf(pf_){
	ti=tr=k={};
	pointSolvedQ=false;
}

speGrid::speGrid(const speGridParameters& gridParam, const laserField LF) : Nx1(gridParam.Nx1), Nz1(gridParam.Nz1)
, Nx2(gridParam.Nx2), Nz2(gridParam.Nz2), pfStart(gridParam.pfStart), pfEnd(gridParam.pfEnd)
, saddlePointEquations(speEqns(gridParam.E01, gridParam.E02, LF)){

	//Set grid step size
	ddpx1 = (pfEnd[0] - pfStart[0])/(Nx1-1);
	ddpz1 = (pfEnd[1] - pfStart[1])/(Nz1-1);
	ddpx2 = (pfEnd[0] - pfStart[0])/(Nx2-1);
	ddpz2 = (pfEnd[1] - pfStart[1])/(Nz2-1);
}

speGrid::~speGrid() {};

void speGrid::populateGrid(){
	/*
	 * Function to populate the grid (grid) with momentum points and set other variables to zero
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

void speGrid::solvePointRandom(spePoint& point){
	/*
	 * This function uses the random solver from the saddle point equations to find all solutions
	 * for a specific point. Using a high number of random numbers will ensure all solutions are found
	 */
	size_t NumRand = 2000;

	std::vector<timePair> tsList = saddlePointEquations.RandSolve(point.pf, NumRand);
	for(auto ts : tsList){
		dcmplx ti = ts[0], tr = ts[1];
		dcmplx k = saddlePointEquations.k(ti, tr);
		point.ti.push_back(ti);
		point.tr.push_back(tr);
		point.k.push_back(k);
	}
	point.pointSolvedQ=true;
	return;
}
void speGrid::solvePointRandom(size_t ix1, size_t iz1, size_t ix2, size_t iz2){
	/*
	 * This function uses the random solver from the saddle point equations to find all solutions
	 * for a specific point. Using a high number of random numbers will ensure all solutions are found
	 */

	solvePointRandom(at(ix1, iz1, ix2, iz2));
	return;
}
void speGrid::solveAllRandom(){
	int i=0;
	for(auto& point: grid){
		++i;
		std::cout<<" solving point randomly: "<<i<<" out of "<<grid.size()<<"\n";
		solvePointRandom(point);
	}
}

void speGrid::solveWithAdjacent(spePoint& pointSolved, spePoint& pointToSolve){
	/*
	 * Function take solutions from one point and uses them as guess for another
	 * This should be used with adjacent solutions on a grid to easily propagate solutions across the grid
	 */
	if(pointSolved.pointSolvedQ==true && pointToSolve.pointSolvedQ == false){
		vec4 pf = pointToSolve.pf;
		std::vector<dcmplx> tiList = pointSolved.ti, trList = pointSolved.tr;
		for(auto i=0; i!=tiList.size(); ++i){
			dcmplx ti = tiList[i], tr = trList[i];
			int errFlag = saddlePointEquations.solveRoot(ti, tr, pf);
			if(errFlag==0){
				dcmplx k = saddlePointEquations.k(ti, tr);
				pointToSolve.ti.push_back(ti);
				pointToSolve.tr.push_back(tr);
				pointToSolve.k.push_back(k);
				pointToSolve.pointSolvedQ=true;
			}
		}
		if(!pointToSolve.pointSolvedQ){solvePointRandom(pointToSolve);}

	}else{
		std::cerr<<"Error either pointSolved is not solved or pointToSolve is already solved\n";
	}
}

void speGrid::propagateSolutionOverGrid(size_t ix1, size_t iz1, size_t ix2, size_t iz2){
	auto& point = at(ix1, iz1, ix2, iz2);
	if(point.pointSolvedQ==false){std::cerr<<"Error input point not solved\n"; return;}

	//solve a z1 row
	std::cout<<"solving row\n";
	for(int i=iz1+1; i!=Nz1; ++i){
//		std::cout<<"i = "<<i<<"\n";
		auto &pointSolved = at(ix1, i-1, ix2, iz2), &pointToSolve = at(ix1, i, ix2, iz2);
		solveWithAdjacent(pointSolved, pointToSolve);
	}
	for(int i=iz1; i!=0; --i){
//		std::cout<<"i = "<<i<<"\n";
		auto &pointSolved = at(ix1, i, ix2, iz2), &pointToSolve = at(ix1, i-1, ix2, iz2);
		solveWithAdjacent(pointSolved, pointToSolve);;
	}

	//solve 2d z1xz2 grid
	std::cout<<"solving 2d grid\n";
	for(int iiz2=iz2+1; iiz2!=Nz2; ++iiz2){
		for(int iiz1=0; iiz1!=Nz1; ++iiz1){
			auto &pointSolved = at(ix1, iiz1, ix2, iiz2-1), &pointToSolve = at(ix1, iiz1, ix2, iiz2);
			solveWithAdjacent(pointSolved, pointToSolve);
		}
	}
	for(int iiz2=iz2; iiz2!=0; --iiz2){
		for(int iiz1=0; iiz1!=Nz1; ++iiz1){
			auto &pointSolved = at(ix1, iiz1, ix2, iiz2), &pointToSolve = at(ix1, iiz1, ix2, iiz2-1);
			solveWithAdjacent(pointSolved, pointToSolve);
		}
	}

	//solve 3d x1xz1xz2 grid
	size_t NThreads=8;
	std::cout<<"solving 3d grid\n";
	#pragma omp parallel num_threads(NThreads)
	#pragma omp for ordered
	for(int iiz1=0; iiz1!=Nz1; ++iiz1){
		#pragma omp task firstprivate(iiz1)
		for(int iix1=ix1+1; iix1!=Nx1; ++iix1){
			for(int iiz2=0; iiz2!=Nz2; ++iiz2){
				auto &pointSolved = at(iix1-1, iiz1, ix2, iiz2), &pointToSolve = at(iix1, iiz1, ix2, iiz2);
				solveWithAdjacent(pointSolved, pointToSolve);
			}
		}
	}
	#pragma omp parallel num_threads(NThreads)
	#pragma omp for ordered
	for(int iiz1=0; iiz1!=Nz1; ++iiz1){
		#pragma omp task firstprivate(iiz1)
		for(int iix1=ix1; iix1!=0; --iix1){
			for(int iiz2=0; iiz2!=Nz2; ++iiz2){
				auto &pointSolved = at(iix1, iiz1, ix2, iiz2), &pointToSolve = at(iix1-1, iiz1, ix2, iiz2);
				solveWithAdjacent(pointSolved, pointToSolve);
			}
		}
	}
	//solve 4d x1xz1xx2xz2 grid
	std::cout<<"solving 4d grid\n";
	#pragma omp parallel num_threads(NThreads)
	#pragma omp for ordered
	for(int iiz1=0; iiz1!=Nz1; ++iiz1){
		#pragma omp task firstprivate(iiz1)
		for(int iix2=ix1+1; iix2!=Nx2; ++iix2){
			for(int iiz2=0; iiz2!=Nz2; ++iiz2){
				for(int iix1=0; iix1!=Nx1; ++iix1){
					auto &pointSolved = at(iix1, iiz1, iix2-1, iiz2), &pointToSolve = at(iix1, iiz1, iix2, iiz2);
					solveWithAdjacent(pointSolved, pointToSolve);
				}
			}
		}
	}
	#pragma omp parallel num_threads(NThreads)
	#pragma omp for ordered
	for(int iiz1=0; iiz1!=Nz1; ++iiz1){
	#pragma omp task firstprivate(iiz1)
		for(int iix2=ix1; iix2!=0; --iix2){
			for(int iiz2=0; iiz2!=Nz2; ++iiz2){
				for(int iix1=0; iix1!=Nx1; ++iix1){
					auto &pointSolved = at(iix1, iiz1, iix2, iiz2), &pointToSolve = at(iix1, iiz1, iix2-1, iiz2);
					solveWithAdjacent(pointSolved, pointToSolve);
				}
			}
		}
	}
	return;
}

void speGrid::printToFile(std::string outputFilename){
	std::ofstream outputFile(outputFilename);
	if(outputFile.fail()){
		std::cerr<<"Output file not opened\n";
	}
	for(const auto& point : grid){
		if(point.pointSolvedQ){
			vec4 pf = point.pf;
			std::vector<dcmplx> tiList = point.ti, trList = point.tr, kList = point.k;
			for(int i=0; i!=tiList.size(); ++i){
				outputFile//<<std::setprecision(std::numeric_limits<double>::digits10-2)
						<<pf[0]<<" "<<pf[1]<<" "<<pf[2]<<" "<<pf[3]<<" "<<tiList[i].real()<<" "
						<<tiList[i].imag()<<" "<<trList[i].real()<<" "<<trList[i].imag()<<" "
						<<kList[i].real()<<" "<<kList[i].imag()<<"\n";
			}
		}
	}
	return;
}

spePoint& speGrid::at(size_t ix1, size_t iz1, size_t ix2, size_t iz2){
	return grid[ix1*Nz1*Nx2*Nz2 + iz1*Nx2*Nz2 + ix2*Nz2 + iz2];
}

//get and set
const speEqns& speGrid::getSaddlePointEquations() const { return saddlePointEquations;}
const std::vector<spePoint>& speGrid::getGrid() const { return grid;}
const std::array<double,2>& speGrid::getPfStart() const { return pfStart;}
const std::array<double,2>& speGrid::getPfEnd() const { return pfEnd;}
const size_t speGrid::getNx1() const { return Nx1;}
const size_t speGrid::getNz1() const { return Nz1;}
const size_t speGrid::getNx2() const { return Nx2;}
const size_t speGrid::getNz2() const { return Nz2;}



