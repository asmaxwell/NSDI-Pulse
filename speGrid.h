/*
 * speGrid.h
 *
 *  Created on: 12 Dec 2023
 *      Author: andy
 *
 *      class to store and compute saddle point equations on a 6D momentum grid
 */

#ifndef SPEGRID_H_
#define SPEGRID_H_

#include <vector>

#include "speEqns.h"

struct spePoint {
	/*
	 * Saddle point equation point
	 * pf --final momentum
	 * ti -- time of ionization
	 * tf -- time of recollision
	 * k -- intermediate momentum
	 */
	spePoint(vec4 pf_);
	vec4 pf;
	std::vector<dcmplx> ti, tr, k;
	bool pointSolvedQ = false;

};

struct speGridParameters{
	int Nx1, Nz1, Nx2, Nz2;
	std::array<double,2> pfStart, pfEnd;
	double E01, E02; //Target parameters
};

class speGrid {
	/*
	 * class to solve saddle point equations over grid of 2x2D momentum
	 */
public:
	speGrid(const speGridParameters& gridParam, const laserField LF);
	virtual ~speGrid();
	void populateGrid();

	void solvePointRandom(size_t ix1, size_t iz1, size_t ix2, size_t iz2);
	void solvePointRandom(spePoint& point);
	void solveAllRandom();
	void solveWithAdjacent(spePoint& pointSolved, spePoint& pointToSolve);
	void propagateSolutionOverGrid(size_t ix1, size_t iz1, size_t ix2, size_t iz2);

	void printToFile(std::string outputFilename);


	spePoint& at(size_t ix1, size_t iz1, size_t ix2, size_t iz2);
	//get and set methods
	const auto& getSaddlePointEquations() const;
	const auto& getGrid() const;
	const auto& getPfStart() const;
	const auto& getPfEnd() const;
	const auto getNx1() const, getNz1() const, getNx2() const, getNz2() const;
private:
	size_t Nx1, Nz1, Nx2, Nz2;
	std::array<double,2> pfStart, pfEnd;
	double ddpx1, ddpz1, ddpx2, ddpz2;
	speEqns saddlePointEquations;
	std::vector<spePoint> grid;
};

#endif /* SPEGRID_H_ */
