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

#include "fields.h"

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
	dcmplx ti, tr, k;

};

struct speGridParameters{
	int Nx1, Nz1, Nx2, Nz2;
	std::array<double,2> pfStart, pfEnd;
};

class speGrid {
public:
	speGrid(const speGridParameters& gridParam);
	virtual ~speGrid();
	void populateGrid();

	//get and set methods
	const auto& getGrid() const;
	const auto& getPfStart() const;
	const auto& getPfEnd() const;
	const auto getNx1() const, getNz1() const, getNx2() const, getNz2() const;
private:
	size_t Nx1, Nz1, Nx2, Nz2;
	std::array<double,2> pfStart, pfEnd;
	double ddpx1, ddpz1, ddpx2, ddpz2;
	std::vector<spePoint> grid;
};

#endif /* SPEGRID_H_ */
