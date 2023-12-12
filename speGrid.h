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

class spePoint {
	vec6 p;
	dcmplx ti, tr, k;

};

class speGrid {
public:
	speGrid();
	virtual ~speGrid();

	std::vector<spePoint> grid;
};

#endif /* SPEGRID_H_ */
