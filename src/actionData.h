/*
 * ActionData.h
 *
 *  Created on: 18 Dec 2023
 *      Author: andy
 */

#ifndef ACTIONDATA_H_
#define ACTIONDATA_H_

#include "speGrid.h"

struct actionDataPoint{
	actionDataPoint();
	actionDataPoint(vec4 pf_, std::vector<dcmplx> SList_);
	vec4 pf;
	std::vector<dcmplx> SList;
};
using actionGrid = std::vector<actionDataPoint>;

class actionData {
public:
	actionData(speGrid &saddlePointGrid);
	virtual ~actionData();

	std::vector<dcmplx> calculateActionPoint(spePoint & saddlePoint);
	void calculateAllAction();
	void printToFile(std::string outputFilename);

private:
	const std::vector<spePoint> &grid;
	const speGrid &saddlePointGrid;
	const laserField &LF;
	const double E01, E02;
	actionGrid gridS;
};

#endif /* ACTIONDATA_H_ */
