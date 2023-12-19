/*
 * ActionData.h
 *
 *  Created on: 18 Dec 2023
 *      Author: andy
 */

#ifndef ACTIONDATA_H_
#define ACTIONDATA_H_

struct ActionDataPoint{
	ActionDataPoint();
	ActionDataPoint(vec3 pf_, dcmplx Stotal_);
	vec3 pf;
	std::vector<dcmplx> SList;
};
using actionGrid = std::vector<ActionDataPoint>;

class ActionData {
public:
	ActionData(speGrid &saddlePointGrid);
	virtual ~ActionData();

	std::vector<dcmplx> calculateActionPoint();
	void calculateAllAction();

private:
	const std::vector<spePoint> &grid;
	const speGrid &saddlePointGrid;
	const laserField &LF;
	actionGrid gridS;
};

#endif /* ACTIONDATA_H_ */
