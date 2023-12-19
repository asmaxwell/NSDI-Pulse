/*
 * ActionData.cpp
 *
 *  Created on: 18 Dec 2023
 *      Author: andy
 */

#include "ActionData.h"

ActionData::ActionData(speGrid &saddlePointGrid_) : saddlePointGrid(saddlePointGrid_), grid(saddlePointGrid_.getGrid())
, LF(saddleGrid_.getSaddlePointEquations().getLaserField()){
	//Fill action grid from saddle point grid

}

std::vector<dcmplx> ActionData::calculateActionPoint(spePoint & saddlePoint){
	std::vector<dcmplx> SList={};

	std::vector<dcmplx> tiList = saddlePoint.ti, trList = saddlePoint.tr, kList = saddlePoint.k;
	vec4 pf = saddlePoint.pf;
	for(auto i=0; i!=tiList.size(); ++i){
		dcmplx ti = tiList[i], tr =trList[i], k = kList[i];
		dcmplx AIti = LF->AIfield(ti), AItr = LF->AIfield(tr);
		dcmplx A2Iti = LF->A2Ifield(ti), A2Itr = LF->A2Ifield(tr);

		dcmplx p1Term = 0.5*((pf[0]*pf[0]+pf[1]*pf[1])*tr + 2*pf[1]*AItr + A2Itr);
		dcmplx p2Term = 0.5*((pf[2]*pf[2]+pf[3]*pf[3])*tr + 2*pf[3]*AItr + A2Itr);

		dcmplx kTerm = -0.5*( k*k*(tr-ti) + 2*k*(AItr-AIti) + A2Itr);

		SList.push_back(p1Term + p2Term + kTerm);
	}

	return SList;

}

ActionData::~ActionData() {
	// TODO Auto-generated destructor stub
}

