/*
 * ActionData.cpp
 *
 *  Created on: 18 Dec 2023
 *      Author: andy
 */
#include "actionData.h"

#include <fstream>
#include <iomanip>
#include <iostream>


actionDataPoint::actionDataPoint(vec4 pf_, std::vector<dcmplx> SList_) : pf(pf_), SList(SList_){};

actionData::actionData(speGrid &saddlePointGrid_) : saddlePointGrid(saddlePointGrid_), grid(saddlePointGrid_.getGrid())
, LF(saddlePointGrid_.getSaddlePointEquations().getLaserField())
, E01(saddlePointGrid_.getSaddlePointEquations().getE01())
, E02(saddlePointGrid_.getSaddlePointEquations().getE02()){
	//Fill action grid from saddle point grid

}

std::vector<dcmplx> actionData::calculateActionPoint(spePoint & saddlePoint){
	std::vector<dcmplx> SList={};

	std::vector<dcmplx> tiList = saddlePoint.ti, trList = saddlePoint.tr, kList = saddlePoint.k;
	vec4 pf = saddlePoint.pf;
	for(auto i=0; i!=tiList.size(); ++i){
		dcmplx ti = tiList[i], tr =trList[i], k = kList[i];
		dcmplx AIti = LF->AIfield(ti), AItr = LF->AIfield(tr);
		dcmplx A2Iti = LF->A2Ifield(ti), A2Itr = LF->A2Ifield(tr);

		dcmplx ETerm = E01*ti + E02*tr;

		dcmplx p1Term = 0.5*((pf[0]*pf[0]+pf[1]*pf[1])*tr + 2*pf[1]*AItr + A2Itr);
		dcmplx p2Term = 0.5*((pf[2]*pf[2]+pf[3]*pf[3])*tr + 2*pf[3]*AItr + A2Itr);

		dcmplx kTerm = -0.5*( k*k*(tr-ti) + 2.*k*(AItr-AIti) + A2Itr);

		SList.push_back(ETerm + p1Term + p2Term + kTerm);
	}

	return SList;

}

void actionData::calculateAllAction(){
	#pragma omp parallel num_threads(8)
	#pragma omp for ordered
	for(auto point : grid){
		vec4 pf = point.pf;
		std::vector<dcmplx> SList = calculateActionPoint(point);
		actionDataPoint actPoint(pf, SList);
		gridS.push_back(actPoint);
	}
}

void actionData::printToFile(std::string outputFilename){
	std::ofstream outputFile(outputFilename);
	if(outputFile.fail()){
		std::cerr<<"Output file not opened\n";
	}
	for(const auto& point : gridS){
		vec4 pf = point.pf;
		std::vector<dcmplx> SList = point.SList;
		outputFile<<pf[0]<<" "<<pf[1]<<" "<<pf[2]<<" "<<pf[3]<<" ";
		for(auto S : SList){
			outputFile<<S.real()<<" "<<S.imag()<<" ";
		}
		outputFile<<"\n";
	}
}

actionData::~actionData() {};

