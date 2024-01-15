/*
 * amplitudeData.cpp
 *
 *  Created on: 19 Dec 2023
 *      Author: andy
 */

#include "amplitudeData.h"

amplitudeData::amplitudeData(const actionData& actData_) actData(actData_) {
	// intialize
}

amplitudeData::~amplitudeData() {};

std::vector<dcmplx> amplitudeData::computeSPAAmplitudes(const actionDataPoint &actPoint){
	for(auto S : actPoint.SList){
		dcmplx expS = std::exp(std::complex(0,S));
		dcmplx detS = ()
	}

}

