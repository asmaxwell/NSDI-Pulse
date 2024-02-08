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
	/*
	 * Method to compute the amplitude exp(i S) and determinant
	 */
	std::vector<dcmplx> AmpOut;
	for(auto S : actPoint.SList){
		dcmplx expS = std::exp(std::complex(0,S));
		AmpOut.push_back(expS);
	}

	return AmpOut;

}

std::vector<dcmplx> amplitudeData::computeSPADet(const actionDataPoint &actPoint){
	/*
	 * Method to compute the amplitude exp(i S) and determinant
	 */
	std::vector<dcmplx> AmpOut;
	for(auto S : actPoint.SList){
		dcmplx expS = std::exp(std::complex(0,S));
		AmpOut.push_back(expS);
	}

	return AmpOut;

}

