/*
 * amplitudeData.h
 *
 *  Created on: 19 Dec 2023
 *      Author: andy
 */

#ifndef AMPLITUDEDATA_H_
#define AMPLITUDEDATA_H_

#include "actionData.h"

class amplitudeData {
public:
	amplitudeData(const actionData& actData);
	virtual ~amplitudeData();

	std::vector<dcmplx> computeSPAAmplitudes(const actionDataPoint &actPoint);
	std::vector<dcmplx> computeSPADet(const actionDataPoint &actPoint);
private:
	const actionData& actData;
};

#endif /* AMPLITUDEDATA_H_ */
