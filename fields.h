/*
 * fields.h
 *
 *  Created on: 12 Dec 2023
 *      Author: andy
 */

#ifndef FIELDS_H_
#define FIELDS_H_

#include <array>
#include <complex>

using vec6 = std::array<double,6>;
using dcmplx = std::complex<double>;

class fields {
public:
	fields();
	virtual ~fields();
};

#endif /* FIELDS_H_ */
