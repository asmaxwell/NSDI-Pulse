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

using vec4 = std::array<double,4>;
using dcmplx = std::complex<double>;

class fields {
public:
	fields();
	virtual ~fields();
};

#endif /* FIELDS_H_ */
