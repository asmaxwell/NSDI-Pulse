/*
 * speEqns.h
 *
 *  Created on: 12 Dec 2023
 *      Author: andy
 */

#ifndef SPEEQNS_H_
#define SPEEQNS_H_

#include <vector>

#include "fields.h"
using timePair = std::array<dcmplx,2>;

class speEqns {
	/*
	 * class to solve the saddle point equations given by SPE_ti and SPE_tr, which can provide
	 * solutions to complex variables ti and tr.
	 */
public:
	speEqns(double E01, double E02, const laserField& LF);
	virtual ~speEqns();

	dcmplx k(dcmplx ti, dcmplx tr) const;
	dcmplx SPE_ti(dcmplx ti, dcmplx tr, const vec4& pf) const;
	dcmplx SPE_tr(dcmplx ti, dcmplx tr, const vec4& pf) const;

	int solveRoot(dcmplx& ti, dcmplx& tr, const vec4& pf);
	std::vector<timePair> RandSolve(const vec4& pf, size_t NumRandGuesses);

	const auto & getLaserField() const;
private:
	const double E01, E02;
	const double rel_err = 1e-12;
	const laserField LF;

};

#endif /* SPEEQNS_H_ */
