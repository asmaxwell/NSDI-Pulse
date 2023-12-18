/*
 * multiroot.cpp
 *
 *  Created on: 15 Dec 2023
 *      Author: andy
 */

#include "multiroot.h"

multiroot::multiroot() {

}

multiroot::~multiroot() {

}

int multiroot::solver(const gsl_vector * x, void *params, gsl_vector * f){
	double x0 = gsl_vector_get (x, 0);
	double x1 = gsl_vector_get (x, 1);
	double x2 = gsl_vector_get (x, 2);
	double x3 = gsl_vector_get (x, 3);

	speEqns* SPEs = ((struct rparams *) params)->SPEs;
	const vec4& pf = ((struct rparams *) params)->pf;

	dcmplx ti = std::complex(x0,x1), tr = std::complex(x2,x3);
	const dcmplx y0 = SPEs->SPE_ti(ti, tr, pf);
	const dcmplx y1 = SPEs->SPE_tr(ti, tr, pf);

	gsl_vector_set (f, 0, y0.real());
	gsl_vector_set (f, 1, y0.imag());
	gsl_vector_set (f, 2, y1.real());
	gsl_vector_set (f, 3, y1.imag());
	return GSL_SUCCESS;
}

