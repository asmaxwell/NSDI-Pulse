/*
 * multiroot.h
 *
 *  Created on: 15 Dec 2023
 *      Author: andy
 */

#ifndef MULTIROOT_H_
#define MULTIROOT_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "speEqns.h"

struct rparams
  {
	rparams();
	rparams(const vec4& pf_in, speEqns* SPEs_in)
	:pf(pf_in), SPEs(SPEs_in){}
    vec4 pf;
    speEqns* SPEs;
  };

class multiroot {
public:
	multiroot();
	virtual ~multiroot();
	static int solver(const gsl_vector * x, void *params, gsl_vector * f);
};

#endif /* MULTIROOT_H_ */
