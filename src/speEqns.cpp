/*
 * speEqns.cpp
 *
 *  Created on: 12 Dec 2023
 *      Author: andy
 */
#include <algorithm>
#include <iostream>
#include <random>

#include "multiroot.h"
#include "speEqns.h"


const double pi = 3.1415926535897932384626433832795028841971693993751;

speEqns::speEqns(double E01_, double E02_, const laserField& LF_) : E01(E01_), E02(E02_), LF(LF_) {};

speEqns::~speEqns() {};

dcmplx speEqns::k(dcmplx ti, dcmplx tr) const{
	if(std::abs(tr-ti)<1e-16){return 0.0;}
	return - ( LF->AIfield(tr) - LF->AIfield(ti) ) / (tr - ti);
}
dcmplx speEqns::SPE_ti(dcmplx ti, dcmplx tr, const vec4& pf) const{
	dcmplx v = k(ti,tr)+LF->Afield(ti);
	return v*v + 2*E01;
}
dcmplx speEqns::SPE_tr(dcmplx ti, dcmplx tr, const vec4& pf) const{
	dcmplx Aftr = LF->Afield(tr);
	dcmplx v1 = pf[1]+Aftr, v2 = pf[3]+Aftr, v3 = k(ti,tr)+Aftr;
	return v1*v1 + pf[0]*pf[0] + v2*v2 + pf[2]*pf[2] - v3*v3 + 2*E02;
}


int speEqns::solveRoot(dcmplx& ti, dcmplx& tr, const vec4& pf){
	/*
	 * Root finding solver for ionization and recollsion times
	 * Solved using GSL
	 */
	timePair timesOut={};
	const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *s;
	int status, errFlag=0;
	size_t iter = 0, max_iter = 100;
	const double Terr = rel_err;
	const size_t n = 4;
	struct rparams params(pf, this);
	gsl_multiroot_function f = {&multiroot::solver, n, &params};
	gsl_vector *x = gsl_vector_alloc(n);

	dcmplx tig = ti, trg = tr;

	gsl_vector_set (x, 0, tig.real());
	gsl_vector_set (x, 1, tig.imag());
	gsl_vector_set (x, 2, trg.real());
	gsl_vector_set (x, 3, trg.imag());
	T = gsl_multiroot_fsolver_hybrid;
	s = gsl_multiroot_fsolver_alloc(T, n);

	gsl_multiroot_fsolver_set (s, &f, x);
	do
	{
	  iter++;
	  status = gsl_multiroot_fsolver_iterate (s);
	  if (status)
		break;
	  status =
		gsl_multiroot_test_residual (s->f, Terr);

	  if(std::isnan(gsl_vector_get(s->f,0))||std::isnan(gsl_vector_get(s->f,1))||std::isnan(gsl_vector_get(s->f,2))
		||std::isnan(gsl_vector_get(s->x,0))||std::isnan(gsl_vector_get(s->x,1))||std::isnan(gsl_vector_get(s->x,2))
	  )
		  {
			  errFlag = -2;
			  throw std::runtime_error("Error NAN encountered solver halted!\n");
		  }
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	double f0 = gsl_vector_get(s->f,0), f1 = gsl_vector_get(s->f,1), f2 = gsl_vector_get(s->f,2), f3 = gsl_vector_get(s->f,3);
	bool residualTest = (std::abs(f0) > rel_err) || (std::abs(f1) > rel_err) || (std::abs(f2) > rel_err) || (std::abs(f3) > rel_err);
//	if(residualTest){std::cout<<"f0 = "<<f0<<", f1 = "<<f1<<", f2 = "<<f2<<", f3 = "<<f3<<"\n";}
	if(status==GSL_ENOPROG || residualTest){
		errFlag = -1;
		//throw std::runtime_error("Error GSL multiroot solver status = GSL_ENOPROG!\n");
	}
	if(errFlag==0){
		ti = std::complex(gsl_vector_get(s->x,0),gsl_vector_get(s->x,1));
		tr = std::complex(gsl_vector_get(s->x,2), gsl_vector_get(s->x,3));
	}


	gsl_multiroot_fsolver_free (s);
	gsl_vector_free (x);

	return errFlag;
}

std::vector<timePair> speEqns::RandSolve(const vec4& pf, size_t NumRandGuesses){
	/*
	 * Function to generate random numbers in complex ti and tr to search full space of solutions
	 * to the saddle point equations.
	 */

	//make output time list
	std::vector<timePair> timeListOut;
	//field parameters
	double omega = LF->getOmega(), halfCycle = pi/omega;
	int NCycles = LF->getNCycles();

	std::random_device rd_tir, rd_tii, rd_trr, rd_tri;  // Will be used to obtain a seed for the random number engine
	std::mt19937 gentir(rd_tir()), gentii(rd_tii()), gentrr(rd_trr()), gentri(rd_tri()); // Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<double> tirRand(0., 2*NCycles*halfCycle), tiiRand(0., 0.5 * halfCycle)
										, triRand(-0.25*halfCycle, 0.25*halfCycle);
	for(auto i=0; i!=NumRandGuesses; ++i){
		double tir = tirRand(gentir), tii = tiiRand(gentii), tri = triRand(gentri);
		double trr_start = tir+0.5*halfCycle, trr_end = tir+2*halfCycle;

		if(LF->getFieldType() == fieldTypes::sin2){
			if(trr_start > 2*NCycles*halfCycle){continue;}
			trr_end = std::min(trr_end, 2*NCycles*halfCycle);
		}

		std::uniform_real_distribution<double> trrRand(trr_start, trr_end);
		double trr = trr = trrRand(gentrr);
		dcmplx ti = std::complex(tir, tii), tr = std::complex(trr, trr);
		int errFlag = solveRoot(ti, tr, pf);
		if(errFlag!=0){continue;}

		bool skipRandTimes = (ti.real() < 0) || (ti.real() > tr.real() ) || (ti.imag() < 0 )
				|| (ti.real() > 2*NCycles*halfCycle) || (std::abs(tr.real() - ti.real()) > 2.25*halfCycle)
				|| (std::abs(tr.real() - ti.real()) < 0.5*halfCycle) || (tr.real() > 2*(NCycles+0.75)*halfCycle);
		bool returnTimeInsidePulse = (LF->getFieldType() == fieldTypes::sin2) && (tr.real() > 2*NCycles*halfCycle);


		const timePair testPair = {ti, tr};
		const double err = rel_err;
		//lambda function to determine if a pair of times are a duplicate
		auto times_equal = [testPair, err](timePair tsPair) {
			double tiErr = std::abs(tsPair[0]-testPair[0]);
			double trErr = std::abs(tsPair[1]-testPair[1]);
			return (tiErr + trErr) < 100 * err;
		};
		//find a duplicate if it exists
		auto duplicate = std::find_if(timeListOut.begin(), timeListOut.end(), times_equal);
		//set bool to skip this solution if there is a duplicate
		bool equalExstingSolution = (duplicate!=timeListOut.end());
		if(skipRandTimes || returnTimeInsidePulse || equalExstingSolution || (errFlag!=0)){
			continue;
		}else{
			timeListOut.push_back({ti, tr});
		}

	}
	//decided to sort solutions as it may be useful, but can comment out if it is not required
	//lambda for sorting
	auto sortByIonizationReal = [](timePair x, timePair y){return x[0].real()<y[0].real();};
	std::sort(timeListOut.begin(), timeListOut.end(), sortByIonizationReal);

	return timeListOut;
}
const laserField & speEqns::getLaserField() const {return LF;}
const double speEqns::getE01() const {return E01;}
const double speEqns::getE02() const {return E02;}


