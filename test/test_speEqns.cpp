/*
 * test_speEqns.cpp
 *
 *  Created on: 08 Feb 2024
 *      Author: andy
 */
#include <iostream>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include "../src/speEqns.h"

const double pi = 3.1415926535897932384626433832795028841971693993751;
constexpr dcmplx I(0,1.);
/// --- Testing Suite --- ///
TEST_CASE( "Testing SPEEqns Methods" ,"[SPEEqns]") {
	//set up code
	double omega=0.057, rtUp=std::sqrt(1.2), E01 = 0.79, E02 = 1.51;
	double phi = GENERATE(0.0, 0.33*pi);
	size_t N = GENERATE(1, 2);
	fieldTypes fieldType = GENERATE(fieldTypes::monochromatic, fieldTypes::sin2);
	///check constructor
	laserField LF;
	if(fieldType == monochromatic){
		LF = std::make_shared<fields::monochromaticField>(omega, rtUp, phi, N, fieldType);
	}else{
		LF = std::make_shared<fields::sin2>(omega, rtUp, phi, N, fieldType);
	}
	speEqns saddlePointEqns(E01, E02, LF);
	const std::vector<vec4> pfList = {{0.,0.,0.,0.}, {0., 1., 0., 1}};

	//test equal ti and tr
	SECTION("test equal ti and tr"){
		for(auto pf : pfList){
			auto ti = GENERATE(std::complex(0.,0.), std::complex (1., 0.5));
			auto tr = ti;
			REQUIRE(saddlePointEqns.k(ti,tr) == 0.0);

			REQUIRE(saddlePointEqns.SPE_ti(ti, tr, pf) == std::pow(LF->Afield(ti),2)+2* E01 );
			dcmplx RHS_SPE_tr = std::pow(pf[1]+LF->Afield(tr),2)+pf[0]*pf[0]
							  + std::pow(pf[3]+LF->Afield(tr),2)+pf[2]*pf[2] - pow(LF->Afield(tr),2) + 2*E02;
			REQUIRE(saddlePointEqns.SPE_tr(ti, tr, pf) == RHS_SPE_tr );
		}
	}
	SECTION("test near saddles"){
		for(auto pf : pfList){
			//tyring to be ner saddle points in this test
			auto ti = std::complex(0.25, 0.25)*2.*pi/omega;
			auto tr = std::complex(0.75, 0.05)*2.*pi/omega;
			dcmplx k = saddlePointEqns.k(ti,tr);

			REQUIRE(saddlePointEqns.SPE_ti(ti, tr, pf) == std::pow(k + LF->Afield(ti),2)+2* E01 );
			dcmplx RHS_SPE_tr = std::pow(pf[1]+LF->Afield(tr),2)+pf[0]*pf[0]
							  + std::pow(pf[3]+LF->Afield(tr),2)+pf[2]*pf[2] - pow(k + LF->Afield(tr),2) + 2*E02;
			REQUIRE(saddlePointEqns.SPE_tr(ti, tr, pf) == RHS_SPE_tr );
		}
	}
	//Test one off cases
	if( (fieldType == monochromatic) && (N==1) && std::abs(phi)<1e-12){
		std::cout<<"Mono Test\n";
		dcmplx ti(pi*0.6/omega, pi*0.1/omega), tr(pi*1.8/omega, pi*0.25/omega);
		auto err = saddlePointEqns.solveRoot(ti, tr, {0.0, -5, 0.0, -5});
		std::cout<<omega*ti/(2*pi)<<", "<<omega*tr/(2*pi)<<"\n";
		REQUIRE(std::abs(omega*ti/(2*pi) - std::complex(0.287772 , 0.052751 ))<1e-3 );//values taken from Mathematica
		REQUIRE(std::abs(omega*tr/(2*pi) - std::complex(0.979982 , 0.126481))<1e-3 );
	}
	if( (fieldType == sin2) && (N==2) && std::abs(phi-0.33*pi)<1e-12 ){
		std::cout<<"Sin2 Test\n";
			dcmplx ti(pi*0.6/omega, pi*0.25/omega), tr(pi*2.2/omega, pi*0.034/omega);
			auto err = saddlePointEqns.solveRoot(ti, tr, {0.0, -2, 0.0, -2});
			std::cout<<omega*ti/(2*pi)<<", "<<omega*tr/(2*pi)<<"\n";
			REQUIRE(std::abs(omega*ti/(2*pi) - std::complex(0.255135 , 0.21564 ))<1e-3 );//values taken from Mathematica
			REQUIRE(std::abs(omega*tr/(2*pi) - std::complex(0.903223 , -0.0170268))<1e-3 );
		}

	//Test RandSolve Method
	SECTION("RandSolve"){
		size_t NumRandGuesses = 500;
		for(auto pf : pfList){
			std::vector<timePair> timePairList = saddlePointEqns.RandSolve(pf, NumRandGuesses);
//			std::cout<<"len = "<<timePairList.size()
//								<<", NCycles = "<<N<<", phi = "<<phi
//								<<", pf = ("<<pf[0]<<", "<<pf[1]<<", "<<pf[2]<<", "<<pf[3]<<")\n";
			for(auto tsPair : timePairList){
				dcmplx ti = tsPair[0], tr = tsPair[1];
//				if(N==1)std::cout<<"ti = "<<omega*ti/(2*pi)<<", tr = "<<omega*tr/(2*pi)<<"\n";
				REQUIRE(ti.real() < tr.real());
				REQUIRE(ti.imag()>0);
				REQUIRE(ti.real()>0);
				REQUIRE(tr.real()>0);
				REQUIRE(ti.real()<2*pi*N/omega);

				dcmplx SPE1val = saddlePointEqns.SPE_ti(ti, tr, pf);
				dcmplx SPE2val = saddlePointEqns.SPE_tr(ti, tr, pf);
				REQUIRE(std::abs(SPE1val+SPE2val) < 100 * 1e-12);

			}
			//check for duplicates
			size_t timePairListSize = timePairList.size();
			if(timePairListSize>0){
				int i = rand()%timePairList.size(), j = rand()%timePairList.size();
				if(i==j){
					continue;
				}else{
					double tiErr = std::abs(timePairList[j][0]-timePairList[i][0]);
					double trErr = std::abs(timePairList[j][1]-timePairList[i][1]);
//					std::cout<<"tiErr = "<<tiErr<<", trErr = "<<trErr<<"\n";
					REQUIRE(tiErr + trErr > 10 * 1e-12);
				}
			}
		}

	}


}