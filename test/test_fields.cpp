/*
 * test_fields.cpp
 *
 *  Created on: 08 Feb 2024
 *  Tests for fields.cpp
 *      Author: andy
 */

#include <iostream>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include "../src/fields.h"

using namespace fields;
const double pi = 3.1415926535897932384626433832795028841971693993751;
constexpr dcmplx I(0,1.);

/// --- Testing Suite --- ///
TEST_CASE( "Testing SPEGrid Methods" ,"[fields]") {
	double omega=0.0551, rtUp=std::sqrt(1.2);
	auto phi = 0.33*pi;//0.33*pi;//GENERATE(0, 0.33*pi);
	SECTION("Testing Monochromatic Case"){
		int N = 1;
		fieldTypes fieldType = fieldTypes::monochromatic;
		///check constructor
		fields::monochromaticField LF(omega, rtUp, phi, N, fieldType);
		//check internal parameters
		REQUIRE(LF.getOmega() == omega);
		REQUIRE(LF.getrtUp() == rtUp);
		REQUIRE(LF.getPhi() == phi);
		REQUIRE(LF.getNCycles() == N);
		REQUIRE(LF.getFieldType() == fieldType);

		//check method values
		REQUIRE(LF.Afield(0.0) == 2*rtUp*std::cos(phi));
		REQUIRE(LF.Efield(0.0) == 2*rtUp*omega*std::sin(phi));

		double Up =rtUp*rtUp;
		double out1 = (Up*(2*(1.2566370614359172 + phi) + std::sin(2*(1.2566370614359172 + phi))))/omega;
		REQUIRE(std::abs(LF.A2Ifield(0.4*pi/omega)-out1) <1e-3);

		dcmplx out2 = (Up*(2.*(std::complex<double>(3.141592653589793, 0.7853981633974483) + phi) 
					+ std::sin(2.*(std::complex<double>(3.141592653589793, 0.7853981633974483) + phi))))/omega;
		REQUIRE(std::abs(LF.A2Ifield(std::complex<double>(pi/omega,0.25*pi/omega))-out2) <1e-3 );

		dcmplx out3 = (Up*(2.*(std::complex<double>(4.71238898038469,0.47123889803846897) + phi) 
					+ std::sin(2.*(std::complex<double>(4.71238898038469,0.47123889803846897) + phi))))/omega;
		REQUIRE(std::abs(LF.A2Ifield(std::complex<double>(1.5*pi/omega,0.15*pi/omega))-out3) <1e-3 );
	}
	SECTION("Testing Sin2 Case"){
		auto N = GENERATE(2, 3);
		fieldTypes fieldType = fieldTypes::sin2;
		///check constructor
		fields::sin2 LF(omega, rtUp, phi, N, fieldType);
		//check internal parameters
		REQUIRE(LF.getOmega() == omega);
		REQUIRE(LF.getrtUp() == rtUp);
		REQUIRE(LF.getPhi() == phi);
		REQUIRE(LF.getNCycles() == N);
		REQUIRE(LF.getFieldType() == fieldType);

		//check method values
		REQUIRE(LF.Afield(0.0) == 0.0);
		REQUIRE(LF.Efield(0.0) == 0.0);
		REQUIRE(LF.F(0.0) == 0.0);
		REQUIRE(LF.F(pi*N/omega) == 2*rtUp);

		if(std::abs(phi)<1e-12){
			if(N==2){
				REQUIRE(std::abs(LF.A2Ifield(0.4*pi/omega)-0.053642) <1e-3);
				REQUIRE(std::abs(LF.A2Ifield(std::complex(pi/omega,0.25*pi/omega))-std::complex(1.44196,19.7848)) <1e-3 );
				REQUIRE(std::abs(LF.A2Ifield(std::complex(1.5*pi/omega,0.15*pi/omega))-std::complex(37.2047,-2.33323)) <1e-3 );
			}
		}else if(std::abs(phi-0.33*pi)<1e-12){
			if(N==3){
				REQUIRE(std::abs(LF.A2Ifield(0.4*pi/omega)-0.0362568) <1e-3);
				REQUIRE(std::abs(LF.A2Ifield(std::complex(pi/omega,0.25*pi/omega))-std::complex(3.72291,1.72139)) <1e-3 );
				REQUIRE(std::abs(LF.A2Ifield(std::complex(1.5*pi/omega,0.15*pi/omega))-std::complex(4.52315,7.43183)) <1e-3 );
			}
			if(N==1){
				REQUIRE(std::abs(LF.A2Ifield(0.4*pi/omega)+4.25528) <1e-3);
				REQUIRE(std::abs(LF.A2Ifield(std::complex(pi/omega,0.25*pi/omega))-std::complex(101.636,10.4421)) <1e-3 );
				REQUIRE(std::abs(LF.A2Ifield(std::complex(1.5*pi/omega,0.15*pi/omega))-std::complex(84.4329,8.83016)) <1e-3 );
			}
		}

	}

}