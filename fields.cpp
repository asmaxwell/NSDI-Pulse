/*
 * fields.cpp
 *
 *  Created on: 12 Dec 2023
 *      Author: andy
 */

#include <iostream>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include "fields.h"

using namespace fields;
const double pi = 3.1415926535897932384626433832795028841971693993751;
constexpr dcmplx I(0,1.);

//// Base laser Field Class definitions
fields::laserField::laserField():omega(0.), rtUp(0.), phi(0.), NCycles(0){
Up=0;
}

fields::laserField::laserField(double omega_, double rtUp_, double phi_, int NCycles_, fieldTypes fieldType_):
		omega(omega_), rtUp(rtUp_), phi(phi_), NCycles(NCycles_), fieldType(fieldType_){
	Up=rtUp*rtUp;
}

fields::laserField::~laserField(){};


//Get set methods
const double fields::laserField::getOmega() const {return omega;}
const double fields::laserField::getrtUp() const {return rtUp;}
const double fields::laserField::getPhi() const {return phi;}
const int fields::laserField::getNCycles() const {return NCycles;}
const fieldTypes fields::laserField::getFieldType() const{return fieldType;}

//Generic Laser field functions
//Vector Potential
double fields::laserField::Afield(double t) const{
	return F(t)*std::cos(omega*t + phi);
}
double fields::laserField::A2field(double t) const{
	return 2*rtUp*std::cos(omega*t + phi);
}
dcmplx fields::laserField::Afield(dcmplx t) const{
	return F(t)*std::cos(omega*t + phi);
}
//Derivative of vector potential
dcmplx fields::laserField::dtAfield(dcmplx t) const{
	return -Efield(t);
}
//Electric Field
double fields::laserField::Efield(double t) const{
	return F(t)*omega*std::sin(omega*t+phi) - dF(t)*std::cos(omega*t+phi);
}
dcmplx fields::laserField::Efield(dcmplx t) const{
		return F(t)*omega*std::sin(omega*t+phi) - dF(t)*std::cos(omega*t+phi);
}
//double fields::laserField::AI2field(double t) const{
//	return (2*rtUp/omega)*std::sin(omega*t+phi);
//}
//
//dcmplx fields::laserField::AI2field(dcmplx t) const{
//	return (2*rtUp/omega)*std::sin(omega*t+phi);
//}

//// Monochromatic Field definitions
//double fields::monochromatic::A_end(double tf) const{ return Afield(tf);}
double fields::monochromaticField::AIfield(double t) const{
	return (2*rtUp/omega)*std::sin(omega*t+phi);
}

dcmplx fields::monochromaticField::AIfield(dcmplx t) const{
	return (2*rtUp/omega)*std::sin(omega*t+phi);
}

double fields::monochromaticField::A2Ifield(double t) const{
	return (Up/omega)*( 2*(phi+omega*t)+std::sin(2*(omega*t+phi)) );
}
dcmplx fields::monochromaticField::A2Ifield(dcmplx t) const{
	return (Up/omega)*( 2.*(phi+omega*t)+std::sin(2.*(omega*t+phi)) );
}


////Sin^2 Pulse
double fields::sin2::F(double t) const{
	double return_val = 0.;
		if(t>0 && t< NCycles*2*pi/omega + 1e-15){
			return_val = 2*rtUp*std::pow(std::sin(omega*t/(2*NCycles)),2);
		}
		return return_val;
}
dcmplx fields::sin2::F(dcmplx t) const{
	dcmplx return_val = 0.;
	double tr = t.real();
		if(tr>0 && tr<NCycles*2*pi/omega + 1e-15){
			return_val = 2*rtUp*std::pow(std::sin((omega/(2*NCycles))*t),2);
		}
		return return_val;
}
double fields::sin2::dF(double t) const{
	double return_val = 0.;
		if(t>0 && t<NCycles*2*pi/omega + 1e-15){
			return_val = (2*rtUp*omega/NCycles)*std::sin(omega*t/(2*NCycles))*std::cos(omega*t/(2*NCycles));
		}
		return return_val;
}
dcmplx fields::sin2::dF(dcmplx t) const{
	dcmplx return_val = 0.;
	double tr = t.real();
		if(tr>0 && tr<NCycles*2*pi/omega + 1e-15){
			return_val = (2*rtUp*omega/NCycles)*std::sin((omega/(2*NCycles))*t)*std::cos((omega/(2*NCycles))*t);
		}
		if(tr>NCycles*2*pi/omega){std::cerr<<"Error complex field can not be called with a real time outside of pulse length\n";}
		return return_val;
}

//imported from linear SFA Cython code
double fields::sin2::AIfield( double t) const{
	//if(fieldType==fieldTypes::sin2 && NCycles==1){return 0.0;}
	double factor, a1, a2, a3;//, AI0=0.;
	//Integral of vector potential
	double s1, s2, s3;
	//Turn field off after end of pulse
	if(t>2*NCycles*pi/omega+1e-15){return AIfield(2*NCycles*pi/omega);}
	//Case if N==1 as general expression will diverge
	if(NCycles==1){
		//AI0 = (3*rtUp*std::sin(phi))/(4.*omega);
		return -(rtUp*(2*t*omega*std::cos(phi) - 4*std::sin(t*omega + phi) + std::sin(2*t*omega + phi)))/(4.*std::sqrt(2)*omega);// - AI0;
	}
	else{
		//AI0 = (rtUp*std::sin(phi))/(omega - NCycles*NCycles*omega); // AfI(0) to ensure limits are correct.
		a1 = NCycles*NCycles - 1;
		a2 = (NCycles + 1.)/NCycles;
		a3 = (NCycles - 1.)/NCycles;
		factor = rtUp/(2*omega*a1);

		s1 = std::sin(omega*t + phi);
		s2 = std::sin(phi + a2 *omega*t);
		s3 = std::sin(phi + a3 *omega*t);

		return factor * (2*a1*s1 - NCycles*NCycles *( a3*s2 + a2*s3) );// - AI0;
	}
}

dcmplx fields::sin2::AIfield( dcmplx t) const{
	//if(fieldType==fieldTypes::sin2 && NCycles==1){return 0.0;}
	//Integral of vector potential
	double factor, a1, a2, a3, tr=t.real();
	dcmplx s1, s2, s3;//, AI0;
	//Turn field off after end of pulse
	if(tr>2*NCycles*pi/omega+1e-15){
//		std::cerr<<"Error complex field should not be evaluated outside pulse\n";
		return AIfield(2*NCycles*pi/omega);
	}
	//Case if N==1 as general expression will diverge
	if(NCycles==1){
		//AI0 = (3*rtUp*std::sin(phi))/(4.*omega);
		return -(rtUp*(2.*t*omega*std::cos(phi) - 4.*std::sin(t*omega + phi) + std::sin(2.*t*omega + phi)))/(4.*std::sqrt(2.)*omega);// - AI0;
	}
	else{
		//AI0 = (rtUp*std::sin(phi))/(omega - NCycles*NCycles*omega); // AfI(0) to ensure limits are correct.
		a1 = NCycles*NCycles - 1;
		a2 = (NCycles + 1.)/NCycles;
		a3 = (NCycles - 1.)/NCycles;
		factor = rtUp/(2*omega*a1);

		s1 = std::sin(omega*t + phi);
		s2 = std::sin(phi + a2 *omega*t);
		s3 = std::sin(phi + a3 *omega*t);

		return factor * (2.*a1*s1 - static_cast<double>(NCycles*NCycles) *( a3*s2 + a2*s3) );// - AI0;
	}
}

double fields::sin2::A2Ifield( double t) const{
	//if(fieldType==fieldTypes::sin2 && NCycles==1){return 0.0;}
	if(t>2*NCycles*pi/omega+1e-15){return AIfield(2*NCycles*pi/omega);}
	//N=1 case
	if(NCycles==1){
		double fac = Up/(16*omega), lin = 12*omega*t;
		double term1 = 2*omega*t*std::cos(2*phi);
		double term2 = -16.*std::sin(omega*t);
		double term3 = 2*std::sin(2*omega*t);
		double term4 = 6*std::sin(2*(phi+omega*t));
		double term5 = -8*std::sin(2*phi+omega*t);
		double term6 = 0.5*std::sin(2*phi+2*omega*t);
		double term7 = (8./3.)*std::sin(2*phi+3*omega*t);

		return fac * (lin+term1+term2+term3+term4+term5+term6+term7);
	}else{
		//General case
		double fac = Up/(16*omega), lin = 12*t*omega;
		double term1 = 6*std::cos(2*t*omega)*std::sin(2*phi);
		double term2 = 6*std::cos(2*phi)*std::sin(2*t*omega);
		double term3 = -16*NCycles*std::sin(omega*t/NCycles);
		double term4 = 2*NCycles*std::sin(2*omega*t/NCycles);
		double term5 = -8*NCycles*std::sin(2*phi+(2+1./NCycles)*omega*t)/(2*NCycles+1);
		double term6 = -8*NCycles*std::sin(2*phi+(2-1./NCycles)*omega*t)/(2*NCycles-1);
		double term7 = NCycles*std::sin(2*(phi+(NCycles-1)*omega*t/NCycles))/(NCycles-1);
		double term8 = NCycles*std::sin(2*(phi+(NCycles+1)*omega*t/NCycles))/(NCycles+1);

		std::cout<<fac<<" "<<lin<<" "<<term1<<" "<<term2<<" "<<term3<<" "<<term4<<" "<<term5<<" "
				<<" "<<term7<<" "<<" "<<term8<<" "<<" "<<term6<<"\n";

		return fac * (lin +term1+term2+term3+term4+term5+term6+term7+term8);
	}

}

dcmplx fields::sin2::A2Ifield( dcmplx t) const{
	//if(fieldType==fieldTypes::sin2 && NCycles==1){return 0.0;}
	double tr = t.real();
	if(tr>2*NCycles*pi/omega+1e-15){return AIfield(2*NCycles*pi/omega);}
	if(NCycles==1){
		double fac = Up/(16*omega);
		dcmplx lin = 12.*omega*t;
		dcmplx term1 = 2.*omega*t*std::cos(2.*phi);
		dcmplx term2 = -16.*std::sin(omega*t);
		dcmplx term3 = 2.*std::sin(2.*omega*t);
		dcmplx term4 = 6.*std::sin(2.*(phi+omega*t));
		dcmplx term5 = -8.*std::sin(2.*phi+omega*t);
		dcmplx term6 = 0.5*std::sin(2.*phi+2*omega*t);
		dcmplx term7 = (8./3.)*std::sin(2.*phi+3.*omega*t);

		return fac * (lin+term1+term2+term3+term4+term5+term6+term7);
	}else{
		//General case
		double fac = Up/(16*omega);
		dcmplx lin = 12.*t*omega;
		dcmplx term1 = 6.*std::cos(2.*t*omega)*std::sin(2.*phi);
		dcmplx term2 = 6.*std::cos(2.*phi)*std::sin(2.*t*omega);
		dcmplx term3 = -16.*NCycles*std::sin((omega/NCycles)*t);
		dcmplx term4 = 2.*NCycles*std::sin((2*omega/NCycles)*t);
		dcmplx term5 = -8.*NCycles*std::sin(2*phi+(2+1./NCycles)*omega*t)/(2.*NCycles+1.);
		dcmplx term6 = -8.*NCycles*std::sin(2*phi+(2-1./NCycles)*omega*t)/(2.*NCycles-1.);
		dcmplx term7 = static_cast<double>(NCycles)*std::sin(2.*(phi+(NCycles-1.)/NCycles*omega*t))/(NCycles-1.);
		dcmplx term8 = static_cast<double>(NCycles)*std::sin(2.*(phi+(NCycles+1.)/NCycles*omega*t))/(NCycles+1.);

		std::cout<<fac<<" "<<lin<<" "<<term1<<" "<<term2<<" "<<term3<<" "<<term4<<" "<<term5<<" "
				<<" "<<term7<<" "<<" "<<term8<<" "<<" "<<term6<<"\n";

		return fac * (lin +term1+term2+term3+term4+term5+term6+term7+term8);
	}

}

/// --- Testing Suite --- ///
TEST_CASE( "Testing SPEGrid Methods" ,"[fields]") {
	double omega=0.0551, rtUp=std::sqrt(1.2);
	auto phi = GENERATE(0, 0.33*pi);
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
		if(N==1 && std::abs(phi-0.33*pi)<1e-12){
			REQUIRE(std::abs(LF.A2Ifield(0.4*pi/omega)-33.1287) <1e-3);
			REQUIRE(std::abs(LF.A2Ifield(std::complex(pi/omega,0.25*pi/omega))-std::complex(184.726 ,10.0647)) <1e-3 );
			REQUIRE(std::abs(LF.A2Ifield(std::complex(1.5*pi/omega,0.15*pi/omega))-std::complex(177.051 ,31.9446 )) <1e-3 );
		}

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

