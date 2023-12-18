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
fields::laserField::laserField():omega(0.), rtUp(0.), phi(0.), NCycles(0)
{};

fields::laserField::laserField(double omega_, double rtUp_, double phi_, int NCycles_, fieldTypes fieldType_):
		omega(omega_), rtUp(rtUp_), phi(phi_), NCycles(NCycles_), fieldType(fieldType_)
{};

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
double fields::laserField::AI2field(double t) const{
	return (2*rtUp/omega)*std::sin(omega*t+phi);
}

dcmplx fields::laserField::AI2field(dcmplx t) const{
	return (2*rtUp/omega)*std::sin(omega*t+phi);
}

//// Monochromatic Field definitions
//double fields::monochromatic::A_end(double tf) const{ return Afield(tf);}
double fields::monochromaticField::AIfield(double t) const{
	return (2*rtUp/omega)*std::sin(omega*t+phi);
}

dcmplx fields::monochromaticField::AIfield(dcmplx t) const{
	return (2*rtUp/omega)*std::sin(omega*t+phi);
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
	if(fieldType==fieldTypes::sin2 && NCycles==1){return 0.0;}
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
	if(fieldType==fieldTypes::sin2 && NCycles==1){return 0.0;}
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

/// --- Testing Suite --- ///
TEST_CASE( "Testing SPEGrid Methods" ,"[fields]") {
	double omega=0.057, rtUp=std::sqrt(1.2);
	auto phi = GENERATE(0, 0.25*pi, 0.33*pi, pi, 0.75*pi);
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

	}
	SECTION("Testing Monochromatic Case"){
		auto N = GENERATE(2, 3, 4);
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


	}
}

