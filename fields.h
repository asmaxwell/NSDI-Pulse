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
#include <memory>

using vec4 = std::array<double,4>;
using dcmplx = std::complex<double>;

enum fieldTypes {
    monochromatic,
	sin2
};

namespace fields {
	class laserField {
	public:
		laserField();
		laserField(double omega, double rtUp, double phi, int N, fieldTypes fieldType_);
		virtual ~laserField();
		// Copy constructor
		//laserField(const laserField &LF1);

		virtual double Efield(double t) const;
		virtual dcmplx Efield(dcmplx t) const;
		virtual double Afield(double t) const;
		double A2field(double t) const;
		virtual dcmplx Afield(dcmplx t) const;

		virtual dcmplx dtAfield(dcmplx t) const;
		virtual double AIfield(double t) const=0;
		virtual dcmplx AIfield(dcmplx t) const=0;
		virtual double A2Ifield(double t) const=0;
		virtual dcmplx A2Ifield(dcmplx t) const=0;
		//double AI2field(double t) const;
		//dcmplx AI2field(dcmplx t) const;

		//get methods
		const double getOmega() const;
		const double getrtUp() const;
		const double getPhi() const;
		const int getNCycles() const;
		const fieldTypes getFieldType() const;
		virtual double F(double t) const=0;
	protected:
		virtual dcmplx F(dcmplx t) const=0;
		virtual double dF(double t) const=0;
		virtual dcmplx dF(dcmplx t) const=0;
		double omega, rtUp, phi, Up;//ang. freq, square root of up and CEP
		int NCycles;// number of cycles
		fieldTypes fieldType;

	};
	class monochromaticField : public laserField {
	public:
		using laserField::laserField;

		double AIfield(double t) const;
		dcmplx AIfield(dcmplx t) const;
		double A2Ifield(double t) const;
		dcmplx A2Ifield(dcmplx t) const;
		//double A_end(double tf) const; // not necessary for monochromatic fields if an appropiate final time chosen
	private:
		double F(double t) const {return 2*rtUp;}
		dcmplx F(dcmplx t) const {return 2*rtUp;}
		double dF(double t) const {return 0;}
		dcmplx dF(dcmplx t) const {return 0;}

	};
	class sin2 : public laserField {
	public:
		using laserField::laserField;

		double AIfield(double t) const;
		dcmplx AIfield(dcmplx t) const;
		double A2Ifield(double t) const;
		dcmplx A2Ifield(dcmplx t) const;

		double F(double t) const;
	private:
		dcmplx F(dcmplx t) const;
		double dF(double t) const;
		dcmplx dF(dcmplx t) const;

	};
};
using laserField = std::shared_ptr<fields::laserField>;

#endif /* FIELDS_H_ */
