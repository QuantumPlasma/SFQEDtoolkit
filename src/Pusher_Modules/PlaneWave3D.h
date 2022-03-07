#pragma once

#include <cmath>
#include "Vector3D.h"

#define pi_over_2 1.5707963267948966192

class PlaneWave3D {
private:
	//direction of the electric field
	Vector3D<double> k;
	//used to represent the vector perp to A in the plane perp to k (that is the direction of the magnetic field)
	Vector3D<double> a2;
	Vector3D<double> A;
	Vector3D<double> E;
	Vector3D<double> B;

	//polatization
	double delta;
	//1 - delta^2
	double delta_comp;

	bool constantFields;

	//variable needed to describe a gaussian time pulse
	// static const double pi_over_2 = 1.5707963267948966192;
	bool gaussian_pulse;
	double peak_pos;
	double gauss_a0;
	double deltaphi2;
	double inv_deltaphi2;

	
	//this stores the factor the fields should be multiplied by. (It is the a0 of the field) 
	//Notice it has to be set after the plane wave instantiation, as the constructor
	//provides for its initialization to 1
	double a0;

	//this double represents the wave frequency. its definition is up to the user and I plan to use its value only at the end of the algorithm
	//hence it won't be set in the constructor
	double w;

	//this stands for \lambda/2PI and is related to w by w * l = c (thus it isn't really necessary)
	double l;
	
	void evaluateEInternal() {
		//since we chose A to represent just the direction of the electric field (i.e. its versor)
		//it is no longer necessary to multiply the vector potential by the frequency to obtain the proper E field vector
		//please, notice we are assuming A has been previously normalized
		this->E = this->A; // * this->w;
		this->E.setName("E");
	}

	void evaluateBInternal() {
		//since we chose k x A to represent just the direction of the magnetic field (i.e. its versor)
		//it is no longer necessary to multiply (k x A) by the wave length to obtain the proper B field vector
		//in this case we are forced to renormalize the cross product, as k and A aren't always perpendicular
		this->B = normalized(this->k / this->A); //*l;
		this->B.setName("B");
	}

public:
	//constructors
	PlaneWave3D(Vector3D<double> k, Vector3D<double> A) {
		//we are assuming that the control on the validity of k and A has already been done (that is the module of k and A is not 0)
		this->k = normalized(k);
		this->k.setName("k versor");

		//normalize the cross prod k x A
		this->a2 = normalized(k / A);
		this->a2.setName("k x A versor");

		//before roughly setting A we have to extract its perpendicular component with respect to k: A_perp = (k cross (k cross A))/k^2
		//of course, we just normalized k so k^2 = 1, and, since we are going to use the normalized a2, we don't need to normalize the following cross prod
		this->A = a2 / k;
		this->A.setName("A versor");
		
		this->a0 = 1;

		evaluateEInternal();
		evaluateBInternal();

		this->gaussian_pulse = false;
		this->constantFields = false;
	}
	
	PlaneWave3D() : PlaneWave3D(Vector3D<double>(0, 0, 1), Vector3D<double>(0, 1, 0)) {
		
	}

	//actually we need this constructor to build a combination of constant E and B fields
	PlaneWave3D(Vector3D<double> E, Vector3D<double> B, bool boolean) {
		this->E = E;
		this->E.setName("E");

		this->B = B;
		this->B.setName("B");

		this->A = Vector3D<double>(0, 0, 0);
		this->A.setName("A");

		this->k = Vector3D<double>(0, 0, 0);
		this->k.setName("k");

		this->a0 = 1;

		this->constantFields = boolean;

	}

	//gaussian enveloped plane wave constructor
	PlaneWave3D(Vector3D<double> k_prop_dir, Vector3D<double> A_dir, double a0, double fwhm, double peak_pos) :
		PlaneWave3D(k_prop_dir, A_dir){

			this->gaussian_pulse = true;
			
			this->gauss_a0 = a0;
			this->peak_pos = peak_pos;

			this->deltaphi2 = fwhm * fwhm / log(2.) * 0.5;

			this->inv_deltaphi2 = 1. / this->deltaphi2;
	}

	//setters
	void setk(Vector3D<double> vec) {
		this->k = vec;
	}
	void setA(Vector3D<double> vec) {
		this->A = vec;
	}
	void setE(Vector3D<double> vec) {
		this->E = vec;
	}
	void setB(Vector3D<double> vec) {
		this->B = vec;
	}
	void seta0(double a0) {
		this->a0 = a0;
	}
	void setConstantFields(bool b) {
		this->constantFields = b;
	}
	void setDelta(double delta) {
		this->delta = delta;
		this->delta_comp = sqrt(1 - delta * delta);
	}

	void turnOffElectricF() {
		this->E = Vector3D<double>(0,0,0);
	}

	void turnOffMagneticF() {
		this->B = Vector3D<double>(0, 0, 0);
	}

	//getters
	Vector3D<double> getk() const {
		return k;
	}
	Vector3D<double> getA() const {
		return A;
	}
	Vector3D<double> geta2() const {
		return a2;
	}
	Vector3D<double> getE() const {
		return E;
	}
	Vector3D<double> getB() const {
		return B;
	}
	double geta0() const {
		return a0;
	}
	bool getConstantFields() const {
		return constantFields;
	}

	double computeWaveOscillation(Vector3D<double> position, double evaluationTime) {
		//double arg = this->k * position - evaluationTime;

		return (this->constantFields) ? 1 : -sin((this->k * position) - evaluationTime);
	}

	Vector3D<double> extractEightParticleMomentum(Vector3D<double> position, double evaluationTime, double mass, double charge) {
		double newMultiplier = this->a0;//togli il delta!!!
		double alpha_and_charge_over_mass_and_c = (newMultiplier * charge) / (mass);
		Vector3D<double> vectorPotComponent = getADirection(position, evaluationTime) * (-alpha_and_charge_over_mass_and_c);
		double component = 0.5 * alpha_and_charge_over_mass_and_c * alpha_and_charge_over_mass_and_c;
		Vector3D<double> waveVecComponent = k * (component);//togli il meno!!!
		return vectorPotComponent + waveVecComponent;
	}

	Vector3D<double> extractEightParticleMomentumLinear(Vector3D<double> position, double evaluationTime, double mass, double charge) {
		double alpha_and_charge_over_mass_and_c = (this->a0 * charge) / (mass);
		double squaredF = alpha_and_charge_over_mass_and_c * alpha_and_charge_over_mass_and_c;
		Vector3D<double> vectorPotComponent = getADirection(position, evaluationTime) * (-alpha_and_charge_over_mass_and_c);
		double component = (0.25 * squaredF) / sqrt(1. + 0.5 * squaredF);
		Vector3D<double> waveVecComponent = k * (component);
		return vectorPotComponent + waveVecComponent;
	}

	Vector3D<double> getADirection(Vector3D<double> position, double evaluationTime) {
		double arg = this->k * position - evaluationTime;
		return (A * (delta * (cos(arg))) + a2 * ((-delta_comp) * sin(arg))) * (1.0/this->delta);//togli divisione delta
	}

	Vector3D<double> evolvedElectricDirection(Vector3D<double> position, double evaluationTime) {
		if(gaussian_pulse){
			double tmp = evaluationTime - position.getZ() + peak_pos;
			tmp = -gauss_a0 * exp(-tmp*tmp * inv_deltaphi2) * inv_deltaphi2 * 
        			(deltaphi2 * cos(tmp + pi_over_2) - 2.0 * tmp * sin(tmp + pi_over_2));

			return E * tmp;

		}
		
		if (constantFields) {
			return E;
		}

		double arg = this->k * position - evaluationTime;
		return (A * (delta * (-sin(arg))) + a2 * ((-delta_comp) * cos(arg))) * (1.0 / this->delta);//togli divisione delta
	}

	Vector3D<double> evolvedMagneticDirection(Vector3D<double> position, double evaluationTime) {

		if(gaussian_pulse){
			double tmp = evaluationTime - position.getZ() + peak_pos;
			tmp = -gauss_a0 * exp(-tmp*tmp * inv_deltaphi2) * inv_deltaphi2 * 
        			(deltaphi2 * cos(tmp + pi_over_2) - 2.0 * tmp * sin(tmp + pi_over_2));

			return B * tmp;

		}

		if (constantFields) {
			return B;
		}

		double arg = this->k * position - evaluationTime;
		return  (A * (delta_comp * cos(arg)) + a2 * (delta * (-sin(arg)))) * (1.0 / this->delta);//togli divisione delta
	}

	//normalized version***************************** start
	Vector3D<double> extractEightParticleMomentumNormalized(Vector3D<double> position, double evaluationTime, double mass, double charge) {
		double alpha_and_charge_over_mass_and_c = (this->a0 * charge) / (mass);
		Vector3D<double> vectorPotComponent = getADirectionNormalized(position, evaluationTime) * (-alpha_and_charge_over_mass_and_c);
		double component = 0.5 * alpha_and_charge_over_mass_and_c * alpha_and_charge_over_mass_and_c;
		Vector3D<double> waveVecComponent = k * (component);
		return vectorPotComponent + waveVecComponent;
	}

	Vector3D<double> extractEightParticleMomentumLinearNormalized(Vector3D<double> position, double evaluationTime, double mass, double charge) {
		double alpha_and_charge_over_mass_and_c = (this->a0 * charge) / (mass);
		double squaredF = alpha_and_charge_over_mass_and_c * alpha_and_charge_over_mass_and_c;
		Vector3D<double> vectorPotComponent = getADirectionNormalized(position, evaluationTime) * (-alpha_and_charge_over_mass_and_c);
		double component = (0.25 * squaredF) / sqrt(1. + 0.5 * squaredF);
		Vector3D<double> waveVecComponent = k * (component);
		return vectorPotComponent + waveVecComponent;
	}

	Vector3D<double> getADirectionNormalized(Vector3D<double> position, double evaluationTime) {
		double arg = this->k * position - evaluationTime;
		return normalized(A * (delta * (cos(arg))) + a2 * ((-delta_comp) * sin(arg)));
	}

	Vector3D<double> evolvedElectricDirectionNormalized(Vector3D<double> position, double evaluationTime) {
		if (constantFields) {
			return E;
		}
		double arg = this->k * position - evaluationTime;
		return normalized(A * (delta * (-sin(arg))) + a2 * ((-delta_comp) * cos(arg)));
	}

	Vector3D<double> evolvedMagneticDirectionNormalized(Vector3D<double> position, double evaluationTime) {
		if (constantFields) {
			return B;
		}
		double arg = this->k * position - evaluationTime;
		return  normalized(A * (delta_comp * cos(arg)) + a2 * (delta * (-sin(arg))));
	}
	//normalized version***************************** end

	static Vector3D<double> normalized(const Vector3D<double>& vec) {
		double inverseMod = 1/(sqrt(vec * vec));
		return vec * inverseMod;
	}

	//overloading ostream operator
	friend std::ostream& operator<<(std::ostream& os, const PlaneWave3D& wave) {
		return std::cout << "Current wave values:\nMultiplier = " << wave.geta0() << "\n" << wave.getA() << wave.geta2() << wave.getk() << wave.getE() << wave.getB()
			<< "Constant fields = " << wave.constantFields << '\n'
			<< "delta = " << wave.delta << '\n'
			<< "is gaussian = " << wave.gaussian_pulse << '\n';
	}
};