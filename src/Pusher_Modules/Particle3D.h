#pragma once

#include <iostream>
#include <fstream>

#include "Vector3D.h"
#include "SFQED_Processes.h"
// #include "BLCFA_Object.h" //eventually change this to SFQED_Processes.h

class Particle3D : public BLCFA_Object{
private:
	Vector3D<double> position;
	//this vector represents the reduced momentum or 4-velocity
	Vector3D<double> velocity;
	//the old velocity vector is useful when testing
	Vector3D<double> velocity_old;
	//mass and charge are given in terms of those of the electron
	double mass;
	double charge;

	//useful while testing
	double chi;
	double gamma;


public:
	//constructors
	Particle3D() {
		this->position = Vector3D<double>(0.0, "position");
		this->velocity = Vector3D<double>(0.0, "velocity");
		this->velocity_old = Vector3D<double>(0.0, "velocity_old");

		mass = 1;
		charge = 1;
	}

	Particle3D(Vector3D<double> position, Vector3D<double> velocity) {
		this->position = position;
		this->velocity = velocity;
		this->velocity_old = Vector3D<double>(0.0, "velocity_old");

		mass = 1;
		charge = 1;
	}

	Particle3D(Vector3D<double> position, Vector3D<double> velocity, double m, double q) : Particle3D(position, velocity) {
		mass = m;
		charge = q;
	}

	//setters
	void setPosition(Vector3D<double> pos) {
		this->position = pos;
	}

	void setVelocity(Vector3D<double> vel) {
		this->velocity = vel;
	}

	void setOldVelocity(Vector3D<double> vel_old) {
		this->velocity_old = vel_old;
	}

	void setMass(double m) {
		this->mass = m;
	}

	void setCharge(double q) {
		this->charge = q;
	}

	//getters
	Vector3D<double> getPosition() const {
		return position;
	}

	Vector3D<double> getVelocity() const {
		return velocity;
	}

	Vector3D<double> getOldVelocity() const {
		return velocity_old;
	}
	
	double getMass() const {
		return mass;
	}
	
	double getCharge() const {
		return charge;
	}

	//When a function returns a reference, it returns an implicit
	// pointer to its return value. This way, a function can be
	//used on the left side of an assignment statement. In other 
	//words a functino returning a reference just returns an alias
	//for the returned variable (this means its adress doesn't change)
	double& chi_alias(){
		return this->chi;
	}

	double& gamma_alias(){
		return this->gamma;
	}

	//overloading ostream operator
	friend std::ostream& operator<<(std::ostream& os, const Particle3D& part) {
		return os << "Particle\n" << part.position << part.velocity << part.velocity_old
					<< "Chi = " << part.chi << ", Gamma = " << part.gamma << "\n";
	}

	//overloading ofstream operator for file output
	friend std::ostream& operator<<(std::ofstream& ofs, const Particle3D& part) {
		return ofs << part.position.getX() << ' ' << part.position.getY() << ' ' << part.position.getZ();
	}
};