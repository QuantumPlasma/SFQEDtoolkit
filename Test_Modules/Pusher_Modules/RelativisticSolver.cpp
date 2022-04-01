#include "RelativisticSolver.h"

#include<algorithm>
#include<iomanip>

double RelativisticRungeKuttaSolver_LCFA(ofstream& file, Particle3D& part, PlaneWave3D wave, double timeStep, int iters, SFQED_Processes& proc, double (*rnd_genarator)(void)) {
	//auxiliary positions and velocities to be used in the algorithm
	//notice thar for the first position and velocity we will directly employ that of the particle
	Vector3D<double> pos2, pos3, pos4, tmpP;
	Vector3D<double> vel2, vel3, vel4, tmpV;

	//needed for intermediate cross products
	Vector3D<double> v_x_B1, v_x_B2, v_x_B3, v_x_B4;

	//the four runge-kutta terms in vectorial form
	Vector3D<double> k1, k2, k3, k4, increment, phtn_k, electron_mom;

	//Electric and magnetic fields
	Vector3D<double> B, E;

	double timeDependence,
		actualTime = 0,
		mass = part.getMass(),
		charge = part.getCharge(),
		fieldFactor = wave.geta0(),
		gamma_check,
		chi_check,
		delta,
		tau_over_tauc,
		xiprime;

	//before starting we choose to store the following timeStep depending value in a constant, as it never changes
	//we include in it the fieldFactor multiplier 
	const double constantTerm = mass > 0. ? timeStep * (charge / mass) * fieldFactor : 0.;
	const double halfTimeStep = timeStep * 0.5;
	const double oneSixth = 1. / 6.;
	const double oneThird = 1. / 3.;

	double inv_lf1, inv_lf2, inv_lf3, inv_lf4, inv_lf_final;

	double *reduced_momentum, *pushed_reduced_momentum, *mid_step_momentum, *E_actual, *B_actual;
	double LCFA_thresh, rnd, rnd2, phtn_nrg, p_minus, k_minus;

	double Lorentz_F_t_der[3],
			Lorentz_F_tt_der[3],
			Lorentz_F[3],
			inv_t = 1. / timeStep,
			inv_t_2 = inv_t * inv_t,
			*F_prev,
			*F_prev_prev,
			four_over_threepi = 4. / (3. * 3.1415926535897932384),
			threepi_over_four = 1. / four_over_threepi,
			epsilon_m = 2.22e-16,
            dOpticalD;

	bool can_continue;

	//print particle
	//cout << "At time = " << actualTime << " we have: " <<
		//part;
		//part.getPosition();
	//part.getVelocity() << '\n';

	for (int i = 0; i < iters; i++) {

		//*************FIRST TERM*******************
		//for each iteration and, in particular, for each one of the 4 terms we have to compute the
		// oscillation behaviour based on time and particle position

		//please notice the these vectors do not incorporate the a0 factor (that got absorbed inside constantTerm)
		B = wave.evolvedMagneticDirection(part.getPosition(), actualTime);
		E = wave.evolvedElectricDirection(part.getPosition(), actualTime);

		inv_lf1 = 1. / sqrt(1. + part.getVelocity() * part.getVelocity());

		// where / acts as the cross product
		v_x_B1 = (part.getVelocity() / B) * (inv_lf1);

		//compute the first runge-kutta contribution k1
		k1 = (E + v_x_B1) * constantTerm;

		//*************SECOND TERM******************
		//compute the particle velocity to be used in the second contribution term as v + (k1/2) 
		vel2 = part.getVelocity() + k1 * 0.5;

		inv_lf2 = 1. / sqrt(1. + vel2 * vel2);

		//prepare the new particle position entering the second term (remember that this process must be performed for every term)
		pos2 = part.getPosition() + vel2 * (halfTimeStep * inv_lf2);

		B = wave.evolvedMagneticDirection(pos2, actualTime + halfTimeStep);
		E = wave.evolvedElectricDirection(pos2, actualTime + halfTimeStep);

		v_x_B2 = (vel2 / B) * (inv_lf2);

		//compute the second runge-kutta contribution k2
		k2 = (E + v_x_B2) * constantTerm;

		//*************THIRD TERM*******************
		//compute the particle velocity for the third contribution term as v + (k2/2) 
		vel3 = part.getVelocity() + k2 * 0.5;

		inv_lf3 = 1. / sqrt(1. + vel3 * vel3);

		//particle position for the third term
		//(this is always based on the initial position except that we now have to take in consideration the third velocity and a half timestep)
		pos3 = part.getPosition() + vel3 * (halfTimeStep * inv_lf3);

		B = wave.evolvedMagneticDirection(pos3, actualTime + halfTimeStep);
		E = wave.evolvedElectricDirection(pos3, actualTime + halfTimeStep);

		v_x_B3 = (vel3 / B)* (inv_lf3);

		//compute the third runge-kutta contribution k3
		k3 = (E + v_x_B3) * constantTerm;

		//*************FOURTH TERM******************
		//fourth particle velocity as v + (k3) 
		vel4 = part.getVelocity() + k3;
		
		inv_lf4 = 1. / sqrt(1. + vel4 * vel4);

		//particle position for the fourth term (in this last case we proceed with a full timeStep)
		pos4 = part.getPosition() + vel4 * (timeStep * inv_lf4);

		B = wave.evolvedMagneticDirection(pos4, actualTime + timeStep);
		E = wave.evolvedElectricDirection(pos4, actualTime + timeStep);

		v_x_B4 = (vel4 / B) * (inv_lf4);

		//compute the fourth runge-kutta contribution k4
		k4 = (E + v_x_B4) * constantTerm;

		//update particle with new velocity
		increment = (k1 * oneSixth) + (k2 * oneThird) + (k3 * oneThird) + (k4 * oneSixth);
		tmpV = part.getVelocity() + increment;

		//****************************************************
		//*************PHOTON EMISSION BLCFA******************
		//****************************************************

		//set the old reduced momentum
		reduced_momentum = part.getVelocity();
		//set the lorentz force pushed momentum
		pushed_reduced_momentum = tmpV;
		//set the average of the previous momenta (in both Vector3D and pointer)
		electron_mom = (tmpV + part.getVelocity()) * 0.5;
		mid_step_momentum = electron_mom;
		
		//set the proper em fields
		E_actual = E * fieldFactor;
		B_actual = B * fieldFactor;
		gamma_check = std::sqrt(1. + electron_mom*electron_mom);
		part.gamma_alias() = gamma_check;
		chi_check = proc.compute_quantum_param(gamma_check, mid_step_momentum, E_actual, B_actual);
		part.chi_alias() = chi_check;
		delete[] E_actual;
		delete[] B_actual;
		
		
        //compute the increase to the optical depth for photon emission
        dOpticalD = proc.SFQED_PHTN_emission_rate(part.gamma_alias(), part.chi_alias()) * timeStep;

        ////////////////////////////////////////////////////
		// UNCOMMENT THIS TO USE LOCAL PROBABILITY METHOD //
		////////////////////////////////////////////////////
		rnd = rnd_genarator();
		if(dOpticalD >= rnd){
		////////////////////////////////////////////////////

		////////////////////////////////////////////////
		// UNCOMMENT THIS TO USE OPTICAL DEPTH METHOD //
		////////////////////////////////////////////////
		// //increase the optical depth
		// part.increase_optical_depth(dOpticalD);
		// //check whether we reached the optical depth reference
		// //if this is the case
		// if(part.check_optical_depth_for_emission()){
		////////////////////////////////////////////////

                
            rnd = rnd_genarator();
            phtn_nrg = proc.SFQED_LCFA_emitted_photon_energy(gamma_check, chi_check, rnd);
			//correct electron momentum after emission
            tmpV -= electron_mom * (phtn_nrg / (gamma_check * mass));

			////////////////////////////////////////////////
			// UNCOMMENT THIS TO USE OPTICAL DEPTH METHOD //
			////////////////////////////////////////////////
            // //reset the optical depth once the emission occured
            // part.reset_optical_depth();
            // rnd = rnd_genarator();
            // part.init_optical_depth(-log(1. - rnd));
			////////////////////////////////////////////////

            file << phtn_nrg << "\n";
				
        }


		//delete tmp pointers
		delete[] reduced_momentum;
		delete[] pushed_reduced_momentum;
		delete[] mid_step_momentum;


		//*************UPDATE POSITION AND MOMENTUM******************
		//update old velocity
		part.setOldVelocity(part.getVelocity());
		part.setVelocity(tmpV);

		inv_lf_final = 1. / sqrt(1. + tmpV * tmpV);

		//PUSH particle position
		tmpP = part.getPosition() + part.getVelocity() * (timeStep * inv_lf_final);
		
		//set particle position
		part.setPosition(tmpP);

		//increment time
		actualTime += timeStep;

	}

	//useful in testing mode: return the last LCFA threshold
	return LCFA_thresh;

}


double RelativisticRungeKuttaSolver_BLCFA(ofstream& file, Particle3D& part, PlaneWave3D wave, double timeStep, int iters, SFQED_Processes& proc, double (*rnd_genarator)(void)) {
	//auxiliary positions and velocities to be used in the algorithm
	//notice thar for the first position and velocity we will directly employ that of the particle
	Vector3D<double> pos2, pos3, pos4, tmpP;
	Vector3D<double> vel2, vel3, vel4, tmpV;

	//needed for intermediate cross products
	Vector3D<double> v_x_B1, v_x_B2, v_x_B3, v_x_B4;

	//the four runge-kutta terms in vectorial form
	Vector3D<double> k1, k2, k3, k4, increment, phtn_k, electron_mom;

	//Electric and magnetic fields
	Vector3D<double> B, E;

	double timeDependence,
		actualTime = 0,
		mass = part.getMass(),
		charge = part.getCharge(),
		fieldFactor = wave.geta0(),
		gamma_check,
		chi_check,
		delta,
		tau_over_tauc,
		xiprime;

	//before starting we choose to store the following timeStep depending value in a constant, as it never changes
	//we include in it the fieldFactor multiplier 
	const double constantTerm = mass > 0. ? timeStep * (charge / mass) * fieldFactor : 0.;
	const double halfTimeStep = timeStep * 0.5;
	const double oneSixth = 1. / 6.;
	const double oneThird = 1. / 3.;

	double inv_lf1, inv_lf2, inv_lf3, inv_lf4, inv_lf_final;

	double *reduced_momentum, *pushed_reduced_momentum, *mid_step_momentum, *E_actual, *B_actual;
	double LCFA_thresh, rnd, rnd2, phtn_nrg, p_minus, k_minus;

	double Lorentz_F_t_der[3],
			Lorentz_F_tt_der[3],
			Lorentz_F[3],
			inv_t = 1. / timeStep,
			inv_t_2 = inv_t * inv_t,
			*F_prev,
			*F_prev_prev,
			four_over_threepi = 4. / (3. * 3.1415926535897932384),
			threepi_over_four = 1. / four_over_threepi,
			epsilon_m = 2.22e-16,
            dOpticalD;

	bool can_continue;

	//print particle
	//cout << "At time = " << actualTime << " we have: " <<
		//part;
		//part.getPosition();
	//part.getVelocity() << '\n';

	for (int i = 0; i < iters; i++) {

		//*************FIRST TERM*******************
		//for each iteration and, in particular, for each one of the 4 terms we have to compute the oscillation behaviour based on time and particle position
		//timeDependence = wave.computeWaveOscillation(part.getPosition(), actualTime);

		//please notice the these vectors do not incorporate the a0 factor (that got absorbed inside constantTerm)
		B = wave.evolvedMagneticDirection(part.getPosition(), actualTime);
		E = wave.evolvedElectricDirection(part.getPosition(), actualTime);

		inv_lf1 = 1. / sqrt(1. + part.getVelocity() * part.getVelocity());

		// where / acts as the cross product
		v_x_B1 = (part.getVelocity() / B) * (inv_lf1);

		//compute the first runge-kutta contribution k1
		k1 = (E + v_x_B1) * constantTerm;

		//*************SECOND TERM******************
		//compute the particle velocity to be used in the second contribution term as v + (k1/2) 
		vel2 = part.getVelocity() + k1 * 0.5;

		inv_lf2 = 1. / sqrt(1. + vel2 * vel2);

		//prepare the new particle position entering the second term (remember that this process must be performed for every term)
		pos2 = part.getPosition() + vel2 * (halfTimeStep * inv_lf2);

		B = wave.evolvedMagneticDirection(pos2, actualTime + halfTimeStep);
		E = wave.evolvedElectricDirection(pos2, actualTime + halfTimeStep);

		v_x_B2 = (vel2 / B) * (inv_lf2);

		//compute the second runge-kutta contribution k2
		k2 = (E + v_x_B2) * constantTerm;

		//*************THIRD TERM*******************
		//compute the particle velocity for the third contribution term as v + (k2/2) 
		vel3 = part.getVelocity() + k2 * 0.5;

		inv_lf3 = 1. / sqrt(1. + vel3 * vel3);

		//particle position for the third term
		//(this is always based on the initial position except that we now have to take in consideration the third velocity and a half timestep)
		pos3 = part.getPosition() + vel3 * (halfTimeStep * inv_lf3);

		B = wave.evolvedMagneticDirection(pos3, actualTime + halfTimeStep);
		E = wave.evolvedElectricDirection(pos3, actualTime + halfTimeStep);

		v_x_B3 = (vel3 / B)* (inv_lf3);

		//compute the third runge-kutta contribution k3
		k3 = (E + v_x_B3) * constantTerm;

		//*************FOURTH TERM******************
		//fourth particle velocity as v + (k3) 
		vel4 = part.getVelocity() + k3;
		
		inv_lf4 = 1. / sqrt(1. + vel4 * vel4);

		//particle position for the fourth term (in this last case we proceed with a full timeStep)
		pos4 = part.getPosition() + vel4 * (timeStep * inv_lf4);

		B = wave.evolvedMagneticDirection(pos4, actualTime + timeStep);
		E = wave.evolvedElectricDirection(pos4, actualTime + timeStep);

		v_x_B4 = (vel4 / B) * (inv_lf4);

		//compute the fourth runge-kutta contribution k4
		k4 = (E + v_x_B4) * constantTerm;

		//update particle with new velocity
		increment = (k1 * oneSixth) + (k2 * oneThird) + (k3 * oneThird) + (k4 * oneSixth);
		tmpV = part.getVelocity() + increment;

		//****************************************************
		//*************PHOTON EMISSION BLCFA******************
		//****************************************************

		//set the old reduced momentum
		reduced_momentum = part.getVelocity();
		//set the lorentz force pushed momentum
		pushed_reduced_momentum = tmpV;
		//set the average of the previous momenta (in both Vector3D and pointer)
		electron_mom = (tmpV + part.getVelocity()) * 0.5;
		mid_step_momentum = electron_mom;
		
		
		can_continue = part.SFQED_BLCFA_update_entities_quantities(proc, pushed_reduced_momentum, reduced_momentum,
                                                Lorentz_F_t_der, Lorentz_F_tt_der,
                                                part.gamma_alias(), part.chi_alias());										


		if(can_continue){
            //if we enter this block it means that the particle isn't new
            //and an emission could actually happen

            //compute the increase to the optical depth for photon emission
            dOpticalD = proc.SFQED_PHTN_emission_rate(part.gamma_alias(), part.chi_alias()) * timeStep;

			//////////////////////////////////////////////////////////////////
			// FOR THE BLCFA ONLY THE LOCAL PROBABILITY METHOD IS AVAILABLE //
			//////////////////////////////////////////////////////////////////
            rnd = rnd_genarator();
            if(dOpticalD >= rnd){
			//////////////////////////////////////////////////////////////////
                

                LCFA_thresh = part.SFQED_BLCFA_find_energy_threshold(proc, Lorentz_F_t_der, Lorentz_F_tt_der, part.gamma_alias(), part.chi_alias());
                
			
                //emission section 
                // if(((proc.SFQED_PHTN_emission_rate(gamma_check, chi_check) * timeStep) >= rnd) && LCFA_thresh >= 0.){
                if(LCFA_thresh >= 0.){

                    //recreate random vars
                    rnd = rnd_genarator();
                    rnd2 = rnd_genarator();
                    // file << ", rnd_2 = " << rnd << ", rnd_3 = " << rnd2 << std::flush;
                    
                    phtn_nrg = proc.SFQED_BLCFA_emitted_photon_energy(LCFA_thresh, part.gamma_alias(), part.chi_alias(), rnd, rnd2);

                    if(phtn_nrg > 0.){

                        //write energy on file
                        file << phtn_nrg << "\n"; //<< ' ' << LCFA_thresh << ' ' << part.gamma_alias() << ' ' << part.chi_alias() << ' ' << i << "\n";

                        //adjust the reduced moment of the emitting particle
                        tmpV -= electron_mom * (phtn_nrg/(part.gamma_alias()*mass));

                    }
                    
                }
                
                //////////////////////////////////////////
            }

		}

		// file << "\n" << std::flush;

		//delete tmp pointers
		delete[] reduced_momentum;
		delete[] pushed_reduced_momentum;
		delete[] mid_step_momentum;


		//*************UPDATE POSITION AND MOMENTUM******************
		//update old velocity
		part.setOldVelocity(part.getVelocity());
		part.setVelocity(tmpV);

		inv_lf_final = 1. / sqrt(1. + tmpV * tmpV);

		//PUSH particle position
		tmpP = part.getPosition() + part.getVelocity() * (timeStep * inv_lf_final);
		
		//set particle position
		part.setPosition(tmpP);

		//increment time
		actualTime += timeStep;
	}

	//useful in testing mode: return the last LCFA threshold
	return LCFA_thresh;

}


double RelativisticRungeKuttaSolverPhotons(ofstream& file, Particle3D& part, PlaneWave3D wave, double timeStep, int iters, SFQED_Processes& proc, double (*rnd_genarator)(void)) {
	//auxiliary positions and velocities to be used in the algorithm
	//notice thar for the first position and velocity we will directly employ that of the particle
	Vector3D<double> pos2, pos3, pos4, tmpP;
	Vector3D<double> vel2, vel3, vel4, tmpV;

	//needed for intermediate cross products
	Vector3D<double> v_x_B1, v_x_B2, v_x_B3, v_x_B4;

	//the four runge-kutta terms in vectorial form
	Vector3D<double> k1, k2, k3, k4, increment, phtn_k, electron_mom;

	//Electric and magnetic fields
	Vector3D<double> B, E;

	double timeDependence,
		actualTime = 0,
		mass = part.getMass(),
		charge = part.getCharge(),
		fieldFactor = wave.geta0(),
		gamma_check,
		chi_check,
		delta,
		tau_over_tauc,
		xiprime;

	//before starting we choose to store the following timeStep depending value in a constant, as it never changes
	//we include in it the fieldFactor multiplier 
	const double constantTerm = 0.;
	const double halfTimeStep = timeStep * 0.5;
	const double oneSixth = 1. / 6.;
	const double oneThird = 1. / 3.;

	double inv_lf1, inv_lf2, inv_lf3, inv_lf4, inv_lf_final;

	double *reduced_momentum, *pushed_reduced_momentum, *mid_step_momentum, *E_actual, *B_actual;
	double LCFA_thresh, rnd, rnd2, phtn_nrg, p_minus, k_minus;

	double Lorentz_F_t_der[3],
			Lorentz_F_tt_der[3],
			Lorentz_F[3],
			inv_t = 1. / timeStep,
			inv_t_2 = inv_t * inv_t,
			*F_prev,
			*F_prev_prev,
			four_over_threepi = 4. / (3. * 3.1415926535897932384),
			threepi_over_four = 1. / four_over_threepi,
			epsilon_m = 2.22e-16,
            dOpticalD;

	bool can_continue;

	for (int i = 0; i < iters; i++) {

		//*************half timestep TERM******************
		//compute the particle velocity to be used in the second contribution term as v + (k1/2) 
		vel2 = part.getVelocity();

		inv_lf2 = 1. / sqrt(vel2 * vel2);

		//prepare the new particle position entering the second term (remember that this process must be performed for every term)
		pos2 = part.getPosition() + vel2 * (halfTimeStep * inv_lf2);

		B = wave.evolvedMagneticDirection(pos2, actualTime + halfTimeStep);
		E = wave.evolvedElectricDirection(pos2, actualTime + halfTimeStep);
		
		//update particle with new velocity
		tmpV = part.getVelocity();

		//*************************************************
		//*************PAIR CREATION LCFA******************
		//*************************************************

		//set the old reduced momentum
		reduced_momentum = part.getVelocity();
		//set the lorentz force pushed momentum
		pushed_reduced_momentum = tmpV;
		//set the average of the previous momenta (in both Vector3D and pointer)
		electron_mom = (tmpV + part.getVelocity()) * 0.5;
		mid_step_momentum = electron_mom;
											
        //Pair emission
		E_actual = E * fieldFactor;
		B_actual = B * fieldFactor;
		gamma_check = std::sqrt(electron_mom*electron_mom);
		part.gamma_alias() = gamma_check;
		chi_check = proc.compute_quantum_param(gamma_check, mid_step_momentum, E_actual, B_actual);
		part.chi_alias() = chi_check;
		delete[] E_actual;
		delete[] B_actual;

        //compute the increase to the optical depth for photon emission
        dOpticalD = proc.SFQED_PAIR_creation_rate(part.gamma_alias(), part.chi_alias()) * timeStep;

        ////////////////////////////////////////////////////
		// UNCOMMENT THIS TO USE LOCAL PROBABILITY METHOD //
		////////////////////////////////////////////////////
        rnd = rnd_genarator();
        if(dOpticalD >= rnd){
		////////////////////////////////////////////////////

		////////////////////////////////////////////////
		// UNCOMMENT THIS TO USE OPTICAL DEPTH METHOD //
		////////////////////////////////////////////////
		// //increase the optical depth
        // part.increase_optical_depth(dOpticalD);
		// //check whether we reached the optical depth reference
        // //if this is the case
        // if(part.check_optical_depth_for_emission()){
		////////////////////////////////////////////////

            
			rnd = rnd_genarator();
            phtn_nrg = proc.SFQED_PAIR_created_electron_energy(gamma_check, chi_check, rnd);
			// 	the photon correction is not needed in the 
            // tmpV -= electron_mom * (phtn_nrg / (gamma_check * mass));

			////////////////////////////////////////////////
			// UNCOMMENT THIS TO USE OPTICAL DEPTH METHOD //
			////////////////////////////////////////////////
            //reset the optical depth once the emission occured
            // part.reset_optical_depth();
            // rnd = rnd_genarator();
            // part.init_optical_depth(-log(1. - rnd));
			////////////////////////////////////////////////

            file << phtn_nrg << "\n";

			//once the pair is created abort the simulation for that photon
			break;
                
        }

		// file << "\n" << std::flush;

		//delete tmp pointers
		delete[] reduced_momentum;
		delete[] pushed_reduced_momentum;
		delete[] mid_step_momentum;


		//*************UPDATE POSITION AND MOMENTUM******************
		//update old velocity
		part.setOldVelocity(part.getVelocity());
		part.setVelocity(tmpV);

		inv_lf_final = 1. / sqrt(tmpV * tmpV);

		//PUSH particle position
		tmpP = part.getPosition() + part.getVelocity() * (timeStep * inv_lf_final);
		
		//set particle position
		part.setPosition(tmpP);

		//increment time
		actualTime += timeStep;

	}

	// file << "\n" << std::flush;

	//useful in testing mode: return the last LCFA threshold
	return LCFA_thresh;

}