! This program shows how to call the main routines provided by SFQEDtoolkit.
! The code at hand mimics what happens inside a generic particle pusher after
! the fields acting on a particle and its new momentum have been established.
! These quantities, which will not be specified in the following, can easily
! be used to determine the particle's quantum parameter chi and gamma, that
! represent our starting point. Therefore before actually moving the particle,
! electron/photon, the user could resort to this suite of functiuons to
! ascertain whether the electron/photon emits/decays /into a synchrophoton/pair.
! Of course, in order to implement these stochastic mechanism we will need some 
! random numbers, which in this case will be fixed.


PROGRAM test

   USE SFQEDtoolkit_Interface
   use iso_c_binding

   IMPLICIT NONE

   REAl(numreal64), parameter :: PI = 3.141592653589793238462643383279503_numreal64
   REAl(numreal64), parameter :: c = 299792458._numreal64 ![m/s]
   REAl(numreal64), parameter :: electron_mass = 9.10938356e-31_numreal64 ! (*kg*)
   REAl(numreal64), parameter :: electron_charge = 1.60217662e-19_numreal64 !(*coulombs*)
   REAl(numreal64), parameter :: coulombConst = 8.987551792313e9_numreal64
   REAl(numreal64), parameter :: HBar = 1.0545718e-34_numreal64 ! (*m^2kg/s*)
   REAl(numreal64), parameter :: eps_0 = 8.85418782e-12_numreal64;
   REAl(numreal64), parameter :: Schwinger_E_field = (electron_mass*electron_mass*c*c*c)/(HBar*electron_charge);

   ! this is the 'random" number we should use to ascertain whether
   ! a SFQED event occurs or not. It has been fixed to zero in such a way
   ! that the emission is forced. However, this should be a proper random number
   ! comprised between 0 and 1
   REAL(numreal64), parameter :: rnd_zero = 0._numreal64

   ! wave length associated 
   REAl(numreal64) :: wave_length ![m]
   ! simulation reference length
   REAl(numreal64) :: reference_length
   ! reference angular frequency
   REAl(numreal64) :: w_r ![1/s]
   ! time step of the simulation
	REAl(numreal64) :: time_step
   ! PIC electric field normalization
   REAl(numreal64) :: Enorm


   ! LCFA variables
   REAl(numreal64) :: gamma, chi, rnd
   LOGICAL :: toolkit_init
   ! photon emission
   REAl(numreal64) :: emission_rate, phtn_energy
   !pair production
   REAl(numreal64) :: production_rate, electron_energy
   ! values used for theoretical inspections
   REAl(numreal64) :: expected_chi
   
   ! LCFA fields and momenta useful
   REAL(numreal64), DIMENSION(1:3) :: BB, EE
   REAL(numreal64), DIMENSION(1:2,1:3) :: parts_momenta
   REAL(numreal64), DIMENSION(1:2,1:3) :: product_momenta
   ! index used insde while cycling over particles
   INTEGER :: part_num
   ! tmp helper variables
   REAL(numreal64) :: tmp_parent, tmp_child

   ! additional variables for BLCFA
   TYPE(C_PTR), DIMENSION(2) :: blcfa_entity
   REAL(numreal64) :: rnd2
   REAL(numreal64), DIMENSION(1:2,1:3) :: momentum
   REAL(numreal64), DIMENSION(1:2,1:3) :: pushed_momentum
   REAL(numreal64), DIMENSION(1:2,1:3) :: Lorentz_F_t_der
   REAL(numreal64), DIMENSION(1:2,1:3) :: Lorentz_F_tt_der
   REAL(numreal64), DIMENSION(1:2) :: delta_blcfa
   REAL(numreal64) :: LCFA_threshold
   LOGICAL :: keep_going_blcfa

   ! CHARACTER(:), ALLOCATABLE :: log_message

   print *, 'Starting Fortran TEST of SFQEDtoolkit!'

   !///////////////////
   !// SET VARIABLES //
   !///////////////////

   ! set main simualtions quantities
   time_step = 0.1_numreal64
   rnd = 2.1432683788011112e-01_numreal64
   rnd2 = 4.6667110607838297e-01_numreal64
   wave_length = 0.8e-6_numreal64
   reference_length = wave_length / (2. * PI)
   w_r = c / reference_length
   Enorm = electron_mass*w_r*c/electron_charge

   !prepare e.m. fields
   EE(1) = 0._numreal64
   EE(2) = 0._numreal64
   EE(3) = 0._numreal64

   BB(1) = 0._numreal64
   BB(2) = 270._numreal64
   BB(3) = 0._numreal64

   ! parts momenta
   parts_momenta(1,1) = -1.2137280955271491e+01_numreal64
   parts_momenta(1,2) = 0.0_numreal64
   parts_momenta(1,3) = -1.9570511118175094e+04_numreal64

   parts_momenta(2,1) = 5.9499124858705343e+00_numreal64
   parts_momenta(2,2) = 0.0_numreal64
   parts_momenta(2,3) = -1.0000512547779021e+04_numreal64


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! SFQEDtoolkit INITIALIZATION !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! init SFQEDtoolkit by providing the reference length
   ! and the timestep of the simulation
   toolkit_init = SFQED_INIT_ALL_ref_len(reference_length, time_step)
   ! toolkit_init = SFQED_INIT_ALL_ref_freq(w_r, time_step)
   PRINT *, 'SFQEDtoolkit initialization = ', toolkit_init
   ! PRINT *, log_message

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! LCFA synchrophoton emission !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   print *, ''
   PRINT *, 'LCFA SECTION:'

   ! please notice that in both the photon emission and pair production sections
   ! we are cycling over the exact same array of momenta (i.e. the parts_momenta
   ! variable defined above). Of course you are free to redefine its components
   ! as you please. The same holds true for the magnetic and electric fields.
   ! If you do so, mind that the expected_chi variablle is generally no longer
   ! reliable, as it was meant to represent the theoretical value of the quantum
   ! parameter chi when the particles' momenta are perpendicular to the magnetic
   ! field and the electric one is turned off
   do part_num = 1,2

      !compute chi and gamma
      gamma = SQRT(1. + scalar_prod(parts_momenta(part_num, 1:3), parts_momenta(part_num, 1:3)))
      chi = compute_chi_with_vectors(gamma, parts_momenta(part_num, 1:3), EE, BB)
      
      !theoretical chi expected when the electric field is zero
      !and the magnetic field is perpendicular to the momenta
      !(this must be compared with the chi obtained from the routine) 
      expected_chi = gamma * Enorm * SQRT(scalar_prod(BB,BB)) / Schwinger_E_field;
      
      PRINT *, '#particle=', part_num, ', gamma=', gamma, ', chi=', chi, ' (', expected_chi, ')'

      !actual SFQED routines
      emission_rate = SFQED_INV_COMPTON_rate(gamma, chi)
      PRINT *, 'Emission rate = ', emission_rate

      IF (emission_rate * time_step >= rnd_zero) THEN
         phtn_energy = SFQED_LCFA_INV_COMPTON_PHOTON_energy(gamma, chi, rnd) 
         PRINT *, 'Emitted photon energy = ', phtn_energy

         !compute new particle's momentum
         call SFQED_build_collinear_momentum(phtn_energy, parts_momenta(part_num,1:3), product_momenta(part_num,1:3))

         !output momenta for checking
         tmp_parent = sqrt(scalar_prod(parts_momenta(part_num,1:3),parts_momenta(part_num,1:3)))
         tmp_child = sqrt(scalar_prod(product_momenta(part_num,1:3),product_momenta(part_num,1:3)))
         print *, 'Parent momentum = ', parts_momenta(part_num,1), parts_momenta(part_num,2), parts_momenta(part_num,3)
         print *, 'Parent direction = ', parts_momenta(part_num,1) / tmp_parent, parts_momenta(part_num,2) / tmp_parent, &
                                          parts_momenta(part_num,3) / tmp_parent
         print *, 'Child momentum = ', product_momenta(part_num,1), product_momenta(part_num,2), product_momenta(part_num,3)
         print *, 'Child direction = ', product_momenta(part_num,1) / tmp_child, product_momenta(part_num,2) / tmp_child, &
                                          product_momenta(part_num,3) / tmp_child
      END IF

      print *, ''

   end do 

   !!!!!!!!!!!!!!!!!!!
   ! PAIR production !
   !!!!!!!!!!!!!!!!!!!
   
   print *, ''
   print *, 'PAIRS SECTION:'

   do part_num = 1,2
      !compute chi and gamma
      gamma = SQRT(1. + scalar_prod(parts_momenta(part_num, 1:3), parts_momenta(part_num, 1:3)))
      chi = compute_chi_with_vectors(gamma, parts_momenta(part_num, 1:3), EE, BB)

      !theoretical chi expected when the electric field is zero
      !and the magnetic field is perpendicular to the momenta
      !(this must be compared with the chi obtained from the routine) 
      expected_chi = gamma * Enorm * SQRT(scalar_prod(BB,BB)) / Schwinger_E_field;

      PRINT *, '#particle=', part_num, ', gamma=', gamma, ', chi=', chi, ' (', expected_chi, ')'

      !actual SFQED routines
      production_rate = SFQED_BREIT_WHEELER_rate(gamma, chi)
      PRINT *, 'production rate = ', production_rate

      IF (production_rate * time_step >= rnd_zero) THEN
         electron_energy = SFQED_BREIT_WHEELER_ELECTRON_energy(gamma, chi, rnd)
         PRINT *, 'Created electron energy = ', electron_energy

         !compute new particle's momentum
         call SFQED_build_collinear_momentum(electron_energy, parts_momenta(part_num,1:3), product_momenta(part_num,1:3))

         !output momenta for checking
         tmp_parent = sqrt(scalar_prod(parts_momenta(part_num,1:3),parts_momenta(part_num,1:3)))
         tmp_child = sqrt(scalar_prod(product_momenta(part_num,1:3),product_momenta(part_num,1:3)))
         print *, 'Parent momentum = ', parts_momenta(part_num,1), parts_momenta(part_num,2), parts_momenta(part_num,3)
         print *, 'Parent direction = ', parts_momenta(part_num,1) / tmp_parent, parts_momenta(part_num,2) / tmp_parent, &
                                          parts_momenta(part_num,3) / tmp_parent
         print *, 'Child momentum = ', product_momenta(part_num,1), product_momenta(part_num,2), product_momenta(part_num,3)
         print *, 'Child direction = ', product_momenta(part_num,1) / tmp_child, product_momenta(part_num,2) / tmp_child, &
                                          product_momenta(part_num,3) / tmp_child
      END IF

      print *, ''

   end do

   !!!!!!!!!!!!!!!!!!!!!!!
   ! BLCFA synchrophoton !
   !!!!!!!!!!!!!!!!!!!!!!!

   ! set BLCFA quantities
   blcfa_entity(1) = SFQED_CREATE_BLCFA_OBJECT()
   blcfa_entity(2) = SFQED_CREATE_BLCFA_OBJECT()

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!! PRELIMINARY BLCFA ACTIONS !!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! when using the blcfa approach we need some preliminary 
   ! actions to properly set the lorentz force derivatives
   ! and chi and gamma. We will perform the same action 3
   ! times, since we are interested in getting the second time
   ! derivatives of the force

   !///////// STEP 1 /////////

   !old momentum
   momentum(1,1) = 2.0923765116618465e+00_numreal64
   momentum(1,2) = 0.0000000000000000e+00_numreal64
   momentum(1,3) = -1.9570512944081507e+04_numreal64
   !new pushed momentum
   pushed_momentum(1,1) = 5.9499124858705343e+00_numreal64
   pushed_momentum(1,2) =  0.0000000000000000e+00_numreal64
   pushed_momentum(1,3) = -1.9570512547779021e+04_numreal64
   !update blcfa obj
   keep_going_blcfa = SFQED_BLCFA_OBJECT_update(blcfa_entity(1), pushed_momentum(1,1:3), momentum(1,1:3), &
                              delta_blcfa(1), gamma, chi)

   !old momentum
   momentum(2,1) = -1.7718145323355337e+01_numreal64
   momentum(2,2) = 0.0000000000000000e+00_numreal64
   momentum(2,3) = -1.9570508989729224e+04_numreal64
   !new pushed momentum
   pushed_momentum(2,1) = -1.6542717780968854e+01_numreal64
   pushed_momentum(2,2) =  0.0000000000000000e+00_numreal64
   pushed_momentum(2,3) = -1.9570509504165970e+04_numreal64
   !update blcfa obj
   keep_going_blcfa = SFQED_BLCFA_OBJECT_update(blcfa_entity(2), pushed_momentum(2,1:3), momentum(2,1:3), &
                                                   delta_blcfa(2), gamma, chi);

   !///////// STEP 2 /////////
   
   !old momentum                           
   momentum(1,1) = 5.9499124858705343e+00_numreal64
   momentum(1,2) =  0.0000000000000000e+00_numreal64
   momentum(1,3) = -1.9570512547779021e+04_numreal64
   !new pushed momentum
   pushed_momentum(1,1) = 9.6034048452081251e+00_numreal64
   pushed_momentum(1,2) = 0.0000000000000000e+00_numreal64
   pushed_momentum(1,3) = -1.9570511821892545e+04_numreal64
   !update blcfa obj
   keep_going_blcfa = SFQED_BLCFA_OBJECT_update(blcfa_entity(1), pushed_momentum(1,1:3), momentum(1,1:3), &
                              delta_blcfa(1), gamma, chi)

   !old momentum
   momentum(2,1) = -1.6542717780968854e+01_numreal64
   momentum(2,2) = 0.0000000000000000e+00_numreal64
   momentum(2,3) = -1.9570509504165970e+04_numreal64
   !new pushed momentum
   pushed_momentum(2,1) = -1.4660387531605418e+01_numreal64
   pushed_momentum(2,2) =  0.0000000000000000e+00_numreal64
   pushed_momentum(2,3) = -1.9570510254460005e+04_numreal64
   !update blcfa obj
   keep_going_blcfa = SFQED_BLCFA_OBJECT_update(blcfa_entity(2), pushed_momentum(2,1:3), momentum(2,1:3), &
                                                   delta_blcfa(2), gamma, chi);

   !///////// STEP 3 /////////
   
   !old momentum 
   momentum(1,1) = 9.6034048452081251e+00_numreal64
   momentum(1,2) = 0.0000000000000000e+00_numreal64
   momentum(1,3) = -1.9570511821892545e+04_numreal64
   !new pushed momentum
   pushed_momentum(1,1) = 1.2895899733826127e+01_numreal64
   pushed_momentum(1,2) = 0.0000000000000000e+00_numreal64
   pushed_momentum(1,3) = -1.9570510875586107e+04_numreal64

   !old momentum
   momentum(2,1) = -1.4660387531605418e+01_numreal64
   momentum(2,2) =  0.0000000000000000e+00_numreal64
   momentum(2,3) = -1.9570510254460005e+04_numreal64
   !new pushed momentum
   pushed_momentum(2,1)= -1.2137280955271491e+01_numreal64
   pushed_momentum(2,2) =  0.0000000000000000e+00_numreal64
   pushed_momentum(2,3) = -1.9570511118175094e+04_numreal64

   !////////////  ACTUAL BLCFA LOOP  ////////////
   ! at this point everything should be fine and the actual blcfa section can now begin
   print *, ''
   print *, 'BLCFA SECTION:'

   do part_num = 1,2

      ! update the blcfa objects. Remember that this update operation must be 
      ! done at each timestep of the simulation for every particle. The update
      ! requires you to know the initial momentum possessed by the given particle,
      ! together with the pushed momentum, which would result by normally applying
      ! the fields and forces on the particle. The determination of these quantities,
      ! which should be performed in this same loop has been here omitted, and replaced
      ! with the 3 preliminary operations above
      keep_going_blcfa = SFQED_BLCFA_OBJECT_update(blcfa_entity(part_num), pushed_momentum(part_num,1:3), momentum(part_num,1:3), &
                                             delta_blcfa(part_num), gamma, chi)


      print *, '#particle = ', part_num
      print *, 'chi = ', chi
      print *, 'gamma = ', gamma
      ! print *, 'first derivative = ', Lorentz_F_t_der(part_num,1), Lorentz_F_t_der(part_num,2), Lorentz_F_t_der(part_num,3)
      ! print *, 'second derivative = ', Lorentz_F_tt_der(part_num,1), Lorentz_F_tt_der(part_num,2), Lorentz_F_tt_der(part_num,3)

      emission_rate = SFQED_INV_COMPTON_rate(gamma, chi)
      PRINT *, 'Emission rate = ', emission_rate

      IF (keep_going_blcfa .AND. emission_rate * time_step >= rnd_zero) THEN

         LCFA_threshold = SFQED_BLCFA_INV_COMPTON_PHOTON_threshold(blcfa_entity(part_num), delta_blcfa(part_num), &
                                       gamma, chi)
         print *, 'LCFA threshold = ', LCFA_threshold

         phtn_energy = SFQED_BLCFA_INV_COMPTON_PHOTON_energy(LCFA_threshold, gamma, chi, rnd, rnd2)
         PRINT *, 'Emitted photon energy = ', phtn_energy

         !now you should check if phtn_energy > 0. In that case the emission occurs, otherwise
         ! nothing happens

         !compute new particle's momentum
         call SFQED_build_collinear_momentum(phtn_energy, pushed_momentum(part_num,1:3), product_momenta(part_num,1:3));

         !output momenta for checking
         tmp_parent = sqrt(scalar_prod(pushed_momentum(part_num,1:3),pushed_momentum(part_num,1:3)));
         tmp_child = sqrt(scalar_prod(product_momenta(part_num,1:3),product_momenta(part_num,1:3)));
         print *, 'Parent momentum = ', pushed_momentum(part_num,1), pushed_momentum(part_num,2), pushed_momentum(part_num,3)
         print *, 'Parent direction = ', pushed_momentum(part_num,1) / tmp_parent, pushed_momentum(part_num,2) / tmp_parent, &
                                          pushed_momentum(part_num,3) / tmp_parent
         print *, 'Child momentum = ', product_momenta(part_num,1), product_momenta(part_num,2), product_momenta(part_num,3)
         print *, 'Child direction = ', product_momenta(part_num,1) / tmp_child, product_momenta(part_num,2) / tmp_child, &
                                          product_momenta(part_num,3) / tmp_child

      END IF

      print *, ''

   end do


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! SFQEDtoolkit FINALIZATION !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call SFQED_FINALIZE_ALL()
    
   contains

      function scalar_prod(arr1, arr2) result (scalar)
         REAl(numreal64) :: scalar
         REAl(numreal64), DIMENSION(3), INTENT(IN) :: arr1
         REAl(numreal64), DIMENSION(3), INTENT(IN) :: arr2
         scalar = arr1(1) * arr2(1) + arr1(2) * arr2(2) + arr1(3) * arr2(3)
      end function 
  
  END PROGRAM test