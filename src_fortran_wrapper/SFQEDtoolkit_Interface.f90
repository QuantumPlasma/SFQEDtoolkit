MODULE SFQEDtoolkit_Interface

    USE, INTRINSIC :: ISO_C_Binding, ONLY: C_INT, C_DOUBLE, C_CHAR, C_NULL_CHAR, C_SIZE_T, C_BOOL, C_PTR, C_NULL_ptr
  
    IMPLICIT NONE

    ! type BLCFA_Object 
    !   private
    !   type(C_ptr) :: object = C_NULL_ptr
    ! end type BLCFA_Object

    INTEGER, PARAMETER :: numreal64 = KIND(1.d0)
  
    INTERFACE
      
      SUBROUTINE SFQED_INIT_PROCS() BIND(C, NAME="SFQED_INIT_PROCS")
      END SUBROUTINE SFQED_INIT_PROCS

      !!!!!!!!!!!!!!!!!!!!!!
      ! SIMULATION SETTERS !
      !!!!!!!!!!!!!!!!!!!!!!
      
      !!!!!  COMPLETE SETTERS AND FINALIZERS  !!!!!

      ! transforming a C array of char into a fortran on it is not trivial at all! Look
      ! https://stackoverflow.com/questions/41247242/c-and-fortran-interoperability-for-strings
      ! for inspiration

      FUNCTION SFQED_INIT_ALL_ref_len(ref_len, ts) RESULT(state) BIND(C, NAME="SFQED_INIT_ALL_ref_len")
        IMPORT :: C_DOUBLE, C_BOOL
        LOGICAL(C_BOOL) :: state
        REAL(C_DOUBLE), INTENT(IN) :: ref_len
        REAL(C_DOUBLE), INTENT(IN) :: ts
        ! CHARACTER(:), ALLOCATABLE, INTENT(OUT) :: log_message
      END FUNCTION SFQED_INIT_ALL_ref_len

      FUNCTION SFQED_INIT_ALL_ref_freq(ref_freq, ts) RESULT(state) BIND(C, NAME="SFQED_INIT_ALL_ref_freq")
        IMPORT :: C_DOUBLE, C_BOOL
        LOGICAL(C_BOOL) :: state
        REAL(C_DOUBLE), INTENT(IN) :: ref_freq
        REAL(C_DOUBLE), INTENT(IN) :: ts
        ! CHARACTER(:), ALLOCATABLE, INTENT(OUT) :: log_message
      END FUNCTION SFQED_INIT_ALL_ref_freq

      SUBROUTINE SFQED_INIT_ALL_debug() BIND(C, NAME="SFQED_INIT_ALL_debug")
        ! IMPORT :: C_CHAR
        ! CHARACTER(:), ALLOCATABLE, INTENT(OUT) :: log_message
      END SUBROUTINE SFQED_INIT_ALL_debug

      SUBROUTINE SFQED_FINALIZE_ALL() BIND(C, NAME="SFQED_FINALIZE_ALL")
      END SUBROUTINE SFQED_FINALIZE_ALL

      !!!!!  PARTIAL SETTERS  !!!!!
      
      SUBROUTINE SFQED_set_ref_len(ref_len) BIND(C, NAME="SFQED_set_ref_len")
        IMPORT :: C_DOUBLE
        REAL(C_DOUBLE), INTENT(IN) :: ref_len
      END SUBROUTINE SFQED_set_ref_len

      SUBROUTINE SFQED_set_ref_freq(ref_freq) BIND(C, NAME="SFQED_set_ref_freq")
        IMPORT :: C_DOUBLE
        REAL(C_DOUBLE), INTENT(IN) :: ref_freq
      END SUBROUTINE SFQED_set_ref_freq

      SUBROUTINE SFQED_set_for_debug() BIND(C, NAME="SFQED_set_for_debug")
      END SUBROUTINE SFQED_set_for_debug

      SUBROUTINE SFQED_set_sim_tstep(ts) BIND(C, NAME="SFQED_set_sim_tstep")
        IMPORT :: C_DOUBLE
        REAL(C_DOUBLE), INTENT(IN) :: ts
      END SUBROUTINE SFQED_set_sim_tstep

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! COEFFICIENT INITIALIZERS AND FINALIZERS (OBSOLETE)! 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! SUBROUTINE SFQED_init_INV_COMPTON(path) BIND(C, NAME="SFQED_init_INV_COMPTON")
      !   IMPORT :: C_CHAR
      !   CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: path
      ! END SUBROUTINE SFQED_init_INV_COMPTON

      ! SUBROUTINE SFQED_finalize_INV_COMPTON() BIND(C, NAME="SFQED_finalize_INV_COMPTON")
      ! END SUBROUTINE SFQED_finalize_INV_COMPTON

      ! SUBROUTINE SFQED_init_BREIT_WHEELER(path) BIND(C, NAME="SFQED_init_BREIT_WHEELER")
      !   IMPORT :: C_CHAR
      !   CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: path
      ! END SUBROUTINE SFQED_init_BREIT_WHEELER

      ! SUBROUTINE SFQED_finalize_BREIT_WHEELER() BIND(C, NAME="SFQED_finalize_BREIT_WHEELER")
      ! END SUBROUTINE SFQED_finalize_BREIT_WHEELER


      !*********************************!
      !* QUANTUM PARAMETER CALCULATION *!
      !*********************************!
      FUNCTION compute_chi_with_vectors(gamma, p_in, EE, BB) RESULT(chi) BIND(C, NAME="compute_chi_with_vectors")
        IMPORT :: C_DOUBLE
        REAL(C_DOUBLE) :: chi
        REAL(C_DOUBLE), INTENT(IN) :: gamma
        REAL(C_DOUBLE), DIMENSION(3), INTENT(IN) :: p_in, EE, BB
      END FUNCTION compute_chi_with_vectors

      FUNCTION compute_chi_with_components(gamma, p_in_x, p_in_y, p_in_z, Ex, Ey, Ez, Bx, By, Bz) &
        RESULT(chi) BIND(C, NAME="compute_chi_with_components")
        IMPORT :: C_DOUBLE
        REAL(C_DOUBLE) :: chi
        REAL(C_DOUBLE), INTENT(IN) :: gamma, p_in_x, p_in_y, p_in_z, Ex, Ey, Ez, Bx, By, Bz
      END FUNCTION compute_chi_with_components

      !*****************!
      !* BUILD MOMENTA *!
      !*****************!
      SUBROUTINE SFQED_build_collinear_momentum(gamma_out, p_in, p_out) BIND(C, NAME="SFQED_build_collinear_momentum")
        IMPORT :: C_DOUBLE
        REAL(C_DOUBLE), DIMENSION(3), INTENT(IN) :: p_in
        REAL(C_DOUBLE), DIMENSION(3), INTENT(OUT) :: p_out
        REAL(C_DOUBLE), INTENT(IN) :: gamma_out
      END SUBROUTINE SFQED_build_collinear_momentum


      !********************************!
      !* PHOTON EMISSION SECTION LCFA *!
      !********************************!
      FUNCTION SFQED_INV_COMPTON_rate(gamma, chi) RESULT(rate) BIND(C, NAME="SFQED_INV_COMPTON_rate")
        IMPORT :: C_DOUBLE
        REAL(C_DOUBLE) :: rate
        REAL(C_DOUBLE), INTENT(IN) :: gamma, chi
      END FUNCTION SFQED_INV_COMPTON_rate

      FUNCTION SFQED_LCFA_INV_COMPTON_PHOTON_energy(gamma, chi, rnd) &
         RESULT(nrg) BIND(C, NAME="SFQED_LCFA_INV_COMPTON_PHOTON_energy")
        IMPORT :: C_DOUBLE
        REAL(C_DOUBLE) :: nrg
        REAL(C_DOUBLE), INTENT(IN) :: gamma, chi, rnd
      END FUNCTION SFQED_LCFA_INV_COMPTON_PHOTON_energy

      !*********************************!
      !* PHOTON EMISSION SECTION BLCFA *!
      !*********************************!

      !IMPORTANT
      ! when dealing with C_PTR passed as argument remember to 
      ! specify the VALUE attribute in the interface. This variable
      ! stores the value of the c++ pointer (an address), and if you do not 
      ! use this attribute the fortran function will by default pass the
      ! argument by reference.
      function SFQED_CREATE_BLCFA_OBJECT() result(this) bind(C,name="SFQED_CREATE_BLCFA_OBJECT")
        import C_PTR
        type(C_PTR) :: this
      end function SFQED_CREATE_BLCFA_OBJECT

      SUBROUTINE SFQED_FINALIZE_BLCFA_OBJECT(entity) bind(C,name="SFQED_FINALIZE_BLCFA_OBJECT")
        import C_PTR
        type(C_PTR), INTENT(IN), VALUE :: entity
      end SUBROUTINE SFQED_FINALIZE_BLCFA_OBJECT

      function SFQED_BLCFA_OBJECT_update(entity, pushed_momentum, momentum, delta, &
        part_gamma, part_chi) result(goon) bind(C,name="SFQED_BLCFA_OBJECT_update")
        import C_PTR, C_DOUBLE, C_BOOL
        LOGICAL(C_BOOL) :: goon
        TYPE(C_PTR), INTENT(IN), VALUE :: entity
        REAL(C_DOUBLE), DIMENSION(3), INTENT(IN) :: pushed_momentum
        REAL(C_DOUBLE), DIMENSION(3), INTENT(IN) :: momentum
        ! REAL(C_DOUBLE), DIMENSION(3), INTENT(INOUT) :: Lorentz_F_t_der
        ! REAL(C_DOUBLE), DIMENSION(3), INTENT(INOUT) :: Lorentz_F_tt_der
        REAL(C_DOUBLE), INTENT(INOUT) :: delta
        REAL(C_DOUBLE), INTENT(INOUT) :: part_gamma
        REAL(C_DOUBLE), INTENT(INOUT) :: part_chi
      end function SFQED_BLCFA_OBJECT_update

      function SFQED_BLCFA_OBJECT_update_raw(pushed_momentum, momentum, Lorentz_F_old, Delta_Lorentz_F_old, &
        just_created, delta, part_gamma, part_chi) result(goon) bind(C,name="SFQED_BLCFA_OBJECT_update_raw")
        import C_PTR, C_DOUBLE, C_BOOL
        LOGICAL(C_BOOL) :: goon
        REAL(C_DOUBLE), DIMENSION(3), INTENT(IN) :: pushed_momentum
        REAL(C_DOUBLE), DIMENSION(3), INTENT(IN) :: momentum
        REAL(C_DOUBLE), DIMENSION(3), INTENT(INOUT) :: Lorentz_F_old
        REAL(C_DOUBLE), DIMENSION(3), INTENT(INOUT) :: Delta_Lorentz_F_old
        LOGICAL(C_BOOL), INTENT(INOUT) :: just_created
        REAL(C_DOUBLE), INTENT(INOUT) :: delta
        REAL(C_DOUBLE), INTENT(INOUT) :: part_gamma
        REAL(C_DOUBLE), INTENT(INOUT) :: part_chi
      end function SFQED_BLCFA_OBJECT_update_raw
      
      function SFQED_BLCFA_INV_COMPTON_PHOTON_threshold(entity, delta, &
        part_gamma, part_chi) result(thresh) bind(C,name="SFQED_BLCFA_INV_COMPTON_PHOTON_threshold")
        import C_PTR, C_DOUBLE
        REAL(C_DOUBLE) :: thresh
        TYPE(C_PTR), INTENT(IN), VALUE :: entity
        ! REAL(C_DOUBLE), DIMENSION(3), INTENT(IN) :: Lorentz_F_t_der
        ! REAL(C_DOUBLE), DIMENSION(3), INTENT(IN) :: Lorentz_F_tt_der
        REAL(C_DOUBLE), INTENT(IN) :: delta
        REAL(C_DOUBLE), INTENT(IN) :: part_gamma
        REAL(C_DOUBLE), INTENT(IN) :: part_chi
      end function SFQED_BLCFA_INV_COMPTON_PHOTON_threshold

      function SFQED_BLCFA_INV_COMPTON_PHOTON_threshold_2(tau, &
        part_gamma, part_chi) result(thresh) bind(C,name="SFQED_BLCFA_INV_COMPTON_PHOTON_threshold_2")
        import C_PTR, C_DOUBLE
        REAL(C_DOUBLE) :: thresh
        ! REAL(C_DOUBLE), DIMENSION(3), INTENT(IN) :: Lorentz_F_t_der
        ! REAL(C_DOUBLE), DIMENSION(3), INTENT(IN) :: Lorentz_F_tt_der
        REAL(C_DOUBLE), INTENT(IN) :: tau
        REAL(C_DOUBLE), INTENT(IN) :: part_gamma
        REAL(C_DOUBLE), INTENT(IN) :: part_chi
      end function SFQED_BLCFA_INV_COMPTON_PHOTON_threshold_2

      function SFQED_BLCFA_INV_COMPTON_PHOTON_threshold_raw(Lorentz_F, delta, &
        part_gamma, part_chi) result(thresh) bind(C,name="SFQED_BLCFA_INV_COMPTON_PHOTON_threshold_raw")
        import C_PTR, C_DOUBLE
        REAL(C_DOUBLE) :: thresh
        REAL(C_DOUBLE), DIMENSION(3), INTENT(IN) :: Lorentz_F
        REAL(C_DOUBLE), INTENT(IN) :: delta
        REAL(C_DOUBLE), INTENT(IN) :: part_gamma
        REAL(C_DOUBLE), INTENT(IN) :: part_chi
      end function SFQED_BLCFA_INV_COMPTON_PHOTON_threshold_raw

      function SFQED_BLCFA_INV_COMPTON_PHOTON_energy(LCFA_limit, gamma, chi, rnd, rnd2) &
        result(energy) bind(C,name="SFQED_BLCFA_INV_COMPTON_PHOTON_energy")
        import C_DOUBLE
        REAL(C_DOUBLE) :: energy
        REAL(C_DOUBLE), INTENT(IN) :: LCFA_limit
        REAL(C_DOUBLE), INTENT(IN) :: gamma
        REAL(C_DOUBLE), INTENT(IN) :: chi
        REAL(C_DOUBLE), INTENT(IN) :: rnd
        REAL(C_DOUBLE), INTENT(IN) :: rnd2
      end function SFQED_BLCFA_INV_COMPTON_PHOTON_energy

      function SFQED_BLCFA_INV_COMPTON_PHOTON_energy_2(LCFA_limit, gamma_photon, gamma, chi, rnd2) &
        result(energy) bind(C,name="SFQED_BLCFA_INV_COMPTON_PHOTON_energy_2")
        import C_DOUBLE
        REAL(C_DOUBLE) :: energy
        REAL(C_DOUBLE), INTENT(IN) :: LCFA_limit
        REAL(C_DOUBLE), INTENT(IN) :: gamma_photon
        REAL(C_DOUBLE), INTENT(IN) :: gamma
        REAL(C_DOUBLE), INTENT(IN) :: chi
        REAL(C_DOUBLE), INTENT(IN) :: rnd2
      end function SFQED_BLCFA_INV_COMPTON_PHOTON_energy_2

      !!!!!DEBUG!!!!!
      SUBROUTINE SFQED_BLCFA_get_F_old(entity, F_old) bind(C,name="SFQED_BLCFA_get_F_old")
        import C_PTR, C_DOUBLE
        REAL(C_DOUBLE), DIMENSION(3), INTENT(OUT) :: F_old
        TYPE(C_PTR), INTENT(IN), VALUE :: entity
      end SUBROUTINE SFQED_BLCFA_get_F_old

      SUBROUTINE SFQED_BLCFA_update_F_old(entity, F_old_new) bind(C,name="SFQED_BLCFA_update_F_old")
        import C_PTR, C_DOUBLE
        REAL(C_DOUBLE), DIMENSION(3), INTENT(IN) :: F_old_new
        TYPE(C_PTR), INTENT(IN), VALUE :: entity
      end SUBROUTINE SFQED_BLCFA_update_F_old

      SUBROUTINE SFQED_BLCFA_get_delta_F_old(entity, delta_F_old) bind(C,name="SFQED_BLCFA_get_delta_F_old")
        import C_PTR, C_DOUBLE
        REAL(C_DOUBLE), DIMENSION(3), INTENT(OUT) :: delta_F_old
        TYPE(C_PTR), INTENT(IN), VALUE :: entity
      end SUBROUTINE SFQED_BLCFA_get_delta_F_old

      SUBROUTINE SFQED_BLCFA_printALL(entity) bind(C,name="SFQED_BLCFA_printALL")
        import C_PTR
        TYPE(C_PTR), INTENT(IN), VALUE :: entity
      end SUBROUTINE SFQED_BLCFA_printALL


      !***************************!
      !* PAIR PRODUCTION SECTION *!
      !***************************!
      FUNCTION SFQED_BREIT_WHEELER_rate(gamma, chi) RESULT(rate) BIND(C, NAME="SFQED_BREIT_WHEELER_rate")
        IMPORT :: C_DOUBLE
        REAL(C_DOUBLE) :: rate
        REAL(C_DOUBLE), INTENT(IN) :: gamma, chi
      END FUNCTION SFQED_BREIT_WHEELER_rate

      FUNCTION SFQED_BREIT_WHEELER_rate_fast(gamma, chi) RESULT(rate) BIND(C, NAME="SFQED_BREIT_WHEELER_rate_fast")
        IMPORT :: C_DOUBLE
        REAL(C_DOUBLE) :: rate
        REAL(C_DOUBLE), INTENT(IN) :: gamma, chi
      END FUNCTION SFQED_BREIT_WHEELER_rate_fast

      FUNCTION SFQED_BREIT_WHEELER_ELECTRON_energy(gamma, chi, rnd) RESULT(nrg) BIND(C, NAME="SFQED_BREIT_WHEELER_ELECTRON_energy")
        IMPORT :: C_DOUBLE
        REAL(C_DOUBLE) :: nrg
        REAL(C_DOUBLE), INTENT(IN) :: gamma, chi, rnd
      END FUNCTION SFQED_BREIT_WHEELER_ELECTRON_energy

      FUNCTION SFQED_BREIT_WHEELER_ELECTRON_energy_fast(gamma, chi, rnd) RESULT(nrg) &
         BIND(C, NAME="SFQED_BREIT_WHEELER_ELECTRON_energy_fast")
        IMPORT :: C_DOUBLE
        REAL(C_DOUBLE) :: nrg
        REAL(C_DOUBLE), INTENT(IN) :: gamma, chi, rnd
      END FUNCTION SFQED_BREIT_WHEELER_ELECTRON_energy_fast


    END INTERFACE

    contains

    ! FUNCTION GET_PATH_TO_SFQEDtoolkit_COEFFICIENTS(path) RESULT(success)
    !   LOGICAL :: success
    !   CHARACTER(len=255), INTENT(OUT) :: path
    !   INTEGER :: len
    !   CALL get_environment_variable("SFQED_TOOLKIT_USER", path)
    !   path = TRIM(path) // "/coefficients/"
    !   ! you need to add the null character to your fortran string
    !   ! if you want to pass it to C function
    !   len = len_trim(path)
    !   path(len+1:len+1) = C_NULL_CHAR
    !   WRITE (*,*) 'Loading SFQED coefficients from ', path
    !   success = .TRUE.
    ! END FUNCTION GET_PATH_TO_SFQEDtoolkit_COEFFICIENTS
  
  END MODULE SFQEDtoolkit_Interface











