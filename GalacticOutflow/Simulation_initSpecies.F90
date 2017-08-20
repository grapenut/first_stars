!!****f* source/Simulation/Simulation_initSpecies
!!
!! NAME
!!
!!  Simulation_initSpecies
!!
!! SYNOPSIS
!!
!!  Simulation_initSpecies()
!!
!! DESCRIPTION
!!
!!  This routine will initialize the species and species values needed for a 
!!  given setup.   The user should add the 
!!  implementation of this routine to the setups directory of a simulation 
!!  that needs to use the multispecies capabilities of the code.
!!
!!  There two general purpose implementations available in the code, one which sets standard  
!!  isotope properties for the nuclear burning source terms, and another one for the 
!!  Ionization source term.
!!
!!  This routine is called from Multispecies_init, and is called BEFORE
!!  the call to Simulation_init.  
!!
!! SEE ALSO
!!  Multispecies_init
!!  Simulation/SimulationComposition/Simulation_initSpecies
!!
!!***

subroutine Simulation_initSpecies()

  use Multispecies_interface, ONLY : Multispecies_setProperty
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

implicit none
  
#include "Flash.h"
#include "Multispecies.h"
#include "constants.h"
  
  real :: emass, pmass, e_abar, metal_mass
  
  real :: gamma53, gamma75, nprotons
   
  gamma53 = 5.0/3.0
  gamma75 = 7.0/5.0
  
  call PhysicalConstants_get("proton mass",pmass)
  call PhysicalConstants_get("electron mass",emass)
  e_abar=emass/pmass

#ifdef Z_SPEC
  call Multispecies_setProperty(Z_SPEC, A, 14.0)
  call Multispecies_setProperty(Z_SPEC, Z, 7.0)
  call Multispecies_setProperty(Z_SPEC, GAMMA, gamma53)
  call Multispecies_setProperty(Z_SPEC, E, 7.0)
#endif

  ! glover's reduced chemistry network

!!$
  ! electrons:
  call Multispecies_setProperty(ELEC_SPEC, A, e_abar)
  call Multispecies_setProperty(ELEC_SPEC, Z, 0.0)
  call Multispecies_setProperty(ELEC_SPEC, GAMMA, gamma53)
  call Multispecies_setProperty(ELEC_SPEC, E, 1.0)
!!$
  ! neutral hydrogen
  call Multispecies_setProperty(H_SPEC, A, 1.0 + e_abar)
  call Multispecies_setProperty(H_SPEC, Z, 1.0)
  call Multispecies_setProperty(H_SPEC, GAMMA, gamma53)
  call Multispecies_setProperty(H_SPEC, E, 1.0)

  ! protons
  call Multispecies_setProperty(HPLU_SPEC, A, 1.0)
  call Multispecies_setProperty(HPLU_SPEC, Z, 1.0)
  call Multispecies_setProperty(HPLU_SPEC, GAMMA, gamma53)
  call Multispecies_setProperty(HPLU_SPEC, E, 0.0)
!!$  
#ifdef HMIN_SPEC
  ! H-minus
  call Multispecies_setProperty(HMIN_SPEC, A, 1.0 + 2.0*e_abar)
  call Multispecies_setProperty(HMIN_SPEC, Z, 1.0)
  call Multispecies_setProperty(HMIN_SPEC, GAMMA, gamma53)
  call Multispecies_setProperty(HMIN_SPEC, E, 2.0)
#endif
#ifdef HTWO_SPEC
  ! molecular hydrogen
  call Multispecies_setProperty(HTWO_SPEC, A, 2.0 + 2.0*e_abar)
  call Multispecies_setProperty(HTWO_SPEC, Z, 2.0)
  call Multispecies_setProperty(HTWO_SPEC, GAMMA, gamma75)
  call Multispecies_setProperty(HTWO_SPEC, E, 2.0)
#endif
!!$  
  ! neutral helium
  call Multispecies_setProperty(HEL_SPEC, A, 4.0 + 2.0*e_abar)
  call Multispecies_setProperty(HEL_SPEC, Z, 2.0)
  call Multispecies_setProperty(HEL_SPEC, GAMMA, gamma53)
  call Multispecies_setProperty(HEL_SPEC, E, 2.0)
  call Multispecies_setProperty(HEL_SPEC, N, 2.0)
!!$  
  ! singly ionized helium
  call Multispecies_setProperty(HEP_SPEC, A, 4.0 + e_abar)
  call Multispecies_setProperty(HEP_SPEC, Z, 2.0)
  call Multispecies_setProperty(HEP_SPEC, GAMMA, gamma53)
  call Multispecies_setProperty(HEP_SPEC, E, 1.0) 
  call Multispecies_setProperty(HEP_SPEC, N, 2.0)
  
  ! doubly ionized helium
  call Multispecies_setProperty(HEPP_SPEC, A, 4.0)
  call Multispecies_setProperty(HEPP_SPEC, Z, 2.0)
  call Multispecies_setProperty(HEPP_SPEC, GAMMA, gamma53)
  call Multispecies_setProperty(HEPP_SPEC, E, 0.0)
  call Multispecies_setProperty(HEPP_SPEC, N, 2.0)
!!$ 
#ifdef HTWP_SPEC
  ! molecular hyrogen plus
  call Multispecies_setProperty(HTWP_SPEC, A, 2.0 + e_abar)
  call Multispecies_setProperty(HTWP_SPEC, Z, 2.0)
  call Multispecies_setProperty(HTWP_SPEC, GAMMA, gamma75)
  call Multispecies_setProperty(HTWP_SPEC, E, 1.0)
#endif
!!$
#ifdef DEUT_SPEC
  ! deuterium
  call Multispecies_setProperty(DEUT_SPEC, A, 2.0 + e_abar)
  call Multispecies_setProperty(DEUT_SPEC, Z, 1.0)
  call Multispecies_setProperty(DEUT_SPEC, GAMMA, gamma53)
  call Multispecies_setProperty(DEUT_SPEC, E, 1.0)
  call Multispecies_setProperty(DEUT_SPEC, N, 1.0)
#endif
#ifdef DPLU_SPEC
  ! deuterium plus
  call Multispecies_setProperty(DPLU_SPEC, A, 2.0)
  call Multispecies_setProperty(DPLU_SPEC, Z, 1.0)
  call Multispecies_setProperty(DPLU_SPEC, GAMMA, gamma53)
  call Multispecies_setProperty(DPLU_SPEC, E, 0.0)
  call Multispecies_setProperty(DPLU_SPEC, N, 1.0)
#endif
#ifdef HD_SPEC
 ! HD
  call Multispecies_setProperty(HD_SPEC, A, 3.0 + 2.0*e_abar)
  call Multispecies_setProperty(HD_SPEC, Z, 2.0)
  call Multispecies_setProperty(HD_SPEC, GAMMA, gamma75)
  call Multispecies_setProperty(HD_SPEC, E, 2.0)
  call Multispecies_setProperty(HD_SPEC, N, 1.0)
#endif  

end subroutine Simulation_initSpecies
