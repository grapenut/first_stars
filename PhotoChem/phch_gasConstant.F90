

subroutine phch_gasConstant(y, dens, gascon)

  use PhotoChem_data, ONLY : phch_spec_list, phch_particle_mass
  use Eos_data, ONLY : eos_gasConstant
  use eos_mgammaData, ONLY : eos_gammam1j
  use Multispecies_interface, ONLY : Multispecies_getSumInv

implicit none

#include "Flash.h"
#include "Multispecies.h"
#include "constants.h"

  real, dimension(SPECIES_BEGIN:SPECIES_END+1), intent(in) :: y
  real, intent(in) :: dens
  
  real, intent(out) :: gascon

  integer :: n
  real, dimension(NSPECIES) :: massFrac
  real :: abarInv, rt, gamma


  do n=SPECIES_BEGIN,SPECIES_END
    massFrac(1-SPECIES_BEGIN+n) = y(n)*phch_particle_mass(n) / dens
  end do
  massFrac(1-SPECIES_BEGIN+Z_SPEC) = y(Z_SPEC)
  
  !if (dr_myPE .eq. MASTER_PE) print *,'weights(',size(massFrac,1),') = ', massFrac
  
  call Multispecies_getSumInv(A, abarInv, massFrac(1:NSPECIES))
 
  massFrac(1:NSPECIES) = massFrac(1:NSPECIES)*eos_gammam1j(1:NSPECIES)
  call Multispecies_getSumInv(A, rt, massFrac(1:NSPECIES))

  gamma = 1.0e0 + abarInv/rt
  
  gascon = (gamma-1.0) / (abarInv*eos_gasConstant)


end subroutine

