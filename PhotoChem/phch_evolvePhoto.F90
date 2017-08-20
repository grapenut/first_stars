subroutine phch_evolvePhoto(y, photoeq, dens, eint, dt)

  use PhotoChem_data, ONLY : phch_iTemp, phch_num_xi, phch_xi_list, phch_usePhoto, &
                             phch_num_ximetal, phch_ximetal_list, phch_spec_list, phch_xi_table

implicit none

#include "Flash.h"
#include "constants.h"

  real, INTENT(INOUT) :: y(SPECIES_BEGIN:SPECIES_END+1), eint
  real, INTENT(IN) :: photoeq, dt, dens

  real :: nH, metal
  
  integer :: ix, jx
  
  real :: xi0, metal0
  
  real, dimension(SPECIES_BEGIN:SPECIES_END+1) :: y_eq
  
  real :: t_eq
  
  if (.not.phch_usePhoto) return

  ! number density of hydrogen atoms
  nH = y(H_SPEC) + y(HPLU_SPEC)
  
  ! mass fraction of metals
  metal = y(Z_SPEC)
   
  ! table indices
  ix = max(1,min(phch_num_xi-1,maxloc(phch_xi_list, 1, phch_xi_list.le.photoeq)))
  jx = max(1,min(phch_num_ximetal-1,maxloc(phch_ximetal_list, 1, phch_ximetal_list.le.metal)))

  ! interpolation scalars
  xi0 = (photoeq - phch_xi_list(ix)) / (phch_xi_list(ix+1) - phch_xi_list(ix))
  metal0 = (metal - phch_ximetal_list(jx)) / (phch_ximetal_list(jx+1) - phch_ximetal_list(jx))
  
  if (metal0.lt.0.0) metal0 = 0.0
  if (metal0.gt.1.0) metal0 = 1.0
  
  ! interpolate equilibrium values at the given coordinates
  call phch_interpPhoto(y_eq, ix, jx, xi0, metal0)
  
  ! rescale equilibrium abundances to local abundances
  y_eq(phch_spec_list) = y_eq(phch_spec_list) * nH

  ! recover the equilibrium electron abundance by conservation
  call phch_electronEquilibrium(y_eq)

  ! calc equilibrium timescale
  !t_eq = abs(eint / ( y_eq(phch_iTemp) * nH * nH / dens ))
  t_eq = dt
  
  ! relax exponentially towards the equilibrium solution
  ! if dt >= t_eq then go directly to y_eq
  ! else if dt < t_eq then scale as exp(-dt/t_eq)
  if (dt.ge.t_eq) then
    y(phch_spec_list) = y_eq(phch_spec_list)
    y(phch_iTemp) = y_eq(phch_iTemp)
    
    !print *,'PHCHDEBUG: ', y(phch_iTemp), phch_xi_table(ix,jx,phch_iTemp), photoeq
  else
    y(phch_spec_list) = y_eq(phch_spec_list) + exp(-dt/t_eq) * (y(phch_spec_list) - y_eq(phch_spec_list))
    y(phch_iTemp) = y_eq(phch_iTemp) + exp(-dt/t_eq) * (y(phch_iTemp) - y_eq(phch_iTemp))
  end if

  ! renormalize the relaxed solution

  ! calculate the electron abundance (without metal electron contribution)
  y(ELEC_SPEC) = y(HPLU_SPEC) + y(HEP_SPEC) + 2.0*y(HEPP_SPEC)
  !metal_electrons = y_eq(

end subroutine !! ch_evolvePhoto

