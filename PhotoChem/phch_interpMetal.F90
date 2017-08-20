
! for a given set of lower indices (kx, jx) and the scaled parameter values relative to those indices (xi0, metal0)
! return y as the interpolated state vector for the true parameters xi = xi(ix) + xi0 and metal = metal(jx) + metal0

subroutine phch_interpMetal(y, ix, jx, kx, dens0, metal0, temp0)

  use PhotoChem_data, ONLY : phch_iTemp, phch_cool_table, phch_helium_abundance

implicit none

#include "Flash.h"
#include "constants.h"

  real, INTENT(INOUT) :: y(SPECIES_BEGIN:SPECIES_END+1)
  integer, INTENT(IN) :: ix, jx, kx
  real, INTENT(IN) :: dens0, metal0, temp0

  real, dimension(SPECIES_BEGIN:SPECIES_END+1) :: y_eq, c0, c1, c00, c01, c10, c11
  
  real :: Hy_eq, Hel_eq
  
  ! trilinear interpolation (final values in y_eq)
  c00 = phch_cool_table(ix,jx,kx,:) * (1.0-dens0) + phch_cool_table(ix+1,jx,kx,:) * dens0
  c10 = phch_cool_table(ix,jx+1,kx,:) * (1.0-dens0) + phch_cool_table(ix+1,jx+1,kx,:) * dens0
  c01 = phch_cool_table(ix,jx,kx+1,:) * (1.0-dens0) + phch_cool_table(ix+1,jx,kx+1,:) * dens0
  c11 = phch_cool_table(ix,jx+1,kx+1,:) * (1.0-dens0) + phch_cool_table(ix+1,jx+1,kx+1,:) * dens0
  c0 = c00*(1.0-metal0) + c10*metal0
  c1 = c01*(1.0-metal0) + c11*metal0
  y_eq = c0*(1.0-temp0) + c1*temp0
  
  Hy_eq = y_eq(H_SPEC) + y_eq(HPLU_SPEC)
  Hel_eq = y_eq(HEL_SPEC) + y_eq(HEP_SPEC) + y_eq(HEPP_SPEC)

  y(H_SPEC) = y_eq(H_SPEC) / Hy_eq
  y(HPLU_SPEC) = y_eq(HPLU_SPEC) / Hy_eq
  y(HEL_SPEC) = y_eq(HEL_SPEC) * phch_helium_abundance / Hel_eq
  y(HEP_SPEC) = y_eq(HEP_SPEC) * phch_helium_abundance / Hel_eq
  y(HEPP_SPEC) = y_eq(HEPP_SPEC) * phch_helium_abundance / Hel_eq
  
  ! return the heating/cooling rate for timescale estimation
  y(phch_iTemp) = y_eq(phch_iTemp)

end subroutine

