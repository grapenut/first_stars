

subroutine phch_electronEquilibrium(y)

implicit none

#include "Flash.h"

  real, dimension(SPECIES_BEGIN:SPECIES_END+1), intent(INOUT) :: y
  
  y(ELEC_SPEC) = y(HPLU_SPEC) + y(HEP_SPEC) + 2.0*y(HEPP_SPEC)

end subroutine


  