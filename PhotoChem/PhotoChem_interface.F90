
Module PhotoChem_interface

implicit none

#include "constants.h"
#include "Flash.h"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Pixel functions
  interface
    subroutine PhotoChem(blockCount, blockList, dt)
      implicit none
      real, INTENT(IN) :: dt
      integer, INTENT(IN) :: blockCount
      integer, dimension(blockCount), INTENT(IN) :: blockList
    end subroutine PhotoChem
  end interface

  interface
    subroutine PhotoChem_init()
    implicit none
    end subroutine PhotoChem_init
  end interface

end Module PhotoChem_interface

