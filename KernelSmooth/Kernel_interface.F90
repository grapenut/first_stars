module Kernel_interface

  implicit none

#include "constants.h"

  interface
    subroutine Kernel_mapParticles(partType, zero_gridVar)
      integer, intent(in) :: partType
      logical, intent(in), optional :: zero_gridVar
    end subroutine
  end interface

  interface
    subroutine krn_calculateWeights()
    end subroutine
  end interface

  interface
    recursive subroutine krn_checkBounds(px, py, pz, len, bndBox, overlap, again)
      real, intent(in) :: px, py, pz, len
      real, intent(in), dimension(LOW:HIGH, MDIM) :: bndBox
      logical, intent(inout) :: overlap
      logical, intent(in) :: again
    end subroutine
  end interface

  interface
    subroutine krn_divideDensity()
    end subroutine
  end interface

  interface
    subroutine krn_findParticles()
    end subroutine
  end interface

  interface
    recursive subroutine krn_getCellWeight(px, py, pz, cx, cy, cz, len, wgt, again)
      real, intent(in) :: px, py, pz, cx, cy, cz, len
      real, intent(inout) :: wgt
      logical, intent(in) :: again  !check to prevent multiple levels of recursion
    end subroutine
  end interface

  interface
    subroutine krn_gridTotals()
    end subroutine
  end interface

  interface
    subroutine krn_initBlocks()
    end subroutine
  end interface

  interface
    subroutine krn_init(partType)
      integer, intent(in) :: partType
    end subroutine
  end interface

  interface
    subroutine krn_kernel(q,w)
      real, intent(in) :: q
      real, intent(inout) :: w
    end subroutine
  end interface

  interface
    subroutine krn_makeParticle(x, y, z, mass, px, py, pz, temp, len)
      real, intent(in) :: x, y, z, mass, px, py, pz, temp, len
    end subroutine
  end interface

  interface
    subroutine krn_moveParticles()
    end subroutine
  end interface

  interface
    subroutine krn_partTotals()
    end subroutine
  end interface

  interface
    subroutine krn_smoothParticles()
    end subroutine
  end interface

  interface
    subroutine krn_sortParticles()
    end subroutine
  end interface

  interface
    subroutine krn_timestamp(msg)
      character (len=*), intent(in) :: msg
    end subroutine
  end interface

  interface
    subroutine krn_zeroGridVar(zeroVal)
      real, intent(in), optional :: zeroVal
    end subroutine
  end interface

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine krn_free()
    
    use Kernel_data, ONLY : krn_particles, krn_recvBuf

    implicit none
    
    deallocate(krn_particles)
    deallocate(krn_recvBuf)
    
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module

