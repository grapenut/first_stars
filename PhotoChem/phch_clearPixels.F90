
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine phch_clearPixels(numblocks, blocklist)

  use Grid_interface, ONLY : Grid_putPointData, Grid_getBlkIndexLimits
  
implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(in) :: numblocks
  integer, intent(in), dimension(numblocks) :: blocklist
  
  integer :: b, blockId, i, j, k
  integer, dimension(MDIM) :: axis
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  
  do b=1,numblocks
    blockId = blocklist(b)

    call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

    do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)

      axis(KAXIS) = k
       
      do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)

        axis(JAXIS) = j

        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

          axis(IAXIS) = i

          call Grid_putPointData(blockId, CENTER, FLUX_VAR, EXTERIOR, axis, 0.0)
          !call Grid_putPointdata(blockId, CENTER, XI_VAR, EXTERIOR, axis, 0.0)

        enddo

      enddo

    enddo

  end do

  return

end subroutine phch_clearPixels

