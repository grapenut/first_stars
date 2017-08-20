
! initialize the grid data to zero, or the optionally supplied zeroVal

subroutine krn_zeroGridVar(zeroVal)

  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPtr, Grid_releaseBlkPtr

  implicit none

#include "constants.h"
  real, intent(in), optional :: zeroVal
  
  integer, dimension(MAXBLOCKS) :: blockList
  integer :: blockCount, b, blockID
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  real :: zero

  zero = 0.0
  if (present(zeroVal)) zero = zeroVal
  
  call Grid_getListOfBlocks(LEAF, blockList, blockCount)
  
  do b = 1, blockCount
    blockID = blockList(b)
    call Grid_getBlkPtr(blockID, solnVec, CENTER)
    solnVec = zero
    call Grid_releaseBlkPtr(blockID, solnVec, CENTER)
  end do

end subroutine

