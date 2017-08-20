
! call Simulation_initBlock and Eos_wrapped to setup the other block variables

subroutine krn_initBlocks()

  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkIndexLimits
  use Simulation_interface, ONLY : Simulation_initBlock
  use Eos_interface, ONLY : Eos_wrapped

  implicit none
  
#include "constants.h"
  
  integer, dimension(MAXBLOCKS) :: blockList
  integer :: blockCount, b, blockID
  integer, dimension(2,MDIM) :: blockLimits
  integer, dimension(2,MDIM) :: blockLimitsGC

  call Grid_getListOfBlocks(LEAF, blockList, blockCount)
  do b=1,blockCount
    blockID = blockList(b)
    call Grid_getBlkIndexLimits(blockID, blockLimits, blockLimitsGC)
    call Simulation_initBlock(blockID)
    call Eos_wrapped(MODE_DENS_TEMP, blockLimits, blockID)
  end do
  
end subroutine

