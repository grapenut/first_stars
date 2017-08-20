

subroutine smooth_prepareNodes(level, cautious)
  
  use Grid_data, ONLY : gr_meshNumProcs, gr_meshMe
  use tree, ONLY : nodetype, lrefine, lnblocks, newchild
  use smooth_Data
  Use paramesh_mpi_interfaces, only : amr_morton_process
  use physicaldata, only : surr_blks_valid
  
#include "constants.h"  
#include "Flash.h" 
  
  
  implicit none
  
  integer, intent(IN) :: level
  logical, intent(IN) :: cautious
  
  integer :: lb
  
  
  surr_blks_valid = .false.

  newchild(:) = .false.
  nodetype(1:lnblocks) = smooth_saveNodeType(1:lnblocks)
 
  do lb = 1, lnblocks
     if (lrefine(lb) > level)  nodetype(lb) = -1 
     if (lrefine(lb) == level) then 
        nodetype(lb) = LEAF
        newchild(lb) = .TRUE.
     endif
     if ((lrefine(lb) == level-1) .and. (nodetype(lb) /= LEAF)) & 
          nodetype(lb) = PARENT_BLK
  enddo
  

  if(cautious) then
     !if (gr_meshMe.eq.MASTER_PE) print*, "smooth_prepareNodes: calling morton process..."
     call amr_morton_process()

  else
     ! not sure if the following actually works...
     call amr_get_new_nodetypes(gr_meshNumProcs,gr_meshMe,level)
  endif
  
  
  return
  
end subroutine smooth_prepareNodes
