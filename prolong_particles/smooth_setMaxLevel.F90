
! not currently used
subroutine smooth_prepareNodes(level)
  
  use Grid_data, ONLY : gr_meshNumProcs, gr_meshMe
  use tree, ONLY : nodetype, lrefine, lnblocks, newchild
  use smooth_Data

#include "constants.h"  
#include "Flash.h" 
  
  
  implicit none
  
  integer, intent(IN) :: level
  
  integer :: lb
  
  
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
  
  call amr_get_new_nodetypes(gr_meshNumProcs,gr_meshMe,level)
  
  
  return
  
end subroutine smooth_prepareNodes
