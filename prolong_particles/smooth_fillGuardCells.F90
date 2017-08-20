
subroutine smooth_fillGuardCells(level, ivar)
  

  use Grid_data, ONLY : gr_meshMe
  use physicaldata, ONLY : unk
  use workspace, ONLY : work
  use tree, ONLY : nodetype, lrefine, lnblocks, newchild
  use smooth_Data
  use Grid_interface, ONLY: Grid_getBlkPtr, Grid_releaseBlkPtr
  use gr_interface, ONLY : gr_setGcFillNLayers
  use paramesh_interfaces, ONLY: amr_guardcell


  implicit none
  
#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: level, ivar
  
  real,pointer :: solnData(:,:,:,:)
  integer ::  lb, i, j, k
  integer,dimension(MDIM) :: layers
  
  
  do lb = 1, lnblocks
     call Grid_getBlkPtr(lb,solnData)
     if (((level == 0) .and. (nodetype(lb) == LEAF)) .or. &
          (lrefine(lb) == level) .or. &
          (lrefine(lb) == level+1) .or. (lrefine(lb) == level-1)) then
        do k = hg_kli, hg_kui
           do j = hg_jli, hg_jui
              do i = hg_ili, hg_iui
                 work(i,j,k,lb,1) = solnData(ivar,i,j,k)
              enddo
           enddo
        enddo
     endif
     call Grid_releaseBlkPtr(lb,solnData)
  enddo
  
  ! guard cell part:
  

  call gr_setGcFillNLayers(layers,ALLDIR,NGUARD,minLayers=NGUARD)
  call amr_guardcell(gr_meshMe, 2, NGUARD, layers(IAXIS), layers(JAXIS), layers(KAXIS))
  
  
  do lb = 1, lnblocks
     if (  ( (level == 0) .and. (nodetype(lb) == LEAF) ) .or. &
          (lrefine(lb) == level) .or. &
          (lrefine(lb) == level+1) .or. & 
          (lrefine(lb) == level-1) ) then
        
        unk(ivar,hg_ile:hg_ili-1,:,:,lb) = work(hg_ile:hg_ili-1,:,:,lb,1)
        unk(ivar,hg_iui+1:hg_iue,:,:,lb) = work(hg_iui+1:hg_iue,:,:,lb,1)
        
        unk(ivar,:,hg_jle:hg_jli-1,:,lb) = work(:,hg_jle:hg_jli-1,:,lb,1)
        unk(ivar,:,hg_jui+1:hg_jue,:,lb) = work(:,hg_jui+1:hg_jue,:,lb,1)
        
        unk(ivar,:,:,hg_kle:hg_kli-1,lb) = work(:,:,hg_kle:hg_kli-1,lb,1)
        unk(ivar,:,:,hg_kui+1:hg_kue,lb) = work(:,:,hg_kui+1:hg_kue,lb,1)
        
     endif
  enddo
  

  return

end subroutine smooth_fillGuardCells
