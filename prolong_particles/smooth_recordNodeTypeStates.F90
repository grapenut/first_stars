
subroutine smooth_recordNodeTypeStates 
  
  use Grid_data, ONLY : gr_meshComm, gr_meshMe
  use tree, ONLY : lnblocks, newchild,nodetype, lrefine
  use smooth_Data
  
  implicit none

#include "Flash_mpi.h"

  integer :: lb, ierr
  integer :: mylrefmin, mylrefmax

  
  ! save node types
  do lb = 1, lnblocks
     smooth_saveNodeType(lb) = nodetype(lb)
     smooth_saveNewChild(lb) = newchild(lb)
  end do
  
  
  ! get max/min refinement level on grid
  mylrefmin = 100
  mylrefmax = 1
  
  do lb = 1, lnblocks
     mylrefmin = min(mylrefmin,lrefine(lb))
     mylrefmax = max(mylrefmax,lrefine(lb))
  end do
  
  call mpi_allreduce(mylrefmin, minMeshRefineLevel, 1, MPI_INTEGER, &
       MPI_MIN, gr_meshComm, ierr)
  call mpi_allreduce(mylrefmax, maxMeshRefineLevel, 1, MPI_INTEGER, & 
       MPI_MAX, gr_meshComm, ierr)
  
  return
  
end subroutine smooth_recordNodeTypeStates

  
