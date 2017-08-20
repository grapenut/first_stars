

subroutine smooth_restoreNodeTypes(cautious)

  
  use smooth_Data
  use Grid_data, ONLY : gr_meshNumProcs, gr_meshMe
  use tree, ONLY : lnblocks, newchild, nodetype
  Use paramesh_mpi_interfaces, only : amr_morton_process
  use physicaldata, only : surr_blks_valid

  implicit none
  
  logical, intent(IN) :: cautious
  
  integer :: lb
  
  surr_blks_valid = .false.

  do lb = 1, lnblocks
     nodetype(lb) = smooth_saveNodeType(lb)
     newchild(lb) = smooth_saveNewChild(lb)
  enddo

  if(cautious) then
     ! meant for use with particle smoothing
     call amr_morton_process()
  else
     ! meant for use with acceleration fetching
     call amr_get_new_nodetypes (gr_meshNumProcs, gr_meshMe, maxMeshRefineLevel)
  endif
  
  return
  
end subroutine smooth_restoreNodeTypes
