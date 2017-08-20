
subroutine sim_markInitial()

  use Simulation_data, ONLY : sim_ismaster, sim_init_refine_done
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use tree, ONLY : nodetype, bsize, coord, refine, lrefine, lnblocks

implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  logical, save :: first_call = .true.
  real, save :: center_x, center_y, center_z, max_distance
  real, save :: xmin, xmax, ymin, ymax, zmin, zmax
  integer, save :: center_refine_level

  real, dimension(MDIM) :: blockCenter, blockSize
  real :: bxl, bxr, byl, byr, bzl, bzr, xmindist2, ymindist2, zmindist2, radius

  integer :: b, ierr
  logical :: refine_local, refine_global


  if(first_call) then
     
     call RuntimeParameters_get("xmin", xmin)
     call RuntimeParameters_get("xmax", xmax)
     call RuntimeParameters_get("ymin", ymin)
     call RuntimeParameters_get("ymax", ymax)
     call RuntimeParameters_get("zmin", zmin)
     call RuntimeParameters_get("zmax", zmax)
     
     center_x = 0.5*(xmax-xmin)
     center_y = 0.5*(ymax-ymin)
     center_z = 0.5*(zmax-zmin)
     center_refine_level = 10
     
     max_distance = 3.0e21
  end if
  
  ! Get bounding box for this block (relative to (center_x,center_y,center_z))
  do b = 1, lnblocks
    if (nodetype(b) == LEAF .and. lrefine(b) .lt. center_refine_level) then
      blockCenter = coord(:,b)
      blockSize  = 0.5 * bsize(:,b)
       
      bxl = blockCenter(1) - blockSize(1) - center_x
      bxr = blockCenter(1) + blockSize(1) - center_x
      if (NDIM > 1) then
        byl = blockCenter(2) - blockSize(2) - center_y
        byr = blockCenter(2) + blockSize(2) - center_y
      else
        byl = 0.
        byr = 0.
      endif
      if (NDIM == 3) then
        bzl = blockCenter(3) - blockSize(3) - center_z
        bzr = blockCenter(3) + blockSize(3) - center_z
      else
        bzl = 0.
        bzr = 0.
      endif
       
      ! Find minimum and maximum distance from (center_x,center_y,center_z) for each dimension.
      if (bxl*bxr > 0.) then
        xmindist2 = min( bxl**2, bxr**2 )
      else
        xmindist2 = 0.
      endif

      if (byl*byr > 0.) then
        ymindist2 = min( byl**2, byr**2 )
      else
        ymindist2 = 0.
      endif
       
      if (bzl*bzr > 0.) then
        zmindist2 = min( bzl**2, bzr**2 )
      else
        zmindist2 = 0.
      endif

      ! restrict max refinement by radius
      radius = sqrt(max(xmindist2 + ymindist2 + zmindist2,1.0))
      if (radius .lt. max_distance) then
        refine(b) = .true.
      end if

    end if ! LEAF
  end do

  ! check if the grid is ready
  refine_local = any(refine(:))
  call MPI_ALLREDUCE(refine_local, refine_global, 1, FLASH_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)

  if (.not.refine_global) then
    if (sim_ismaster) print *,'SIM: The grid is zoomed in on the densest particle.'
    sim_init_refine_done = .true.
  end if

end subroutine


