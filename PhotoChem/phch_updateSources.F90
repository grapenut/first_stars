


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!
!! This function sets up the point source coordinates. Below is my method, but you will certainly have your own method.
!! It would probably be a good idea to make a separate STAR_PART_TYPE to use for star particles, but I was lazy and
!! just reused the TAG_PART_PROP to correspond to the source's ID, but this requires tampering with other routines that set
!! particle tags.
!!
!! Set phch_num_sources to the number of point sources (probably just 1)
!! For each source s, set phch_source_x(s), phch_source_y(s) and phch_source_z(s) accordingly
!!
!! Currently they all use the same source luminosity/spectrum. If you want to use multiple sources, talk to me first!
!!!!!!!!!!!!!!!!!!

subroutine phch_updateSources()

  use Particles_data, ONLY : particles, pt_numLocal
  use PhotoChem_data, ONLY : phch_source_x, phch_source_y, phch_source_z, phch_ismaster
  use tree, ONLY : grid_xmin, grid_xmax, grid_ymin, grid_ymax, grid_zmin, grid_zmax

implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer :: tag, ierr, p
  
  real :: local_source_x, local_source_y, local_source_z

  ! fixed source in the center of the box
  !phch_source_x = (grid_xmax-grid_xmin)*0.5
  !phch_source_y = (grid_ymax-grid_ymin)*0.5
  !phch_source_z = (grid_zmax-grid_zmin)*0.5

  local_source_x = 0.0
  local_source_y = 0.0
  local_source_z = 0.0

  ! iterate through local particles to see if we have a source
  do p = 1, pt_numLocal
  
    if (particles(TYPE_PART_PROP,p).eq.real(DARK_PART_TYPE) .and. particles(TAG_PART_PROP,p).eq.0.0) then
    
      local_source_x = particles(POSX_PART_PROP,p)
      local_source_y = particles(POSY_PART_PROP,p)
      local_source_z = particles(POSZ_PART_PROP,p)
    
    end if
  
  end do
  
  ! I've commented this so you can use the phch_default_source_* method. If you need to use update_source_coords()
  !  you'll want to collect the coordinates of all sources onto all processors, something like what I have here

  call MPI_ALLREDUCE(local_source_x, phch_source_x, 1, FLASH_REAL, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(local_source_y, phch_source_y, 1, FLASH_REAL, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(local_source_z, phch_source_z, 1, FLASH_REAL, MPI_MAX, MPI_COMM_WORLD, ierr)
  
  ! double check that we found it!
  if (phch_source_x .le. 0.0 .or. phch_source_y .le. 0.0 .or. phch_source_z .le. 0.0) then
    ! one of the coordinates is bad!
    if (phch_ismaster) print *,'PHCH: ERROR didnt find a source!', phch_source_x, phch_source_y, phch_source_z
    phch_source_x = (grid_xmax-grid_xmin)*0.5
    phch_source_y = (grid_ymax-grid_ymin)*0.5
    phch_source_z = (grid_zmax-grid_zmin)*0.5
  end if
  
  if (phch_ismaster) print *,'PHCH: Source updated', phch_source_x, phch_source_y, phch_source_z
  
end subroutine phch_updateSources

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





