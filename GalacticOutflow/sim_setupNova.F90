

subroutine sim_setupNova(n)

  use PhotoChem_data, ONLY : phch_source_x, phch_source_y, phch_source_z
  
  use SImulation_data, ONLY : sim_ismaster, sim_starburst_radius, sim_star_coords, sim_mype, &
           sim_nova_radius, sim_nova_velocity, sim_star_mass, M_Sun, sim_newton, sim_pi, &
           sim_nova_r2, sim_zoom_level, sim_box_size, sim_log2, sim_zoom_level_offset

  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkBoundBox, Grid_getPointData, &
           Grid_getBlkIDFromPos, Grid_getBlkIndexLimits
  
implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(IN) :: n
  real, dimension(3) :: coords
  integer :: ierr
  
  integer :: blk, proc
  real :: dist
  
  integer, dimension(MDIM) :: intPos, guard
  real, dimension(MDIM) :: deltas, hfunc
  real, dimension(LOW:HIGH,MDIM) :: bndBox
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  
  real :: density, metallicity
  real :: redshift, invscale
  real :: nova_radius, nova_velocity
  
  !! generate a set of random coordinates on the master process
  if (sim_ismaster) then
    call Cosmology_getRedshift(redshift)
    invscale = 1.0 + redshift

    do
      call random_number(coords)
      coords = sim_starburst_radius * (coords*2.0 - 1.0)
    
      dist = sqrt(sum(coords*coords))

      print *,'SETUPNOVA: Picked coordinates: ', dist, coords(:)

      if (dist.gt.sim_starburst_radius) cycle
      
      ! scale into comoving coordinates and offset from the source position
      coords = coords * invscale
      coords(1) = phch_source_x + coords(1)
      coords(2) = phch_source_y + coords(2)
      coords(3) = phch_source_z + coords(3)

      exit
    end do
    
    print *,'SETUPNOVA: Scaled coordinates: ', n, coords(:)
  end if
  
  !! send the coordinates to everyone else
  call MPI_BCAST(coords, 3, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  
  sim_star_coords(n,:) = coords(:)
  
  !! find the cell density so we can set the nova_radius and velocity
  call Grid_getBlkIDFromPos(coords, blk, proc, MPI_COMM_WORLD)
  
  if (sim_ismaster .and. proc.lt.0) then
    print *,'DEBUG: coordinates not found ', proc, blk
  end if
  
  if (sim_mype.eq.proc) then
    call Grid_getDeltas(blk, deltas)
    call Grid_getBlkBoundBox(blk, bndBox)
    call Grid_getBlkIndexLimits(blk, blkLimits, blkLimitsGC)

    !! Get the integer indices of the cell containing the chosen coords
    !! using hfunc as the temporary storage for offset from the block bdry
    guard=blkLimits(LOW,:)-blkLimitsGC(LOW,:)
    hfunc(1:NDIM) = (coords(:) - bndBox(LOW,1:NDIM)) / deltas(1:NDIM)
    intPos(1:NDIM)= floor(hfunc(1:NDIM)) + 1 + guard(1:NDIM)
  
    call Grid_getPointData(blk, CENTER, DENS_VAR, EXTERIOR, intPos, density)
    !call Grid_getPointData(blk, CENTER, Z_SPEC, EXTERIOR, intPos, metallicity)
    
    ! calculate supernova radius and shock velocity
    nova_radius = (0.1 * sim_star_mass(n) * M_Sun / density * 3.0 / 4.0 / sim_pi)**(1.0/3.0)
    nova_velocity = (10.0 / 3.0 * sim_newton * sim_star_mass(n) * M_Sun / nova_radius)**(0.5)
  end if

  ! distribute the calculated values
  call MPI_BCAST(nova_radius, 1, FLASH_REAL, proc, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nova_velocity, 1, FLASH_REAL, proc, MPI_COMM_WORLD, ierr)
  
  ! store them (keep in comoving coordinates for easy comparison to the grid)
  sim_nova_radius(n) = nova_radius
  sim_nova_r2(n) = nova_radius * nova_radius
  sim_zoom_level(n) = floor(1.0 - log(nova_radius / sim_box_size) / sim_log2) + sim_zoom_level_offset
  sim_nova_velocity(n) = nova_velocity
    
  if (sim_ismaster) print *, 'SETUPNOVA: Radius / Velocity: ', sim_nova_radius(n), sim_nova_velocity(n)
  
end subroutine sim_setupNova

