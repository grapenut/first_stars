!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine PhotoChem(blockCount, blockList, dt)

  use PhotoChem_data, ONLY : phch_usePhoto
  use phch_interface, ONLY : phch_mapPixels, phch_clearPixels, phch_updateSources, phch_evolveChem
  !use Simulation_data, ONLY : sim_nova_started
  
implicit none

#include "Flash.h"

  real, INTENT(IN) :: dt
  integer, INTENT(IN) :: blockCount
  integer, dimension(blockCount), INTENT(IN) :: blockList
  
  integer :: n
  
  ! update source coordinates and intensities
  call phch_updateSources()

  if (phch_usePhoto) then

    ! reset data structures from the last run
    call phch_clearPixels(blockCount, blockList) 

    ! main radiation transport entry point
    ! maps grid density to healpix
    ! calculates fluxes on healpix
    ! maps healpix flux back to the grid
    call phch_mapPixels(blockCount, blockList)
  end if
  
  ! evolve chemistry (both photoionization and cooling)
  call phch_evolveChem(blockCount, blockList, dt)

end subroutine PhotoChem


