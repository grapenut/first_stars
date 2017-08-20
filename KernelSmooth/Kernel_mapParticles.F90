
! external interface to map particles to grid using a kernel

subroutine Kernel_mapParticles(partType, zero_gridVar)

  use Kernel_data, ONLY : krn_numProcs
  use Kernel_interface
  use Particles_interface, ONLY : Particles_updateGridVar
  
  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: partType
  logical, intent(in), optional :: zero_gridVar
  logical :: zero
  integer :: p

  ! bypass smoothing and return early for testing
  !call krn_zeroGridVar(1.0e-40)
  !call krn_initBlocks()
  !return

  ! optional control to clear the grid
  zero = .true.
  if (present(zero_gridVar)) zero = zero_gridVar

  ! allocate memory and setup initial variables
  call krn_init(partType)
  call krn_timestamp('Started mapping.')
  
  ! sort particles and update indices
  call krn_sortParticles()
  call krn_timestamp('Sorted particles.')
  
  ! calculate total quantities to be conserved
  call krn_partTotals()
  
  ! locate smoothed particles and cull them from the regular list
  call krn_findParticles()
  call krn_timestamp('Found particles.')
  
  ! now shuffle particles and do the same
  do p = 1, krn_numProcs
    
    ! find smoothing kernel weights for each smoothed particle
    call krn_calculateWeights()

    ! send particles to the next processor and recv new particles
    call krn_moveParticles()
  end do
  call krn_timestamp('Calculated weights.')
  
  ! clean up the grid before we map to it, if necessary
  if (zero) call krn_zeroGridVar()
  
  ! map the smoothed particles onto the grid
  do p = 1, krn_numProcs
    
    ! map the smoothed particles onto the grid
    call krn_smoothParticles()

    ! send particles to the next processor and recv new particles
    ! don't need to send them back to the original processor
    if (p.ne.krn_numProcs) call krn_moveParticles()
  end do
  call krn_timestamp('Smoothed particles.')

  ! free up the temporary storage
  call krn_free()
  
  ! map the non-smoothed particles onto the grid
  if (partType.eq.PASSIVE_PART_TYPE) then
    call Particles_updateGridVar(MASS_PART_PROP, DENS_VAR, 1)
    call Particles_updateGridVar(VELX_PART_PROP, VELX_VAR, 1)
    call Particles_updateGridVar(VELY_PART_PROP, VELY_VAR, 1)
    call Particles_updateGridVar(VELZ_PART_PROP, VELZ_VAR, 1)
    call Particles_updateGridVar(MASS_SAVE_PART_PROP, TEMP_VAR, 1)
  else
    call Particles_updateGridVar(MASS_PART_PROP, PDE_VAR, 1)
  end if
  call krn_timestamp('Updated gridVar.')
  
  ! retype the smoothed particles and sort them again
  call krn_sortParticles()
  call krn_timestamp('Resorted particles.')
  
  ! divide by density for momentum and energy density to get velocities and temperature
  call krn_divideDensity()
  call krn_timestamp('Normalized blocks.')
  
  ! check that everything was conserved
  call krn_gridTotals()
  
  ! call Simulation_initBlock to initialize the other values and perform EOS
  call krn_initBlocks()
  call krn_timestamp('Done mapping.')

end subroutine

