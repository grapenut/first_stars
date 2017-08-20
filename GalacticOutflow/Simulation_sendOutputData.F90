!!****f* source/Simulation/Simulation_sendOutputData
!!
!! NAME
!!  Simulation_sendOutputData
!!
!! SYNOPSIS
!! 
!!  Simulation_sendOutputData()
!!  
!! DESCRIPTION 
!!
!! This routine sends the scalar variables owned by the Simulation unit
!! to the IO unit, to be written to a checkpoint file.
!!
!!
!!***

subroutine Simulation_sendOutputData()

  use Simulation_data, ONLY : sim_num_stars, sim_start_time_offset

  use IO_interface, ONLY : IO_setScalar

implicit none

  call IO_setScalar("sim_num_saved", sim_num_stars)
  call IO_setScalar("sim_start_time_offset", sim_start_time_offset)

end subroutine Simulation_sendOutputData

