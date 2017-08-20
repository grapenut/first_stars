!!****if* source/Simulation/SimulationMain/StirTurb/IO_readUserArray
!!
!!  NAME
!!    IO_readUserArray
!!
!!  SYNOPSIS
!!    call IO_readUserArray()
!!
!!
!!  DESCRIPTION
!!
!!    This is the supplied interface for users to read in additional
!!    quantities to the checkpoint or plotfile.  This routine should be used
!!    for reading in various types of arrays.  If the array is a global quantity
!!    only the master processor needs to read in the data.  If it is a quantity 
!!    which is different on all processors then each processor must read in its 
!!    own section of the array. (For a serial IO implementation each processor would
!!    need to send its data to the master.)  The specific implementation is left up
!!    to the user.  
!!
!!    In each case the user should make a call to either 
!!    io_h5read_generic_int_arr (hdf5) or io_ncmpi_read_generic_iarr (pnetcdf)  or
!!    io_h5read_generic_real_arr (hdf5) or io_ncmpi_read_generic_rarr (pnetcdf)
!!    depending on the io implementation.
!!  
!!  ARGUMENTS
!!    
!!   
!!
!!  NOTES 
!!
!!    This routine should NOT
!!    be used to read in grid scope data or to read in single scalar
!!    values.  To read in user defined grid scope variables the user should
!!    use the keyword 'GRIDVAR' to declare a grid scope variable in the Config
!!    files.  Then set the runtime parameters plot_grid_var_1, plot_grid_var_2,   
!!    to the name of the grid var to include them in the checkpoint files and
!!    plotfiles.
!!
!!    To read in single scalar quantities the use the IO_setScalar routine to
!!    add a scalar to the scalar output list.
!!
!!
!!  SEE ALSO
!!
!!    io_h5write_generic_int_arr
!!    io_h5read_generic_real_arr
!!    IO_setScalar
!!    
!!    For the pnetcdf implementation see
!!    io_ncmpi_write_generic_iarr
!!    io_ncmpi_read_generic_rarr
!!
!!***



subroutine IO_readUserArray ()
  
  use Simulation_data, ONLY : sim_num_saved, sim_saved_nova_velocity, sim_saved_star_coords, &
          sim_saved_nova_radius, sim_saved_nova_time, sim_saved_nova_started, sim_saved_zoom_status   

  use IO_data, ONLY : io_chkptFileID, io_globalMe
  
  use IO_interface, ONLY : IO_getScalar

  implicit none

#include "constants.h"
#include "Flash_mpi.h"

  integer :: offset, datasetNameLen, ierr
  
  if (io_globalMe.eq.MASTER_PE) print *, 'SAVE: reading checkpoint'
  call IO_getScalar('sim_num_saved', sim_num_saved)
  if (io_globalMe.eq.MASTER_PE) print *, 'SAVE: reading ', sim_num_saved, ' supernovae'
  
  if (sim_num_saved.eq.0) return
  
  allocate(sim_saved_nova_velocity(sim_num_saved))
  allocate(sim_saved_nova_radius(sim_num_saved))
  allocate(sim_saved_nova_time(sim_num_saved))
  allocate(sim_saved_nova_started(sim_num_saved))
  allocate(sim_saved_zoom_status(sim_num_saved))
  allocate(sim_saved_star_coords(sim_num_saved, 3))

  if (io_globalMe.eq.MASTER_PE) then
    print *,'SAVE: num stars saved =', sim_num_saved

    print *,'SAVE: reading sim_nova_velocity'
    offset = 0
    datasetNameLen = 17
    call io_h5read_generic_real_arr( &
         io_chkptFileID, &
         sim_saved_nova_velocity, &
         sim_num_saved, &
         sim_num_saved, &
         offset, &
         "sim_nova_velocity", &
         datasetNameLen)
    
    print *,'SAVE: reading sim_nova_radius'
    datasetNameLen = 15
    call io_h5read_generic_real_arr( &
         io_chkptFileID, &
         sim_saved_nova_radius, &
         sim_num_saved, &
         sim_num_saved, &
         offset, &
         "sim_nova_radius", &
         datasetNameLen)
    
    print *,'SAVE: reading sim_nova_time'
    datasetNameLen = 13
    call io_h5read_generic_real_arr( &
         io_chkptFileID, &
         sim_saved_nova_time, &
         sim_num_saved, &
         sim_num_saved, &
         offset, &
         "sim_nova_time", &
         datasetNameLen)
    
    print *,'SAVE: reading sim_nova_started'
    datasetNameLen = 16
    call io_h5read_generic_int_arr( &
         io_chkptFileID, &
         sim_saved_nova_started, &
         sim_num_saved, &
         sim_num_saved, &
         offset, &
         "sim_nova_started", &
         datasetNameLen)
    
    print *,'SAVE: reading sim_zoom_status'
    datasetNameLen = 15
    call io_h5read_generic_int_arr( &
         io_chkptFileID, &
         sim_saved_zoom_status, &
         sim_num_saved, &
         sim_num_saved, &
         offset, &
         "sim_zoom_status", &
         datasetNameLen)

    print *,'SAVE: reading sim_star_coords'
    datasetNameLen = 15
    call io_h5read_generic_real_arr( &
         io_chkptFileID, &
         sim_saved_star_coords, &
         3*sim_num_saved, &
         3*sim_num_saved, &
         offset, &
         "sim_star_coords", &
         datasetNameLen)
  
  end if
  
  call mpi_bcast(sim_saved_nova_radius, sim_num_saved, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call mpi_bcast(sim_saved_nova_velocity, sim_num_saved, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call mpi_bcast(sim_saved_nova_time, sim_num_saved, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call mpi_bcast(sim_saved_nova_started, sim_num_saved, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call mpi_bcast(sim_saved_zoom_status, sim_num_saved, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call mpi_bcast(sim_saved_star_coords, 3*sim_num_saved, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  
end subroutine IO_readUserArray
