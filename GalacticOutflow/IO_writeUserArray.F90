!!****if* source/Simulation/SimulationMain/StirTurb/IO_writeUserArray
!!
!!
!!  NAME
!!    IO_writeUserArray
!!
!!  SYNOPSIS
!!    call IO_writeUserArray()
!!
!!
!!  DESCRIPTION
!!
!!    This is the supplied interface for users to write out additional
!!    quantities to the checkpoint or plotfile.  This routine should be used
!!    for writing out various types of arrays.  If the array is a global quantity
!!    only the master processor needs to write out the data.  If it is a quantity 
!!    which is different on all processors then each processor must write out its 
!!    own section of the array. (For a serial IO implementation each processor would
!!    need to send its data to the master.)  The specific implementation is left up
!!    to the user.  
!!
!!    In each case the user should make a call to either 
!!    io_h5write_generic_int_arr (hdf5) or io_ncmpi_write_generic_iarr (pnetcdf) or
!!    io_h5write_generic_real_arr (hdf5) or io_ncmpi_write_generic_rarr (pnetcdf)
!!    depending on the io implementation.
!!  
!!  ARGUMENTS
!!    
!!
!!  NOTES 
!!
!!    This routine should NOT
!!    be used to write out grid scope data or to write out single scalar
!!    values.  To write out user defined grid scope variables the user should
!!    use the keyword 'GRIDVAR' to declare a grid scope variable in the Config
!!    files.  Then set the runtime parameters plot_grid_var_1, plot_grid_var_2,   
!!    to the name of the grid var to include them in the checkpoint files and
!!    plotfiles.
!!
!!    To write out single scalar quantities the use the IO_setScalar routine to
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

subroutine IO_writeUserArray ()
  
  use Simulation_data, ONLY : sim_num_stars, sim_nova_velocity, sim_star_coords, &
          sim_nova_radius, sim_nova_time, sim_nova_started, sim_zoom_status   

  use IO_data, ONLY : io_chkptFileID, io_globalMe

  implicit none

#include "constants.h"

  integer :: offset, datasetNameLen
  
  if (io_globalMe.ne.MASTER_PE) return
  
  offset = 0
  datasetNameLen = 17 !this makes things a lot easier because string handling from fortran to c is a big pain
  call io_h5write_generic_real_arr(io_globalMe, &
       io_chkptFileID, &
       sim_nova_velocity, &
       sim_num_stars, &
       sim_num_stars, &
       offset, &
       "sim_nova_velocity", &
       datasetNameLen)
  
  datasetNameLen = 15
  call io_h5write_generic_real_arr(io_globalMe, &
       io_chkptFileID, &
       sim_nova_radius, &
       sim_num_stars, &
       sim_num_stars, &
       offset, &
       "sim_nova_radius", &
       datasetNameLen)
  
  datasetNameLen = 13
  call io_h5write_generic_real_arr(io_globalMe, &
       io_chkptFileID, &
       sim_nova_time, &
       sim_num_stars, &
       sim_num_stars, &
       offset, &
       "sim_nova_time", &
       datasetNameLen)
  
  datasetNameLen = 16
  call io_h5write_generic_int_arr(io_globalMe, &
       io_chkptFileID, &
       sim_nova_started, &
       sim_num_stars, &
       sim_num_stars, &
       offset, &
       "sim_nova_started", &
       datasetNameLen)
  
  datasetNameLen = 15
  call io_h5write_generic_int_arr(io_globalMe, &
       io_chkptFileID, &
       sim_zoom_status, &
       sim_num_stars, &
       sim_num_stars, &
       offset, &
       "sim_zoom_status", &
       datasetNameLen)

  datasetNameLen = 15
  call io_h5write_generic_real_arr(io_globalMe, &
       io_chkptFileID, &
       sim_star_coords, &
       3*sim_num_stars, &
       3*sim_num_stars, &
       offset, &
       "sim_star_coords", &
       datasetNameLen)
  
end subroutine IO_writeUserArray
