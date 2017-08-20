module Kernel_data
  
  implicit none
  
  logical, save :: krn_debug = .false.
  
  ! mpi parameters
  integer, save :: krn_mype, krn_numProcs
  logical, save :: krn_evenProc
  
  ! grid boundaries
  real, save :: krn_xmin, krn_xmax, krn_dx
  real, save :: krn_ymin, krn_ymax, krn_dy
  real, save :: krn_zmin, krn_zmax, krn_dz
  logical, save :: krn_periodic = .false.
  
  ! small value limits
  real, save :: krn_smallrho, krn_smalltemp
  
  ! smoothed particle parameters
  real, dimension(:,:), allocatable, save :: krn_particles
  integer, save :: krn_num_particles
  integer, save :: krn_max_particles
  integer, save :: krn_total_particles
  
  ! particle shuffling buffer and parameters
  real, dimension(:,:), allocatable, save :: krn_recvBuf
  integer, save :: krn_num_recv

  ! smoothed particle properties
  integer, parameter :: krn_prop_x = 1
  integer, parameter :: krn_prop_y = 2
  integer, parameter :: krn_prop_z = 3
  integer, parameter :: krn_prop_mass = 4
  integer, parameter :: krn_prop_px = 5
  integer, parameter :: krn_prop_py = 6
  integer, parameter :: krn_prop_pz = 7
  integer, parameter :: krn_prop_temp = 8
  integer, parameter :: krn_prop_len = 9
  integer, parameter :: krn_prop_wgt = 10
  integer, parameter :: krn_num_props = 10
  integer, save :: krn_smoothType, krn_partType
  real, save :: krn_real_smoothType, krn_real_partType

  ! conservation checking values
  real, save :: krn_local_grid_mass, krn_local_part_mass, krn_global_grid_mass, krn_global_part_mass
  real, save :: krn_local_grid_px, krn_local_part_px, krn_global_grid_px, krn_global_part_px
  real, save :: krn_local_grid_py, krn_local_part_py, krn_global_grid_py, krn_global_part_py
  real, save :: krn_local_grid_pz, krn_local_part_pz, krn_global_grid_pz, krn_global_part_pz
  real, save :: krn_local_grid_temp, krn_local_part_temp, krn_global_grid_temp, krn_global_part_temp
  
  ! particles array quantities
  integer, save :: krn_global_num_particles, krn_pstart, krn_pnum, krn_pend

end module 
