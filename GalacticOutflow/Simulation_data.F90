
module Simulation_data

  implicit none

#include "constants.h"
#include "Flash.h"
  
  !! particle initialization
  character(256), save :: sim_data_path, sim_baryon_init_file, sim_darkmatter_init_file
  logical, save :: sim_baryons_loaded, sim_darkmatter_loaded, sim_init_refine_done, sim_particle_refine_done
  integer, save :: sim_particle_max_limit
  
  !! starting time offset
  real, save :: sim_start_time_offset
  
  !! gas initialization
  real, save :: sim_initial_metallicity

  !! overdensity refinement
  logical, save :: sim_useOverdensity, sim_useUnderdensity
  real, save :: sim_phi
  real, save :: sim_refine_level_offset

  !! galaxy refinment
  real, save :: sim_galaxy_radius, sim_galaxy_r2
  integer, save :: sim_galaxy_refine_level

  !! star list
  character(256), save :: sim_star_list_file
  integer, save :: sim_num_stars
  real, save :: sim_starburst_radius, sim_starburst_r2
  integer, save :: sim_starburst_refine_level
  real, save, dimension(:), allocatable :: sim_star_mass, sim_star_lifetime
  real, save, dimension(:,:), allocatable :: sim_star_coords
  
  !! supernova properties
  real, save :: sim_nova_energy, sim_nova_temp, sim_nova_yield, sim_nova_particles
  real, save, dimension(:), allocatable :: sim_nova_velocity, sim_nova_radius, sim_nova_r2

  !! zoom control variables
  integer, save, dimension(:), allocatable :: sim_nova_started, sim_zoom_level
  real, save, dimension(:), allocatable :: sim_zoom_time, sim_nova_time
  integer, save, dimension(:), allocatable :: sim_zoom_status
  real, save :: sim_derefine_radius_base, sim_box_size, sim_block_usage_threshold
  integer, save :: sim_derefine_min, sim_star_lref, sim_star_lref_offset, sim_derefine_step
  integer, save :: sim_zoom_level_offset
  
  !! saved variables loaded from a checkpoint prior to Simulation_init being called
  integer, save :: sim_num_saved
  real, save, dimension(:), allocatable :: sim_saved_nova_velocity, sim_saved_nova_radius, sim_saved_nova_time
  integer, save, dimension(:), allocatable :: sim_saved_nova_started, sim_saved_zoom_status
  real, save, dimension(:,:), allocatable :: sim_saved_star_coords

  !! Some physical constants
  real, save :: m_H, m_e, k_B, M_Sun, sim_Newton, sim_pi, sim_log2
  real, save :: sim_smallx, sim_smalltemp, sim_smallRho, sim_smallE, sim_smallP
  real, save :: sim_meanden, sim_hubble, sim_lambda, sim_zinitial
  real, save :: sim_omegaMatter, sim_omegaBaryon, sim_bmeanden, sim_dmeanden, sim_critden
  
  
  !! Misc variables
  real, save :: sim_xMax, sim_xMin
  integer, save :: sim_mype, sim_numPEs
  logical, save :: sim_ismaster, sim_restart
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module Simulation_data
