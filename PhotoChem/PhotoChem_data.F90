
Module PhotoChem_data

implicit none

#include "Flash.h"
#include "constants.h"

  !! Control switches
  logical, save :: phch_usePhoto, phch_useChem
  
  !! MPI data
  integer, save :: phch_proc_id
  logical, save :: phch_ismaster
  
  !! data paths
  character(len=256), save :: phch_data_path, phch_photoeq_file, phch_metal_cooling_file

  !! physical constants
  real, save :: phch_proton_mass, phch_electron_mass, phch_newton, phch_lightspeed, phch_smallx
  real, save :: phch_hydrogen_mass, phch_helium_mass, phch_carbon_mass, phch_oxygen_mass
  real, save :: phch_helium_abundance
  real, save, dimension(SPECIES_BEGIN:SPECIES_END) :: phch_particle_mass
  real, save, dimension(NSPECIES-1) :: phch_spec_list

  integer, parameter :: phch_iTemp = SPECIES_END+1
  integer, parameter :: phch_iThermal = Z_SPEC

  !! pixel mapping parameters
  integer, save :: phch_nside, phch_npix, phch_nbin
  real, save, allocatable, dimension(:,:) :: phch_pixel_volume, phch_pixel_number, phch_pixel_density, phch_pixel_lambda, phch_pixel_photons, phch_pixel_hplus
  real, save :: phch_binzero_volume, phch_binzero_number, phch_binzero_density, phch_binzero_lambda, phch_binzero_photons, phch_binzero_hplus
  real, save, allocatable, dimension(:) :: phch_stromgren_radius, phch_stromgren_r2

  real, save, allocatable, dimension(:) :: phch_bin_radius, phch_bin_radius_physical, phch_bin_r3, phch_bin_r3_physical
  real, save :: phch_current_dens, phch_current_alpha, phch_current_invscale3, phch_current_scale2, phch_current_scale3
  real, save :: phch_bin_min, phch_bin_max, phch_bin_min_physical, phch_bin_max_physical
  real, save :: phch_bin_log_max_over_min, phch_bin_max_over_min, phch_bin_min_cubed, phch_bin_min_squared
  real, save :: phch_total_h_fraction, phch_hplus_fraction, phch_num_dens_ceiling
  integer, save :: phch_max_octree_level, phch_octree_level

  !! photo source properties
  real, save :: phch_source_x, phch_source_y, phch_source_z, phch_photon_number

  !! default source parameters
  real, save :: phch_default_source_x, phch_default_source_y, phch_default_source_z, phch_default_photon_number
  
  !! photo equilibrium paramters
  integer, save :: phch_num_xi, phch_num_ximetal
  real, save :: phch_xi_min, phch_xi_max
  real, save :: phch_ximetal_min, phch_ximetal_max
  real, save, dimension(:), allocatable :: phch_xi_list, phch_ximetal_list
  real, save, dimension(:,:,:), allocatable :: phch_xi_table
  
  !! metal cooling parameters
  real, save :: phch_dens_min, phch_dens_max
  real, save :: phch_temp_min, phch_temp_max
  real, save :: phch_metal_min, phch_metal_max
  integer, save :: phch_dens_num, phch_temp_num, phch_metal_num
  real, save, dimension(:), allocatable :: phch_dens_list, phch_temp_list, phch_metal_list
  real, save, dimension(:,:,:,:), allocatable :: phch_cool_table  
  
  
  !! Miscellaneous

end Module PhotoChem_data

