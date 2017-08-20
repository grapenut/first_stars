

subroutine Simulation_init()

  use Simulation_data
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs, Driver_getSimTime
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY :  PhysicalConstants_get
  use Logfile_interface, ONLY : Logfile_stamp
  use Cosmology_data, ONLY : csm_scaleFactor
  use phch_interface, ONLY : phch_readFile
  use IO_interface, ONLY : IO_getScalar

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
#include "zoom.h"

  real :: scale
  character(len=256) :: str

  integer :: l, n
  real :: max_dens
  
  real, pointer, dimension(:,:) :: table
  integer :: rows, cols

  real :: time
 
  !! Setup processor/master for the simulation data
  call Driver_getMype(GLOBAL_COMM, sim_mype)
  call Driver_getNumProcs(GLOBAL_COMM, sim_numPEs)
  call Driver_getSimTime(time)

  if (sim_mype.eq.MASTER_PE) then
    sim_ismaster = .true.
  else
    sim_ismaster = .false.
  endif

  !! check for restarts
  call RuntimeParameters_get('restart', sim_restart)
 
  !! simulation domain size
  call RuntimeParameters_get('xmin',sim_xMin)
  call RuntimeParameters_get('xmax',sim_xMax)
  sim_box_size = sim_xMax - sim_xMin
  
  !! gas initialization
  call RuntimeParameters_get('initial_metallicity', sim_initial_metallicity)

  !! starburst/galaxy
  call RuntimeParameters_get('starburst_radius', sim_starburst_radius)
  call RuntimeParameters_get('starburst_refine_level', sim_starburst_refine_level)
  call RuntimeParameters_get('galaxy_refine_level', sim_galaxy_refine_level)
  call RuntimeParameters_get('galaxy_radius', sim_galaxy_radius)
  sim_galaxy_r2 = sim_galaxy_radius * sim_galaxy_radius
  sim_starburst_r2 = sim_starburst_radius * sim_starburst_radius

  !! supernovae
  call RuntimeParameters_get('nova_energy', sim_nova_energy)
  call RuntimeParameters_get('nova_yield', sim_nova_yield)
  call RuntimeParameters_get('nova_temp', sim_nova_temp)
  call RuntimeParameters_get('nova_particles', sim_nova_particles)
  call RuntimeParameters_get('zoom_level_offset', sim_zoom_level_offset)
  
  !! overdensity refinement 
  call RuntimeParameters_get('useOverdensity', sim_useOverdensity)
  call RuntimeParameters_get('useUnderdensity', sim_useUnderdensity)
  call RuntimeParameters_get('phi', sim_phi)
  call RuntimeParameters_get('refine_level_offset', sim_refine_level_offset)

  !! physical constants
  call PhysicalConstants_get('proton mass', m_H)
  call PhysicalConstants_get('electron mass', m_e)
  call PhysicalConstants_get('Boltzmann', k_B)
  !call RuntimeParameters_get('M_Sun', M_Sun)
  M_Sun = 1.99e33
  call RuntimeParameters_get('smallx', sim_smallx)
  call RuntimeParameters_get('HubbleConstant',sim_hubble)
  call RuntimeParameters_get('CosmologicalConstant',sim_lambda)
  call RuntimeParameters_get('OmegaMatter',sim_omegaMatter)
  call RuntimeParameters_get('OmegaBaryon',sim_omegaBaryon)
  call PhysicalConstants_get('Newton', sim_Newton)
  call PhysicalConstants_get('pi', sim_pi)

  !! constant limits
  call RuntimeParameters_get('smlrho',sim_smallrho)
  call RuntimeParameters_get('smallt',sim_smalltemp)
  call RuntimeParameters_get('smalle',sim_smalle)
  call RuntimeParameters_get('smallx',sim_smallx)
  call RuntimeParameters_get('smallp',sim_smallp)

  !! initialization data paths
  call RuntimeParameters_get('flash_data_path', str)
  write (sim_data_path, '(A)') trim(str)

  call RuntimeParameters_get('baryon_init_file', str)
  sim_baryon_init_file = trim(trim(sim_data_path)//"/"//trim(str))
  
  call RuntimeParameters_get('darkmatter_init_file', str)
  sim_darkmatter_init_file = trim(trim(sim_data_path)//"/"//trim(str))

  call RuntimeParameters_get('star_list_file', str)
  sim_star_list_file = trim(trim(sim_data_path)//"/"//trim(str))
  
  !! particle loading parameters
  call RuntimeParameters_get('pt_maxPerProc', sim_particle_max_limit)

  ! Compute some derived quantities that we will need.
  scale = csm_scaleFactor
  sim_zinitial = (1.0/scale) - 1.0

  ! Comoving critical density.
  sim_critden  = (3.*sim_hubble**2 / (8.*sim_pi*sim_Newton))
  
  ! Comoving total matter density.
  sim_meanden  = sim_omegaMatter * sim_critden
  
  ! Comoving baryonic density
  sim_bmeanden = sim_omegaBaryon * sim_critden
  sim_dmeanden = (sim_omegaMatter-sim_omegaBaryon)*sim_critden
  
  ! The above densities will not be correct if extrapolated to the present day
  ! because of the influence of dark energy. But in the range 10 < z < 1000
  ! the above approximations works extremely well.
  
  if (sim_mype.eq.MASTER_PE) then
    do l = 1, 24
      max_dens = 3.0 * sim_bmeanden * 2.0** ((l - sim_refine_level_offset) * 3.0 *  (1.0 + sim_phi) )
      print *,'SIM_INIT: lref=', l, '  max_dens=', max_dens/m_H*scale*scale*scale
    end do
  end if

  if (.not.sim_restart) then
    ! happens when we are making a new grid from scratch
    sim_darkmatter_loaded = .false.
    sim_baryons_loaded = .false.
    sim_init_refine_done = .false.
    sim_particle_refine_done = .false.
    
    sim_start_time_offset = time
  else
    ! happens most of the time, except for the first time running
    sim_darkmatter_loaded = .true.
    sim_baryons_loaded = .true.
    sim_init_refine_done = .true.
    sim_particle_refine_done = .true.
    
    call IO_getScalar("sim_start_time_offset", sim_start_time_offset)
  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! HII Region & Supernova init

  !! load star list
  call phch_readFile(sim_mype, sim_star_list_file, table, rows, cols)
  
  sim_num_stars = rows
  
  ! these will be read from the file each startup
  allocate(sim_star_mass(sim_num_stars))
  allocate(sim_star_lifetime(sim_num_stars))
  allocate(sim_zoom_time(sim_num_stars))

  ! get table values 
  do n=1, sim_num_stars
    sim_star_lifetime(n) = table(n,1) * 3.154e7
    sim_star_mass(n) = table(n,2)
  end do
  
  deallocate(table)
  
  ! add the simulation initial time at start to the lifetime
  sim_zoom_time = sim_star_lifetime + sim_start_time_offset

  ! these will be stored in the checkpoint on restart
  allocate(sim_nova_radius(sim_num_stars))
  allocate(sim_nova_r2(sim_num_stars))
  allocate(sim_nova_velocity(sim_num_stars))
  allocate(sim_nova_started(sim_num_stars))
  allocate(sim_nova_time(sim_num_stars))
  allocate(sim_zoom_status(sim_num_stars))
  allocate(sim_zoom_level(sim_num_stars))
  allocate(sim_star_coords(sim_num_stars, 3))
  
  !! initialize these properties, though they may be overwritten
  sim_nova_radius = 0.0
  sim_nova_r2 = 0.0
  sim_nova_velocity = 0.0
  sim_nova_started = 0
  sim_nova_time = sim_zoom_time
  sim_zoom_status = ZOOM_NOT_STARTED
  sim_zoom_level = sim_starburst_refine_level
  sim_star_coords = 0.0
  
  ! copy previously loaded values from IO_readUserArray
  if (sim_restart .and. sim_num_saved.gt.0) then
    if (sim_num_saved.gt.sim_num_stars) then
      n = sim_num_stars
    else
      n = sim_num_saved
    end if
  
    sim_nova_radius(1:n) = sim_saved_nova_radius(1:n)
    sim_nova_r2(1:n) = sim_nova_radius(1:n)*sim_nova_radius(1:n)
    sim_nova_velocity(1:n) = sim_saved_nova_velocity(1:n)
    sim_nova_started(1:n) = sim_saved_nova_started(1:n)
    sim_nova_time(1:n) = sim_saved_nova_time(1:n)
    sim_zoom_status(1:n) = sim_saved_zoom_status(1:n)
    sim_star_coords(1:n,1:3) = sim_saved_star_coords(1:n,1:3)

    deallocate(sim_saved_nova_velocity)
    deallocate(sim_saved_nova_radius)
    deallocate(sim_saved_nova_time)
    deallocate(sim_saved_nova_started)
    deallocate(sim_saved_zoom_status)
    deallocate(sim_saved_star_coords)
  end if
  
  sim_log2 = log(2.0)

  where (sim_nova_radius.ne.0.0) 
    sim_zoom_level = floor(1.0 - log(sim_nova_radius / sim_box_size) / sim_log2) + sim_zoom_level_offset
  end where
  
end subroutine Simulation_init



