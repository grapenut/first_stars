!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine PhotoChem_init()

  use PhotoChem_data
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, & 
                                          RuntimeParameters_set
  use PhysicalConstants_interface, ONLY :  PhysicalConstants_get
  use IO_interface, ONLY : IO_getScalar
  use Multispecies_interface, ONLY : Multispecies_getProperty

implicit none

#include "Flash.h"
#include "constants.h"
#include "Multispecies.h"

  integer :: i, n, bin
  real :: pifourn, num_nucleons, num_electrons

  !logical :: restart

  call Driver_getMype(MESH_COMM, phch_proc_id)   

  if (phch_proc_id.eq.MASTER_PE) then
    phch_ismaster = .true.
  else
    phch_ismaster = .false.
  endif

  ! physical constants
  call PhysicalConstants_get('proton mass', phch_proton_mass)
  call PhysicalConstants_get('electron mass', phch_electron_mass)
  call PhysicalConstants_get('newton', phch_newton)
  call PhysicalConstants_get('speed of light', phch_lightspeed)

  ! atomic masses
  phch_hydrogen_mass = 1.6733e-24
  phch_helium_mass = 4.0 * phch_hydrogen_mass
  phch_carbon_mass = 12.0 * phch_hydrogen_mass
  phch_oxygen_mass = 16.0 * phch_hydrogen_mass

  
  !call RuntimeParameters_get('restart', restart)

  call RuntimeParameters_get('usePhoto', phch_usePhoto)
  call RuntimeParameters_get('useChem', phch_useChem)
  
  call RuntimeParameters_get('smallx', phch_smallx)

  call RuntimeParameters_get('phch_num_sides', phch_nside)
  call RuntimeParameters_get('phch_bin_min_radius', phch_bin_min_physical)
  call RuntimeParameters_get('phch_bin_max_radius', phch_bin_max_physical)
  call RuntimeParameters_get('phch_max_octree_level', phch_max_octree_level)

  call RuntimeParameters_get('phch_default_source_x', phch_default_source_x)
  call RuntimeParameters_get('phch_default_source_y', phch_default_source_y)
  call RuntimeParameters_get('phch_default_source_z', phch_default_source_z)
  call RuntimeParameters_get('phch_default_photon_number', phch_default_photon_number)
  
  ! go ahead and use those defaults for now
  phch_source_x = phch_default_source_x
  phch_source_y = phch_default_source_y
  phch_source_z = phch_default_source_z
  phch_photon_number = phch_default_photon_number
  
  ! paths to text files containing our data tables
  call RuntimeParameters_get('flash_data_path', phch_data_path)
  call RuntimeParameters_get('photoeq_file', phch_photoeq_file)
  call RuntimeParameters_get('metal_cooling_file', phch_metal_cooling_file)

  ! photoeq table abundances
  call RuntimeParameters_get('helium_abundance', phch_helium_abundance)

  ! SPECIES PARTICLE MASSES
  n = 0
  do i = SPECIES_BEGIN, SPECIES_END
    call Multispecies_getProperty(i, A, num_nucleons)   
    call Multispecies_getProperty(i, E, num_electrons)

    if (i .eq. ELEC_SPEC) then
      num_nucleons = 0.0
      num_electrons = 1.0
    end if

    phch_particle_mass(i) = num_nucleons * phch_proton_mass + num_electrons * phch_electron_mass
    
    if (i.ne.Z_SPEC) then
      n = n + 1
      phch_spec_list(n) = i
    end if
  end do

  !! setup pixel and bin limits
  phch_npix = 12*phch_nside*phch_nside

  !! set some defaults for now
  phch_bin_min = phch_bin_min_physical
  phch_bin_max = phch_bin_max_physical
  phch_bin_min_squared = phch_bin_min * phch_bin_min
  phch_bin_min_cubed = phch_bin_min_squared * phch_bin_min
  
  phch_bin_max_over_min = phch_bin_max/phch_bin_min
  phch_bin_log_max_over_min = log10(phch_bin_max_over_min)

  !! automatically scale nbins so that pixel width is approximately
  !!  equal to bin width
  pifourn = PI / (4.0 * phch_nside)
  phch_nbin = floor(phch_bin_log_max_over_min / log10(-(pifourn+1.0)/(pifourn-1.0)))

  allocate(phch_bin_radius(phch_nbin+1))
  allocate(phch_bin_radius_physical(phch_nbin+1))
  allocate(phch_bin_r3(phch_nbin+1))
  allocate(phch_bin_r3_physical(phch_nbin+1))
  do i=1,phch_nbin+1
    ! radius of inner edge of bin
    phch_bin_radius_physical(i) = phch_bin_min * phch_bin_max_over_min**(real(i-1)/real(phch_nbin))
    phch_bin_r3_physical(i) = phch_bin_radius_physical(i)**3.0
  enddo
  phch_bin_radius(:) = phch_bin_radius_physical(:)
  phch_bin_r3(:) = phch_bin_r3_physical(:)

  phch_current_dens = 0.0
  phch_current_alpha = 0.0
  phch_current_invscale3 = 1.0
  phch_current_scale2 = 1.0
  phch_current_scale3 = 1.0
  phch_octree_level = 0
  phch_total_h_fraction = 0.0
  phch_hplus_fraction = 0.0
  
  allocate(phch_pixel_volume(phch_npix,phch_nbin))
  phch_pixel_volume = 0.0
  
  allocate(phch_pixel_number(phch_npix,phch_nbin))
  phch_pixel_number = 0.0

  allocate(phch_pixel_density(phch_npix,phch_nbin))
  phch_pixel_density = 0.0

  allocate(phch_pixel_lambda(phch_npix,phch_nbin))
  phch_pixel_lambda = 0.0

  allocate(phch_pixel_photons(phch_npix,phch_nbin))
  phch_pixel_photons = 0.0

  allocate(phch_pixel_hplus(phch_npix,phch_nbin))
  phch_pixel_hplus = 0.0

  allocate(phch_stromgren_radius(phch_npix))
  phch_stromgren_radius = 0.0

  allocate(phch_stromgren_r2(phch_npix))
  phch_stromgren_r2 = 0.0

  ! load photo equilibrium and metal cooling tables
  call phch_loadPhoto()
  call phch_loadMetal()

end subroutine PhotoChem_init

