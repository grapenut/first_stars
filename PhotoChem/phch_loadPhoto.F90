subroutine phch_loadPhoto()

  use PhotoChem_data, ONLY : phch_num_xi, phch_xi_min, phch_xi_max, phch_xi_list, &
       phch_num_ximetal, phch_ximetal_min, phch_ximetal_max, phch_ximetal_list, &
       phch_xi_table, phch_data_path, phch_photoeq_file, &
       phch_iThermal, phch_iTemp, phch_smallx, phch_proc_id, phch_ismaster, phch_helium_abundance

  use phch_interface, ONLY : phch_findNumIndices, phch_makeLogIndices, phch_readFile
  
implicit none

#include "Flash.h"

  character(len=256) :: which
  real, pointer, dimension(:,:) :: table

  integer :: rows, cols, i, j, n
  real :: metal, flux
  
  real :: log_xi_coeff, log_ximetal_coeff

  integer, parameter :: iFlux = 1
  integer, parameter :: iMetal = 2
  integer, parameter :: iTemp = 3
  integer, parameter :: ie = 4
  integer, parameter :: iH = 5
  integer, parameter :: iHplus = 6
  integer, parameter :: iHe = 7
  integer, parameter :: iHeplus = 8
  integer, parameter :: iHeplusplus = 9
  integer, parameter :: iC = 10
  integer, parameter :: iCplus = 11
  integer, parameter :: iCplusplus = 12
  integer, parameter :: iThermal = 13
  integer, parameter :: iRecomb = 14
  integer, parameter :: iXmetal = 15
  
  ! get the raw data table from the data file
  which = trim(trim(phch_data_path)//"/"//trim(phch_photoeq_file))
  call phch_readFile(phch_proc_id, which, table, rows, cols)
  
  ! scan through and find the min/max and number of unique indices for each parameter
  call phch_findNumIndices(table(:,iFlux), rows, phch_num_xi, phch_xi_min, phch_xi_max)
  call phch_findNumIndices(table(:,iMetal), rows, phch_num_ximetal, phch_ximetal_min, phch_ximetal_max)
  
  ! allocate some arrays to hold the data table
  allocate(phch_xi_list(phch_num_xi))
  allocate(phch_ximetal_list(phch_num_ximetal))
  allocate(phch_xi_table(phch_num_xi, phch_num_ximetal, SPECIES_BEGIN:SPECIES_END+1))
  
  ! generate log spaced indices based on the above findings
  ! returns the values, not their logs
  call phch_makeLogIndices(phch_xi_list, phch_num_xi, phch_xi_min, phch_xi_max)
  call phch_makeLogIndices(phch_ximetal_list, phch_num_ximetal, phch_ximetal_min, phch_ximetal_max)
  
  ! rescale metal indices in mass fraction Z instead of relative abundance X
  !phch_ximetal_list(:) = phch_metal_scale * phch_ximetal_list(:) / (phch_hyhel_scale + phch_metal_scale * phch_ximetal_list(:))
  
  ! rebuild the metal index min/max
  ! use n to avoid errors that might occur by passing phch_num_ximetal as both INTENT IN and OUT separately
  !n = phch_num_ximetal
  !call phch_findNumIndices(phch_ximetal_list, n, phch_num_ximetal, phch_ximetal_min, phch_ximetal_max)
  
  phch_xi_table = phch_smallx
  log_xi_coeff = 1.0 / log10(phch_xi_max/phch_xi_min)
  log_ximetal_coeff = 1.0 / log10(phch_ximetal_max/phch_ximetal_min)

  if (phch_ismaster) then
    print*,'PHCH: flux = ', phch_xi_min, phch_xi_max, phch_num_xi
    print*,'PHCH: metal = ', phch_ximetal_min, phch_ximetal_max, phch_num_ximetal
  end if

  do n=1,rows
  
    ! the photoequilibrium parameters
    ! flux is really xi=F/nH where nH has been set to 1
    flux = table(n,iFlux)
    ! metal is the relative abundance here still
    metal = max(table(n,iMetal), 1.0e-99)
    !metal = phch_metal_scale * metal / (phch_hyhel_scale + phch_metal_scale * metal)
 
    ! these indices will be in the same order as the rescaled mass fractions
    i = max(1,min(phch_num_xi, floor( 1.1 + real(phch_num_xi-1) * log10(flux / phch_xi_min) * log_xi_coeff ) ))
    j = max(1,min(phch_num_ximetal, floor( 1.1 + real(phch_num_ximetal-1) * log10(metal / phch_ximetal_min) * log_ximetal_coeff ) ))
 
    ! for the purpose of units, consider these ion fractions to be multiplied by nH=1
    ! we do this normalization just in case it is not exactly 1, and to enforce helium abundance
    phch_xi_table(i, j, H_SPEC) = ( table(n,iH) ) / ( table(n,iH) + table(n,iHplus) )
    phch_xi_table(i, j, HPLU_SPEC) = ( table(n,iHplus) ) / ( table(n,iH) + table(n,iHplus) )
    
    phch_xi_table(i, j, HEL_SPEC) = table(n,iHe) * phch_helium_abundance         / ( table(n,iHe) + table(n,iHeplus) + table(n,iHeplusplus) )
    phch_xi_table(i, j, HEP_SPEC) = table(n,iHeplus) * phch_helium_abundance     / ( table(n,iHe) + table(n,iHeplus) + table(n,iHeplusplus) )
    phch_xi_table(i, j, HEPP_SPEC) = table(n,iHeplusplus) * phch_helium_abundance / ( table(n,iHe) + table(n,iHeplus) + table(n,iHeplusplus) )

    phch_xi_table(i, j, ELEC_SPEC) = phch_xi_table(i, j, HPLU_SPEC) + phch_xi_table(i, j, HEP_SPEC) + 2.0*phch_xi_table(i, j, HEPP_SPEC)

    phch_xi_table(i, j, phch_iTemp) = table(n,iTemp)
    phch_xi_table(i, j, phch_iThermal) = table(n,iThermal)
 
    !if (phch_ismaster) print '(A,I,E,I,E,E)','PHCHINIT: ', i, flux, j, metal, table(n,iTemp)
  
  end do

  deallocate(table)
  
  !if (phch_ismaster) then
  !  do n=1,phch_num_xi
  !    print *,'PHCHINIT: ', phch_xi_list(n), phch_xi_table(n,1,phch_iTemp)
  !  end do
  !end if

end subroutine

