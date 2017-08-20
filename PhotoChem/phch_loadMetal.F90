

subroutine phch_loadMetal()

  use PhotoChem_data, ONLY : phch_dens_min, phch_dens_max, phch_dens_num, phch_metal_max, &
       phch_metal_min, phch_metal_num, phch_temp_max, phch_temp_min, phch_temp_num, phch_proc_id, &
       phch_dens_list, phch_metal_list, phch_temp_list, phch_cool_table, phch_data_path, &
       phch_smallx, phch_helium_abundance, phch_metal_cooling_file, phch_ismaster, phch_iTemp
  
  use phch_interface, ONLY : phch_findNumIndices, phch_makeLogIndices, phch_readFile
  
implicit none

#include "Flash.h"

  character(len=256) :: which
  real, pointer, dimension(:,:) :: table

  integer :: rows, cols, i, j, k, n

  real :: dens, temp, metal
  real :: log_dens_coeff, log_temp_coeff, log_metal_coeff

  integer, parameter :: iDens = 1
  integer, parameter :: iMetal = 2
  integer, parameter :: iTemp = 3
  integer, parameter :: ie = 4
  integer, parameter :: iCool = 5
  integer, parameter :: iHeat = 6
  integer, parameter :: iH = 7
  integer, parameter :: iHplus = 8
  integer, parameter :: iHe = 9
  integer, parameter :: iHeplus = 10
  integer, parameter :: iHeplusplus = 11
  integer, parameter :: iC = 12
  integer, parameter :: iCplus = 13
  integer, parameter :: iCplusplus = 14

  ! get the raw data table from the data file
  which = trim(trim(phch_data_path)//"/"//trim(phch_metal_cooling_file))
  call phch_readFile(phch_proc_id, which, table, rows, cols)
  
  ! scan through and find the min/max and number of unique indices for each parameter
  call phch_findNumIndices(table(:,iDens), rows, phch_dens_num, phch_dens_min, phch_dens_max)
  call phch_findNumIndices(table(:,iTemp), rows, phch_temp_num, phch_temp_min, phch_temp_max)
  call phch_findNumIndices(table(:,iMetal), rows, phch_metal_num, phch_metal_min, phch_metal_max)

  ! allocate some arrays to hold the data table
  allocate(phch_dens_list(phch_dens_num))
  allocate(phch_temp_list(phch_temp_num))
  allocate(phch_metal_list(phch_metal_num))
  allocate(phch_cool_table(phch_dens_num, phch_metal_num, phch_temp_num, SPECIES_BEGIN:SPECIES_END+1))
  
  ! generate log spaced indices based on the above findings
  ! returns the values, not their logs
  call phch_makeLogIndices(phch_dens_list, phch_dens_num, phch_dens_min, phch_dens_max)
  call phch_makeLogIndices(phch_temp_list, phch_temp_num, phch_temp_min, phch_temp_max)
  call phch_makeLogIndices(phch_metal_list, phch_metal_num, phch_metal_min, phch_metal_max)
  
  ! rescale metal indices in mass fraction Z instead of relative abundance X
  !phch_metal_list(:) = phch_metal_scale * phch_metal_list(:) / (phch_hyhel_scale + phch_metal_scale * phch_metal_list(:))

  ! rebuild the metal index min/max
  ! use n to avoid errors that might occur by passing phch_num_ximetal as both INTENT IN and OUT separately
  !n = phch_metal_num
  !call phch_findNumIndices(phch_metal_list, n, phch_metal_num, phch_metal_min, phch_metal_max)
  
  phch_cool_table(:,:,:,:) = phch_smallx  
  log_dens_coeff = 1.0 / log10(phch_dens_max/phch_dens_min)
  log_metal_coeff = 1.0 / log10(phch_metal_max/phch_metal_min)
  log_temp_coeff = 1.0 / log10(phch_temp_max/phch_temp_min)
  
  if (phch_ismaster) then
    print*,'PHCH: dens = ', phch_dens_min, phch_dens_max, phch_dens_num
    print*,'PHCH: metal = ', phch_metal_min, phch_metal_max, phch_metal_num
    print*,'PHCH: temp = ', phch_temp_min, phch_temp_max, phch_temp_num
  end if

  do n=1,rows
    dens = table(n,iDens)
    metal = max(table(n,iMetal), 1.0e-99)
    !metal = phch_metal_scale * metal / (phch_hyhel_scale + phch_metal_scale * metal)
    temp = table(n,iTemp)

    i = max(1,min(phch_dens_num, floor( 1.1 + real(phch_dens_num-1) * log10(dens / phch_dens_min) * log_dens_coeff ) ))
    j = max(1,min(phch_metal_num, floor( 1.1 + real(phch_metal_num-1) * log10(metal / phch_metal_min) * log_metal_coeff ) ))
    k = max(1,min(phch_temp_num, floor( 1.1 + real(phch_temp_num-1) * log10(temp / phch_temp_min) * log_temp_coeff ) ))

    phch_cool_table(i,j,k,H_SPEC) = dens * table(n,iH) / (table(n,iH) + table(n,iHplus))
    phch_cool_table(i,j,k,HPLU_SPEC) = dens * table(n,iHplus) / (table(n,iH) + table(n,iHplus))

    phch_cool_table(i,j,k,HEL_SPEC) = phch_helium_abundance * dens * table(n,iHe) / (table(n,iHe) + table(n,iHeplus) + table(n,iHeplusplus))
    phch_cool_table(i,j,k,HEP_SPEC) = phch_helium_abundance * dens * table(n,iHeplus) / (table(n,iHe) + table(n,iHeplus) + table(n,iHeplusplus))
    phch_cool_table(i,j,k,HEPP_SPEC) = phch_helium_abundance * dens * table(n,iHeplusplus) / (table(n,iHe) + table(n,iHeplus) + table(n,iHeplusplus))

    phch_cool_table(i,j,k,ELEC_SPEC) = max(phch_smallx,table(n,ie) - (phch_cool_table(i,j,k,HPLU_SPEC) + phch_cool_table(i,j,k,HEP_SPEC) + 2.0*phch_cool_table(i,j,k,HEPP_SPEC)))
    
    phch_cool_table(i,j,k,phch_iTemp) = table(n,iCool) - table(n,iHeat)

    !if (phch_ismaster) print '(A,I,E,I,E,I,E)','PHCHINIT: ', i, dens, j, metal, k, temp
  
  end do
  
  deallocate(table)

end subroutine

