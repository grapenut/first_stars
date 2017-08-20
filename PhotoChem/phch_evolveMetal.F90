
subroutine phch_evolveMetal(y, dens, eint, dt)

  use PhotoChem_data, ONLY : phch_metal_num, phch_dens_num, phch_temp_num, &
       phch_metal_list, phch_dens_list, phch_temp_list, phch_spec_list, phch_temp_min, &
       phch_temp_max, phch_iTemp, phch_useChem
  use phch_interface, ONLY : phch_interpMetal, phch_gasConstant, phch_electronEquilibrium
  use Multispecies_interface, ONLY : Multispecies_getSumInv

implicit none

#include "Flash.h"
#include "Multispecies.h"

  real, INTENT(INOUT) :: y(SPECIES_BEGIN:SPECIES_END+1), eint
  real, INTENT(IN) :: dens, dt

  real :: current_time

  real :: nH, metal, temperature
  
  integer :: ix, jx, kx
  
  real :: dens0, metal0, temp0
  
  real, dimension(SPECIES_BEGIN:SPECIES_END+1) :: y_eq
  
  real :: t_eq
  real :: gascon

  real,parameter :: safety_coeff = 0.1

  if (.not.phch_useChem) return

  current_time = 0.0
  
  ! number density of hydrogen atoms
  nH = y(H_SPEC) + y(HPLU_SPEC)
  
  ! metallicity has been rescaled from cloudy units to mass fraction, so this is no longer necessary
  metal = y(Z_SPEC)
  
  temperature = max(phch_temp_min, min(phch_temp_max,y(phch_iTemp)))

  ! table indices
  ix = max(1,min(phch_dens_num-1,maxloc(phch_dens_list, 1, phch_dens_list.le.nH)))
  jx = max(1,min(phch_metal_num-1,maxloc(phch_metal_list, 1, phch_metal_list.le.metal)))
  kx = max(1,min(phch_temp_num-1,maxloc(phch_temp_list, 1, phch_temp_list.le.temperature)))

  ! interpolation scalars
  dens0 = (nH - phch_dens_list(ix)) / (phch_dens_list(ix+1) - phch_dens_list(ix))
  metal0 = (metal - phch_metal_list(jx)) / (phch_metal_list(jx+1) - phch_metal_list(jx))
  temp0 = (temperature - phch_temp_list(kx)) / (phch_temp_list(kx+1) - phch_temp_list(kx))
  
  if (metal0.lt.0.0) metal0 = 0.0
  if (metal0.gt.1.0) metal0 = 1.0

  ! interpolate equilibrium values at the given coordinates
  call phch_interpMetal(y_eq, ix, jx, kx, dens0, metal0, temp0)

  ! rescale equilibrium abundances to local abundances
  y_eq(phch_spec_list) = y_eq(phch_spec_list) * nH
  
  ! recover the electron abundance by conservation (minus metal contribution)
  call phch_electronEquilibrium(y_eq)
  
  !! cooling subcycle loop
  do

    !!!!!!!!!!!!!!

    t_eq = safety_coeff * abs(eint / ( y_eq(phch_iTemp) * nH * nH / dens))

    if (current_time + t_eq.lt.dt) then
      current_time = current_time + t_eq
    else
      t_eq = dt - current_time
      current_time = dt
    end if

    ! remove the radiated energy from internal energy, recalculate temperature
    eint = eint - y_eq(phch_iTemp)*nH*nH*t_eq/dens

    call phch_gasConstant(y_eq, dens, gascon)
    
    y(phch_iTemp) = max(phch_temp_min, min(phch_temp_max, eint*gascon))
    temperature = y(phch_iTemp)

    !kx = max(1,min(phch_temp_num-1, 1 + floor( real(phch_temp_num-1) * log10(temperature / phch_temp_min) / log_temp_coeff) ))
    kx = max(1,min(phch_temp_num-1,maxloc(phch_temp_list, 1, phch_temp_list.le.temperature)))
    
    ! interpolation scalars
    temp0 = (temperature - phch_temp_list(kx)) / (phch_temp_list(kx+1) - phch_temp_list(kx))
    
    ! interpolate equilibrium values at the given coordinates
    call phch_interpMetal(y_eq, ix, jx, kx, dens0, metal0, temp0)

    ! rescale equilibrium abundances to local abundances
    y_eq(phch_spec_list) = y_eq(phch_spec_list) * nH
    
    ! recover the electron abundance by conservation (minus metal contribution)
    call phch_electronEquilibrium(y_eq)
    
    if (current_time.ge.dt) exit
    
  end do

  !!!  end of metal cooling  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine !! ch_evolveMetal

