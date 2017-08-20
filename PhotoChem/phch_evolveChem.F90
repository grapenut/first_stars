
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine phch_evolveChem(blockCount, blockList, dt)
  
  use PhotoChem_data, only : phch_useChem, phch_ismaster, phch_iTemp, phch_particle_mass, &
       phch_xi_min, phch_spec_list, phch_usePhoto, phch_xi_max, phch_smallx
  use phch_interface, ONLY : phch_evolvePhoto, phch_evolveMetal
  use Cosmology_interface, only : Cosmology_getRedshift
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
       Grid_releaseBlkPtr
  use Eos_interface, ONLY : Eos_wrapped

implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: blockCount
  integer, intent(in), dimension(blockCount)  :: blocklist
  real, intent(in) :: dt

  integer :: blk, blockID 
  integer :: i, j, k ! zone indices
 
  integer,dimension(2,MDIM)::blkLimits,blkLimitsGC
  real,pointer, dimension(:,:,:,: ) :: solnData
  
  real :: redshift, scale, scale2, invscale2, invscale3

  real, dimension(SPECIES_BEGIN:SPECIES_END+1) :: y

  real :: dens, sum_mass_fractions, eint, metal, photoeq, flux, nH

  real :: chem_timestep

  character (len=40) :: time_string

  !real :: jeans_delta_sq, jeans_numb_sq, cs2, temp_floor

  !real :: newton
  !real :: pmass, emass

  if (.not.(phch_useChem.or.phch_usePhoto)) return

  call current_date_time(time_string)
  if (phch_ismaster) print '(3A)', "PHCH: ", trim(time_string), " Evolving chemistry"

  call Cosmology_getRedshift(redshift)
  scale = 1.0 / (redshift + 1.0)
  scale2 = scale * scale
  invscale2 = 1.0 / scale2
  invscale3 = 1.0 / (scale * scale * scale)

  do blk = 1, blockCount
    blockID = blockList(blk)

    call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
    call Grid_getBlkPtr(blockID, solnData)

    do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
      do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

          !! chemical timestep for one zone
          
          !! basic properties in physical units
          dens = solnData(DENS_VAR,i,j,k) * invscale3
          eint = solnData(EINT_VAR,i,j,k) * scale2
          metal = solnData(Z_SPEC,i,j,k)
         
          !! convert mass fractions into number densities
          y(phch_spec_list) = dens * solnData(phch_spec_list, i, j, k) / phch_particle_mass(phch_spec_list)
          y(Z_SPEC) = metal
          y(phch_iTemp) = solnData(TEMP_VAR, i, j, k) * scale2
          
          !do m = SPECIES_BEGIN, SPECIES_END
          !  if (m.ne.Z_SPEC) y(m) = dens * solnData(m, i, j, k) / phch_particle_mass(m)
          !end do
          
          ! xi = F / n (cloudy has U = F / n / c but we multiply that out in the data table
          nH = y(H_SPEC) + y(HPLU_SPEC)
          flux = solnData(FLUX_VAR,i,j,k)
          
          ! let's disperse dense gas by blasting it with radiation
          if (nH.gt.1.0e3) then
            photoeq = phch_xi_max
          else
            photoeq = flux / nH
          end if
          
          ! Jeans floor
          !cs2 = solnData(PRES_VAR,i,j,k) / solnData(DENS_VAR,i,j,k) * scale
          !jeans_numb_sq = PI * cs2 / phch_newton / dens / jeans_delta_sq
          !temp_floor = temp * (jeans_ref_sq / jeans_numb_sq)

          if (photoeq.ge.phch_xi_min) then
            !if (photoeq.le.phch_xi_min) print *,'PHCHDEBUG: low photoeq ', photoeq, y
            call phch_evolvePhoto(y, photoeq, dens, eint, dt)
          else if (phch_useChem) then
            call phch_evolveMetal(y, dens, eint, dt)
          end if
          
          !if (chem_timestep .eq. 0.0) chem_timestep = 1.0e50
          solnData(CHEM_VAR, i, j, k) = 1.0e50 !chem_timestep
          
          !! convert number densities back to unitless mass fractions in solnData
          solnData(phch_spec_list, i, j, k) = y(phch_spec_list) * phch_particle_mass(phch_spec_list) / dens
          solnData(TEMP_VAR, i, j, k) = y(phch_iTemp) * invscale2
          
          !do m = SPECIES_BEGIN, SPECIES_END
          !  if (m.ne.Z_SPEC) solnData(m, i, j, k) = y(m) * phch_particle_mass(m) / dens
          !end do

          !! enforce crude mass conservation
          sum_mass_fractions = sum(solnData(SPECIES_BEGIN:SPECIES_END, i, j, k))
          solnData(SPECIES_BEGIN:SPECIES_END, i, j, k) = solnData(SPECIES_BEGIN:SPECIES_END, i, j, k) / sum_mass_fractions
          
        enddo
      enddo
    enddo

    call Grid_releaseBlkPtr(blockID, solnData)
    call Eos_wrapped(MODE_DENS_TEMP, blkLimits, blockID)
  end do

  call current_date_time(time_string)
  if (phch_ismaster) print '(3A)', "PHCH: ", trim(time_string), " Done with chemistry"

  return

end subroutine phch_evolveChem







