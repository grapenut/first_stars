! Initializes baryons - if they exist
subroutine Simulation_initBlock (blockId)

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords, & 
                             Grid_putPointData, &
                             Grid_getPointData
  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_interface, ONLY : Eos, Eos_wrapped
  use Multispecies_interface, ONLY : Multispecies_getProperty
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Particles_data, ONLY : useParticles
  use Simulation_data
  use PhotoChem_data, ONLY : phch_helium_abundance

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Multispecies.h"
#include "Particles.h"

  integer,intent(IN) ::  blockId

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, parameter :: particleTypes = 1
  integer, dimension(MDIM) :: axis

  real, dimension(SPECIES_END) :: massFrac, numDens
  
  real :: particleMass, particleMass1
  real, save :: emass, pmass, abhe, abd
  
  integer :: i, j, k, n
  real :: rhodm, rho  
  logical, save :: first_call = .true.
  
   if(first_call) then
     call PhysicalConstants_get("proton mass", pmass)
     call PhysicalConstants_get("electron mass", emass)

  
     first_call = .false.
   end if

   ! Get the coordinate information for the current block
   call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

   do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
      do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
         do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
            
            axis(IAXIS) = i
            axis(JAXIS) = j
            axis(KAXIS) = k
            
            call Grid_getPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
            
            if (rho.ne.rho.or.rho.le.0.0) then
              print *, 'SIMERROR: density is whack', rho
              rho = sim_smallrho
            end if
            !call Grid_getPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, vx)
            !call Grid_getPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vy)
            !call Grid_getPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, vz)
            !call Grid_getPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, temp_zone)
            
            if (NSPECIES == 1) then 
               massFrac(1) = 1.0
               call Grid_putPointData(blockId, CENTER, SPECIES_BEGIN, EXTERIOR, axis, & 
                    massFrac(1))
            else if (NSPECIES > 1) then
               
               ! Initial number densities of species...
               ! (references from Anninos &  Norman (1996), Abel, Peebles...)
               
               numDens(H_SPEC)          = rho / (m_H * (1.0 + 4.0 * phch_helium_abundance))
               numDens(HPLU_SPEC)       = numDens(H_SPEC) * sim_smallx

               numDens(HEL_SPEC)        = numDens(H_SPEC) * phch_helium_abundance
               numDens(HEP_SPEC)        = numDens(H_SPEC) * sim_smallx
               numDens(HEPP_SPEC)       = numDens(H_SPEC) * sim_smallx

               numDens(ELEC_SPEC)	= numDens(HPLU_SPEC) + numDens(HEP_SPEC) + 2.0*numDens(HEPP_SPEC)
               numDens(Z_SPEC)		= 0.0
               ! set all the initial mass fractions:
               
               do n=SPECIES_BEGIN, SPECIES_END 
                  
                  call Multispecies_getProperty(n, A, particleMass)
                  call Multispecies_getProperty(n, E, particleMass1)
                  
                  if(n.eq.ELEC_SPEC) then
                     massFrac(n) = emass * particleMass1 * numDens(n) / rho
                  else
                     massFrac(n) = ((pmass * particleMass) + (emass * particleMass1)) * numDens(n) / rho
                  end if
                  
               end do
               
               ! normalize mass fractions and set the initial metallicity
               massFrac(Z_SPEC) = 0.0
               massFrac = massFrac / (sum(massFrac(:)) + sim_initial_metallicity)
               massFrac(Z_SPEC) = sim_initial_metallicity
               
               do n=SPECIES_BEGIN, SPECIES_END
                  call Grid_putPointData(blockId, CENTER, n, EXTERIOR, axis, massFrac(n))
               end do
               
            endif
            
            
            
            if (useParticles) then
               rhodm = sim_smallRho
               call Grid_putPointData(blockId, CENTER, PDEN_VAR, EXTERIOR, axis, rhodm)
            endif
            
            
         enddo
      enddo
   enddo

   ! No use or eos_wrapped when dealing with multispecies
   !call Eos_wrapped(MODE_DENS_TEMP, blklimitsGC, blockId)

   return
end subroutine Simulation_initBlock
