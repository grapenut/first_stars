!!****if* source/Grid/GridMain/paramesh/gr_markVarThreshold
!!
!! NAME
!!  gr_markVarThreshold
!!
!! SYNOPSIS
!!  gr_markVarThreshold(  integer(in) :: Var
!!                         real(in)   :: var_th, 
!!                         integer(in) :: icmp, 
!!                         integer(in) :: lref )
!!
!! PURPOSE
!!  Refine all blocks for which a given variable (Var) exceeds or falls
!!  below some threshold (var_th).  The direction of the threshold is
!!  controlled by the parameter icmp.  Either blocks are brought
!!  up to a specific level of refinement or each block is refined once.
!!
!! ARGUMENTS
!!  Var -    the variable of interest
!!
!!  var_th  -     the limit on the variable
!! 
!!  icmp  -  icmp < 0  refine if the variable is less than var_th
!!         icmp >= 0 refine if the variable is greater then var_th
!! 
!!   lref -       If > 0, bring all qualifying blocks to this level of refinement.
!!
!!               If <= 0, refine qualifying blocks once.
!!
!! NOTES
!! 
!!  This routine has not yet been tested and should be used only as a guideline for
!!  a user's implementation.
!!
!!***

!!REORDER(5):unk

subroutine sim_markOverdensity()

  use tree, ONLY : refine, derefine, lrefine, lnblocks, nodetype, lrefine_min, stay
  use Grid_data, ONLY : gr_maxRefine
  use physicaldata, ONLY : unk
  
  use Simulation_data, ONLY : sim_bmeanden, sim_ismaster, sim_phi, sim_refine_level_offset, &
           sim_useOverdensity, sim_useUnderdensity
  use PhotoChem_data, ONLY : phch_source_x, phch_source_y, phch_source_z
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_getBlkPhysicalSize
  implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"

  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------

  ! Local data

  integer :: b, ii, jj, kk, count
  logical :: Grid_mark
  real :: factor_1, max_dens
  real :: dens_upper_threshold, dens_lower_threshold

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real :: redshift, oneplusz, oneplusz_neg2, oneplusz_cu
  
  real, allocatable, dimension(:) :: cell_x, cell_y, cell_z
  !  integer, dimension(2,MDIM)    :: blkLimits, blkLimitsGC
  integer                       :: sizeX, sizeY, sizeZ
  
  real			:: dist2, cx, cy, cz
  
  real :: refine_dist2, ic, jc, kc
  logical, save :: gcell = .true.
 
  !-------------------------------------------------------------------------------

  call Cosmology_getRedshift(redshift)
  
  oneplusz = redshift + 1.0
  oneplusz_neg2 = oneplusz**(-2.0)
  oneplusz_cu = oneplusz**3.0
  
  refine_dist2 = (0.5*3.0856775807e25)**2.0
  
  ic = phch_source_x
  jc = phch_source_y
  kc = phch_source_z
  
  
  if (sim_useOverdensity) then
  
    if (sim_ismaster) print *,'REFINE: Searching for Overdensity'
     
     ! loop over all leaf blocks: 
     do b = 1, lnblocks
        if (nodetype(b) .eq. LEAF) then
           
          ! Setup block limits with/without guard cells
          call Grid_getBlkIndexLimits(b,blkLimits,blkLimitsGC)
          !call Grid_getBlkPtr(b, solnData)

          ! Number of cells w/guardcells
          sizeX = blkLimitsGC(HIGH,IAXIS)
          sizeY = blkLimitsGC(HIGH,JAXIS)
          sizeZ = blkLimitsGC(HIGH,KAXIS)

          ! Allocate coordinate arrays
          allocate(cell_x(sizeX))
          allocate(cell_y(sizeY))
          allocate(cell_z(sizeZ))

          ! Initialize coordinate arrays
          cell_x = 0.0
          cell_y = 0.0
          cell_z = 0.0

          ! Setup coordinate arrays
          if (NDIM == 3) call Grid_getCellCoords(KAXIS, b, CENTER, gcell, cell_z, sizeZ)
          if (NDIM >= 2) call Grid_getCellCoords(JAXIS, b, CENTER, gcell, cell_y, sizeY)
          call Grid_getCellCoords(IAXIS, b, CENTER, gcell, cell_x, sizeX)
          
          cx = (cell_x(blkLimits(HIGH,IAXIS)) + cell_x(blkLimits(LOW,IAXIS)))*0.5 - ic
          cy = (cell_y(blkLimits(HIGH,JAXIS)) + cell_y(blkLimits(LOW,JAXIS)))*0.5 - jc
          cz = (cell_z(blkLimits(HIGH,KAXIS)) + cell_z(blkLimits(LOW,KAXIS)))*0.5 - kc
          
          dist2 = cx*cx + cy*cy + cz*cz
          
!          if ((oneplusz_neg2*dist2).lt.refine_dist2) then
           
           dens_upper_threshold = 3.0 * sim_bmeanden * 2.0**((lrefine(b) - sim_refine_level_offset) * 3.0 * (1.0 + sim_phi))
           
           ! Compute maximum density in block
           
           factor_1 = 0.0
           count = 0.0
           max_dens = 0.0
           
           do kk = NGUARD*K3D+1,NGUARD*K3D+NZB
              do jj = NGUARD*K2D+1,NGUARD*K2D+NYB
                 do ii = NGUARD+1,NGUARD+NXB 
                    
                    if(unk(DENS_VAR,ii,jj,kk,b) .gt. max_dens) then
                       max_dens = unk(DENS_VAR,ii,jj,kk,b)
                    endif
                    
                 end do
              end do
           end do
             
           Grid_mark = (max_dens*oneplusz_cu .gt. dens_upper_threshold)
            
           ! if refinement
           if (Grid_mark) then
              if(lrefine(b) .lt. gr_maxRefine) then
                 !print '(A,G,A,G,A,I)','REFINE: Overdense ', max_dens, ' > ', dens_upper_threshold, 'lref=', lrefine(b)+1
                 refine(b) = .true.
                 derefine(b) = .false.
              else if (lrefine(b) .ge. gr_maxRefine) then
                 derefine(b) = .false.
!                 stay(b) = .true.
              end if
           endif
           
!          endif !refine_dist2
           
!          if (lrefine(b) .ge. gr_maxRefine) refine(b) = .false.
           
          deallocate(cell_x)
          deallocate(cell_y)
          deallocate(cell_z)
           
        endif    ! leaf blocks
     enddo     ! blocks

  end if   
  
  if (sim_useUnderdensity) then
     
     !--- Doesn't work quite right ... probably shouldn't use.
     
     if (sim_ismaster) print *,'REFINE: Searching for Underdensity'
     
     ! loop over all leaf blocks: 
     do b = 1, lnblocks
        if (nodetype(b) .eq. LEAF) then
           
           
           dens_lower_threshold =  0.25 * sim_bmeanden * 2.0** ((lrefine(b) - sim_refine_level_offset) * 3.0 *  (1.0 + sim_phi) ) 
           
           
           ! compute maximum density in block
           factor_1 = 0.0
           count = 0.0
           max_dens = 0.0
           
           do kk = NGUARD*K3D+1,NGUARD*K3D+NZB
              do jj = NGUARD*K2D+1,NGUARD*K2D+NYB
                 do ii = NGUARD+1,NGUARD+NXB
                    factor_1 = factor_1 + unk(DENS_VAR,ii,jj,kk,b)
                    count = count + 1.0 
                    
                    if(unk(DENS_VAR,ii,jj,kk,b) .gt. max_dens) then
                       max_dens = unk(DENS_VAR,ii,jj,kk,b)
                    endif
                    
                 end do
              end do
           end do
           
           
           Grid_mark = (max_dens .lt. dens_lower_threshold)
           
           
           if(Grid_mark .and. (lrefine(b) .gt. lrefine_min) .and. (.not. refine(b)) .and. (.not. stay(b))) then
              derefine(b) =  .true.
              !print *,'REFINE: Underdense ', lrefine(b), ' to ', lrefine(b) - 1
           else
              derefine(b) = .false.
           end if
           
             
              
              !if(.not. Grid_mark) stay(b) = .true.
              
              !if (lrefine(b) .le. lrefine_min) derefine(b) = .false.
              
              
              
        endif   ! leaf blocks
     enddo     ! blocks
     
  end if
  

  return
end subroutine sim_markOverdensity

