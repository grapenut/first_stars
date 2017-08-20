
! 2nd pass sums particle density on cells based on the kernel weighting

subroutine krn_smoothParticles()

  use Kernel_data
  use Kernel_interface, ONLY : krn_checkBounds, krn_getCellWeight
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getCellCoords, &
        Grid_getBlkIndexLimits, Grid_getDeltas, Grid_getBlkPtr, &
        Grid_releaseBlkPtr, Grid_getBlkBoundBox

  implicit none

#include "constants.h"
#include "Flash_mpi.h"
#include "Flash.h"

  real, dimension(LOW:HIGH, MDIM) :: bndBox
  integer, dimension(MAXBLOCKS) :: blockList
  integer :: blockCount, b, blockID
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(MDIM) :: delta
  integer :: sizeX, sizeY, sizeZ, istat
  real, dimension(:), allocatable :: x, y, z
  logical :: gcell=.true.
  real :: cell_vol
  integer :: i, j, k, p
  real :: posx, posy, posz, len, twolen, wgt, norm, factor
  logical :: overlap
  integer :: count
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec

  call Grid_getListOfBlocks(LEAF, blockList, blockCount)
  
  do b = 1, blockCount
  
    blockID = blockList(b)

    call Grid_getBlkBoundBox(blockID, bndBox)
    
    call Grid_getDeltas(blockID,delta)
    cell_vol = delta(IAXIS)*delta(JAXIS)*delta(KAXIS)
    
    call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

    sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
    allocate(x(sizeX), stat=istat)

    sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
    allocate(y(sizeY), stat=istat)

    sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1
    allocate(z(sizeZ), stat=istat)

    call Grid_getCellCoords(KAXIS, blockID, CENTER, gcell, z, sizeZ)
    call Grid_getCellCoords(JAXIS, blockID, CENTER, gcell, y, sizeY)
    call Grid_getCellCoords(IAXIS, blockID, CENTER, gcell, x, sizeX)
 
    call Grid_getBlkPtr(blockID, solnVec, CENTER)

    count = 0
    do p = 1, krn_num_particles
    
      posx = krn_particles(krn_prop_x,p)
      posy = krn_particles(krn_prop_y,p)
      posz = krn_particles(krn_prop_z,p)
      len = krn_particles(krn_prop_len,p)
      twolen = 2.0 * len
      overlap = .false.
      
      call krn_checkBounds(posx, posy, posz, twolen, bndBox, overlap, .true.)
      
      if (.not.overlap) then
        cycle
      end if
      
      count = count + 1

      norm = 1.0/krn_particles(krn_prop_wgt,p)
      
      do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
          do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

            wgt = 0.0
            call krn_getCellWeight(posx, posy, posz, x(i), y(j), z(k), twolen, wgt, .true.)
            
            if (wgt.gt.0.0) then
              factor = wgt*norm
              
              if (krn_partType.eq.PASSIVE_PART_TYPE) then
                solnVec(DENS_VAR,i,j,k) = solnVec(DENS_VAR,i,j,k) + krn_particles(krn_prop_mass,p)*factor
                solnVec(TEMP_VAR,i,j,k) = solnVec(TEMP_VAR,i,j,k) + krn_particles(krn_prop_temp,p)*factor
                solnVec(VELX_VAR,i,j,k) = solnVec(VELX_VAR,i,j,k) + krn_particles(krn_prop_px,p)*factor
                solnVec(VELY_VAR,i,j,k) = solnVec(VELY_VAR,i,j,k) + krn_particles(krn_prop_py,p)*factor
                solnVec(VELZ_VAR,i,j,k) = solnVec(VELZ_VAR,i,j,k) + krn_particles(krn_prop_pz,p)*factor
              else
                solnVec(PDE_VAR,i,j,k) = solnVec(PDE_VAR,i,j,k) + krn_particles(krn_prop_mass,p)*factor
              end if
            end if
            
          end do
        end do
      end do
    end do
    
    deallocate(x)
    deallocate(y)
    deallocate(z)
    
    call Grid_releaseBlkPtr(blockID, solnVec, CENTER)

    if (krn_mype.eq.MASTER_PE .and. krn_debug) then
      if (count.gt.0) then
        print *, 'KERNEL: checked block', blockID, blockCount, count
      else
        print *, 'KERNEL: skipped block', blockID, blockCount
      end if
    end if

  end do
  
  if (krn_debug .and. krn_mype.eq.MASTER_PE) then
    print *, 'KERNEL: done 2nd pass'
  end if
  
end subroutine

