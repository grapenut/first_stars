
! divide density averaged quantities by the total cell density after all mapping is done

subroutine krn_divideDensity()

  use Kernel_data
  use PhotoChem_data, ONLY : phch_default_source_x, phch_default_source_y, phch_default_source_z
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPtr, Grid_getDeltas, Grid_getBlkBoundBox, &
       Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_getBlkIDFromPos, Grid_getPointData

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  
  integer, dimension(MAXBLOCKS) :: blockList
  integer :: blockCount, b, blockID
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  real, dimension(MDIM) :: pos
  integer, dimension(MDIM) :: intPos, guard
  real, dimension(MDIM) :: deltas, hfunc
  real, dimension(LOW:HIGH,MDIM) :: bndBox
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: i, j, k, empty, krn_global_empty, ierr
  
  integer :: centerBlk, centerProc
  real :: vx, vy, vz
  
  real :: cmbtemp
  
  cmbtemp = 2.7*(1.0+10.0)
  
  call Grid_getListOfBlocks(LEAF, blockList, blockCount)

  empty = 0
  do b = 1, blockCount
    blockID = blockList(b)
    call Grid_getBlkPtr(blockID, solnVec, CENTER)
    call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
    do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
      do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
          if (solnVec(DENS_VAR,i,j,k).gt.krn_smallrho) then
            solnVec(VELX_VAR,i,j,k)=solnVec(VELX_VAR,i,j,k)/solnVec(DENS_VAR,i,j,k)
            solnVec(VELY_VAR,i,j,k)=solnVec(VELY_VAR,i,j,k)/solnVec(DENS_VAR,i,j,k)
            solnVec(VELZ_VAR,i,j,k)=solnVec(VELZ_VAR,i,j,k)/solnVec(DENS_VAR,i,j,k)
            solnVec(TEMP_VAR,i,j,k)=solnVec(TEMP_VAR,i,j,k)/solnVec(DENS_VAR,i,j,k)
          else
            empty = empty + 1
            solnVec(DENS_VAR,i,j,k) = krn_smallrho
            solnVec(VELX_VAR,i,j,k) = 0.0
            solnVec(VELY_VAR,i,j,k) = 0.0
            solnVec(VELZ_VAR,i,j,k) = 0.0
            solnVec(TEMP_VAR,i,j,k) = krn_smalltemp
          end if
          
          !if (solnVec(TEMP_VAR,i,j,k).lt.cmbtemp) then
          !  solnVec(TEMP_VAR,i,j,k) = cmbtemp
          !end if
        end do
      end do
    end do
    call Grid_releaseBlkPtr(blockID, solnVec, CENTER)
  end do
  
  call MPI_ALLREDUCE(empty, krn_global_empty, 1, FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  
  if (krn_mype.eq.MASTER_PE .and. krn_global_empty.gt.0) then
    print *,'KERNEL: found empty cells', krn_global_empty
  end if
  
  ! get the velocity of the cell at the source coordinates (center of the domain)
  ! subtract out that velocity from the entire grid
  
  ! first get the block and proc containing the source point
  pos(1) = phch_default_source_x
  pos(2) = phch_default_source_y
  pos(3) = phch_default_source_z
  call Grid_getBlkIDFromPos(pos, centerBlk, centerProc, MPI_COMM_WORLD)

  ! extract the velocity
  if (centerProc.eq.krn_mype) then
    call Grid_getDeltas(centerBlk, deltas)
    call Grid_getBlkBoundBox(centerBlk, bndBox)
    call Grid_getBlkIndexLimits(centerBlk, blkLimits, blkLimitsGC)

    !! Get the integer indices of the cell containing the chosen coords
    !! using hfunc as the temporary storage for offset from the block bdry
    guard=blkLimits(LOW,:)-blkLimitsGC(LOW,:)
    hfunc(1:NDIM) = (pos(:) - bndBox(LOW,1:NDIM)) / deltas(1:NDIM)
    intPos(1:NDIM)= floor(hfunc(1:NDIM)) + 1 + guard(1:NDIM)
  
    call Grid_getPointData(centerBlk, CENTER, VELX_VAR, EXTERIOR, intPos, vx)
    call Grid_getPointData(centerBlk, CENTER, VELY_VAR, EXTERIOR, intPos, vy)
    call Grid_getPointData(centerBlk, CENTER, VELZ_VAR, EXTERIOR, intPos, vz)
  
    pos(1) = vx
    pos(2) = vy
    pos(3) = vz
  end if
  
  ! broadcast the central velocity to all procs
  call MPI_BCAST(pos, 3, FLASH_REAL, centerProc, MPI_COMM_WORLD, ierr)

  ! remove the central velocity from all cells
  do b = 1, blockCount
    blockID = blockList(b)
    call Grid_getBlkPtr(blockID, solnVec, CENTER)
    call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
    do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
      do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
          solnVec(VELX_VAR,i,j,k) = solnVec(VELX_VAR,i,j,k) - pos(1)
          solnVec(VELY_VAR,i,j,k) = solnVec(VELY_VAR,i,j,k) - pos(2)
          solnVec(VELZ_VAR,i,j,k) = solnVec(VELZ_VAR,i,j,k) - pos(3)
        end do
      end do
    end do
    call Grid_releaseBlkPtr(blockID, solnVec, CENTER)
  end do
  
  if (krn_mype.eq.MASTER_PE) then
    print *,'KERNEL: shifted velocities by ', pos(:)
  end if
  
end subroutine

