
! calculate total quantities for the grid and print a summary

subroutine krn_gridTotals()

  use Kernel_data
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPtr, &
       Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_getDeltas

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  
  integer, dimension(MAXBLOCKS) :: blockList
  integer :: blockCount, b, blockID
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(MDIM) :: delta
  integer :: i, j, k, ierr
  real :: cell_vol, cell_mass, momentum, partmom
  
  call Grid_getListOfBlocks(LEAF, blockList, blockCount)
  
  krn_local_grid_mass = 0.0
  krn_local_grid_temp = 0.0
  krn_local_grid_px = 0.0
  krn_local_grid_py = 0.0
  krn_local_grid_pz = 0.0
  
  do b = 1, blockCount
    blockID = blockList(b)
    call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
    call Grid_getBlkPtr(blockID, solnVec, CENTER)
    call Grid_getDeltas(blockID, delta)
    cell_vol = delta(IAXIS)*delta(JAXIS)*delta(KAXIS)
    
    do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
      do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
          cell_mass = solnVec(DENS_VAR,i,j,k) * cell_vol
          krn_local_grid_mass = krn_local_grid_mass + cell_mass
          krn_local_grid_temp = krn_local_grid_temp + cell_mass * solnVec(TEMP_VAR,i,j,k)
          krn_local_grid_px = krn_local_grid_px + cell_mass * solnVec(VELX_VAR,i,j,k)
          krn_local_grid_py = krn_local_grid_py + cell_mass * solnVec(VELY_VAR,i,j,k)
          krn_local_grid_pz = krn_local_grid_pz + cell_mass * solnVec(VELZ_VAR,i,j,k)
        end do
      end do
    end do
    call Grid_releaseBlkPtr(blockID, solnVec, CENTER)
  end do

  call MPI_ALLREDUCE(krn_local_grid_mass, krn_global_grid_mass, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(krn_local_grid_temp, krn_global_grid_temp, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(krn_local_grid_px, krn_global_grid_px, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(krn_local_grid_py, krn_global_grid_py, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(krn_local_grid_pz, krn_global_grid_pz, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  
  momentum = sqrt(krn_global_grid_px*krn_global_grid_px + krn_global_grid_py*krn_global_grid_py + krn_global_grid_pz*krn_global_grid_pz)
  partmom = sqrt(krn_global_part_px*krn_global_part_px + krn_global_part_py*krn_global_part_py + krn_global_part_pz*krn_global_part_pz)
  
  if (krn_mype.eq.MASTER_PE) then
    print *,'KERNEL: total grid mass', krn_global_grid_mass, 1.0 - krn_global_part_mass/krn_global_grid_mass
    if (krn_partType .eq. PASSIVE_PART_TYPE) then
      print *,'KERNEL: total grid temp', krn_global_grid_temp, 1.0 - krn_global_part_temp/krn_global_grid_temp
      print *,'KERNEL: total grid momentum', momentum, 1.0 - partmom/momentum
    end if
  end if
  
end subroutine

