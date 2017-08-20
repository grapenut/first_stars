
! shuffle particles from one processor to the next

subroutine krn_moveParticles()

  use Kernel_data

  implicit none

#include "constants.h"
#include "Flash_mpi.h"

  integer :: destProc, srcProc, ierr

  ! for every pair of processors 0 and 1, 1 and 2, 2 and 3, 3 and 4, etc.
  ! even numbered procs will send data to proc+1 first
  ! odd processors will receive data from proc-1
  ! then reverse, odd procs send to proc+1 while even receive from proc-1
  
  destProc = krn_mype+1
  if (destProc.gt.krn_numProcs-1) destProc = 0
  
  srcProc = krn_mype-1
  if (srcProc.lt.0) srcProc = krn_numProcs - 1
  
  krn_recvBuf(:,:) = 0.0
  krn_num_recv = 0
  if (krn_evenProc) then
    ! send first
    call MPI_SEND(krn_num_particles, 1, FLASH_INTEGER, destProc, 111, MPI_COMM_WORLD, ierr)
    if (krn_num_particles.gt.0) then
      call MPI_SEND(krn_particles, krn_num_particles*krn_num_props, FLASH_REAL, destProc, 222, MPI_COMM_WORLD, ierr)
    end if

    ! recv now
    call MPI_RECV(krn_num_recv, 1, FLASH_INTEGER, srcProc, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    if (krn_num_recv.gt.0) then
      call MPI_RECV(krn_recvBuf, krn_num_recv*krn_num_props, FLASH_REAL, srcProc, 222, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
  else
    ! recv first
    call MPI_RECV(krn_num_recv, 1, FLASH_INTEGER, srcProc, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    if (krn_num_recv.gt.0) then
      call MPI_RECV(krn_recvBuf, krn_num_recv*krn_num_props, FLASH_REAL, srcProc, 222, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    ! send now
    call MPI_SEND(krn_num_particles, 1, FLASH_INTEGER, destProc, 111, MPI_COMM_WORLD, ierr)
    if (krn_num_particles.gt.0) then
      call MPI_SEND(krn_particles, krn_num_particles*krn_num_props, FLASH_REAL, destProc, 222, MPI_COMM_WORLD, ierr)
    end if
  end if

  ! copy recvBuf to krn_particles
  krn_num_particles = krn_num_recv
  if (krn_num_recv.gt.0) then
    krn_particles(:,1:krn_num_recv) = krn_recvBuf(:,1:krn_num_recv)
  end if

  if (krn_debug .and. krn_mype.eq.MASTER_PE) then
    print *,'KERNEL: shuffled particles', krn_num_particles
  end if

end subroutine

