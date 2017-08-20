
! calculate total quantities for particles and print a summary

subroutine krn_partTotals()

  use Kernel_data
  use Particles_data, ONLY : particles, pt_typeInfo

  implicit none

#include "Particles.h"
#include "Flash.h"  
#include "Flash_mpi.h"
#include "constants.h"

  integer :: p, ierr
  real :: momentum
  
  integer, dimension(NPART_TYPES) :: counts
  
  krn_local_part_mass = 0.0
  krn_local_part_temp = 0.0
  krn_local_part_px = 0.0
  krn_local_part_py = 0.0
  krn_local_part_pz = 0.0

  do p = krn_pstart, krn_pend
    krn_local_part_mass = krn_local_part_mass + particles(MASS_PART_PROP,p)
    krn_local_part_temp = krn_local_part_temp + particles(MASS_SAVE_PART_PROP,p)
    krn_local_part_px = krn_local_part_px + particles(VELX_PART_PROP,p)
    krn_local_part_py = krn_local_part_py + particles(VELY_PART_PROP,p)
    krn_local_part_pz = krn_local_part_pz + particles(VELZ_PART_PROP,p)
  end do
  
  call MPI_ALLREDUCE(pt_typeInfo(PART_LOCAL,:), counts, NPART_TYPES, FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  
  call MPI_ALLREDUCE(krn_local_part_mass, krn_global_part_mass, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(krn_local_part_temp, krn_global_part_temp, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(krn_local_part_px, krn_global_part_px, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(krn_local_part_py, krn_global_part_py, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(krn_local_part_pz, krn_global_part_pz, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  
  call MPI_ALLREDUCE(krn_pnum, krn_global_num_particles, 1, FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

  momentum = sqrt(krn_global_part_px*krn_global_part_px + krn_global_part_py*krn_global_part_py + krn_global_part_pz*krn_global_part_pz)
  
  if (krn_mype.eq.MASTER_PE) then
    print *,'KERNEL: total num particles', krn_global_num_particles
    print *,'KERNEL: particle type count', sum(counts), counts(:)
    print *,'KERNEL: total particle mass', krn_global_part_mass
    print *,'KERNEL: total particle temp', krn_global_part_temp
    print *,'KERNEL: total particle momentum', momentum
  end if
  
end subroutine

