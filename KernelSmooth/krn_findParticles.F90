
! locate particles with smooth_len > cell_size
! and copy them to temporary smoothing buffer

subroutine krn_findParticles()

  use Kernel_data
  use Kernel_interface, ONLY : krn_makeParticle
  use Particles_data, ONLY : particles, pt_numLocal
  use tree, ONLY : lrefine, lnblocks
  use Grid_data, ONLY : gr_delta

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"
  
  integer :: p, b, ierr, krn_count, total_count
  real :: x, y, z, mass, px, py, pz, temp, len
  
  krn_num_particles = 0
  
  ! find the local number of particles that need to be smoothed
  if (krn_pnum.gt.0) then
    do p = krn_pstart, krn_pend
      b = int(particles(BLK_PART_PROP,p))
      if (b < 0 .or. b > lnblocks) print *,'KERNEL: invalid block', b
      
      if (krn_partType.eq.DARK_PART_TYPE .or. particles(SMOOTHTAG_PART_PROP,p).gt.gr_delta(IAXIS,lrefine(b))) then
        krn_num_particles = krn_num_particles + 1
      end if
    end do
    !krn_num_particles = count(particles(SMOOTHTAG_PART_PROP,krn_pstart:krn_pend).gt.gr_delta(IAXIS,lrefine(int(particles(BLK_PART_PROP,krn_pstart:krn_pend)))))
  end if

  ! find the largest buffer and resize to match
  call MPI_ALLREDUCE(krn_num_particles, krn_max_particles, 1, FLASH_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(krn_num_particles, krn_total_particles, 1, FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(krn_pnum, krn_count, 1, FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(pt_numLocal, total_count, 1, FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (krn_mype.eq.MASTER_PE) then
    print *,'KERNEL: found particles ', krn_max_particles, krn_total_particles, krn_count, total_count
  end if
  
  ! allocate the buffer to the appropriate size
  allocate(krn_particles(krn_num_props, krn_max_particles))
  allocate(krn_recvBuf(krn_num_props, krn_max_particles))

  krn_particles(:,:) = 0.0
  
  ! nothing here to copy
  if (krn_num_particles.eq.0) return

  ! we'll reset this now and increment when we read in each particle
  krn_num_particles = 0
  
  ! copy smoothed particles to buffer
  do p = krn_pstart, krn_pend
  
    b = int(particles(BLK_PART_PROP,p))
    if (b < 0 .or. b > lnblocks) print *,'KERNEL: invalid block', b

    ! if smooth_len > delta, we need to smooth
    if (krn_partType.eq.DARK_PART_TYPE .or. particles(SMOOTHTAG_PART_PROP,p).gt.gr_delta(IAXIS,lrefine(b))) then
      x = particles(POSX_PART_PROP,p)
      y = particles(POSY_PART_PROP,p)
      z = particles(POSZ_PART_PROP,p)
      mass = particles(MASS_PART_PROP,p)
      px = particles(VELX_PART_PROP,p)
      py = particles(VELY_PART_PROP,p)
      pz = particles(VELZ_PART_PROP,p)
      temp = particles(MASS_SAVE_PART_PROP,p)
      len = particles(SMOOTHTAG_PART_PROP,p)
      
      call krn_makeParticle(x, y, z, mass, px, py, pz, temp, len)
      
      ! we can just change the particle type to avoid mapping later
      particles(TYPE_PART_PROP,p) = krn_real_smoothType

    end if
    
  end do

end subroutine

