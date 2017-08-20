
! sort particles and retype any previously smoothed particles (sorting again if needed)

subroutine krn_sortParticles()

  use Kernel_data
  use Particles_data, ONLY : pt_numLocal, pt_maxPerProc, particles, pt_typeInfo
  use pt_interface, ONLY : pt_updateTypeDS
  use Grid_interface, ONLY : Grid_sortParticles

  implicit none

#include "Flash.h"
#include "Particles.h"
  
  integer, dimension(MAXBLOCKS, NPART_TYPES) :: particlesPerBlk
  integer :: pstart, pend, pnum
  
  ! presort particles, shouldn't take much time if they are already sorted from before, but hey why not?
  call Grid_sortParticles(particles, NPART_PROPS, pt_numLocal, NPART_TYPES, pt_maxPerProc, particlesPerBlk, BLK_PART_PROP, TYPE_PART_PROP)
  call pt_updateTypeDS(particlesPerBlk)

  ! make sure all the particles are reclassified from SMOOTH to PASSIVE
  pstart = pt_typeInfo(PART_TYPE_BEGIN, krn_smoothType)
  pnum = pt_typeInfo(PART_LOCAL, krn_smoothType)
  pend = pstart + pnum - 1
  
  if (pnum.gt.0) then
    particles(TYPE_PART_PROP, pstart:pend) = krn_real_partType
  
    ! sort again now that we've retyped the SMOOTH particles
    call Grid_sortParticles(particles, NPART_PROPS, pt_numLocal, NPART_TYPES, pt_maxPerProc, particlesPerBlk, BLK_PART_PROP, TYPE_PART_PROP)
    call pt_updateTypeDS(particlesPerBlk)
  end if
  
  ! gather the new particle indices
  krn_pstart = pt_typeInfo(PART_TYPE_BEGIN, krn_partType)
  krn_pnum = pt_typeInfo(PART_LOCAL, krn_partType)
  krn_pend = krn_pstart + krn_pnum - 1

end subroutine

