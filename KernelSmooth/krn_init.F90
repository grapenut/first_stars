
! initialize constant parameters

subroutine krn_init(partType)

  use Kernel_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs

  implicit none
  
  integer, intent(in) :: partType

#include "constants.h"
#include "Flash.h"

  call Driver_getMype(MESH_COMM, krn_mype)
  call Driver_getNumProcs(MESH_COMM, krn_numProcs)
  krn_evenProc = mod(krn_mype,2).eq.0

  ! check that shuffle algorithm will work
  ! needs an even number of procs for alternating send/recv pairs
  if (krn_mype.eq.MASTER_PE .and. mod(krn_numProcs,2).ne.0) then
    print *,"KERNEL: There are an odd number of processors! krn_moveParticles() won't work."
  end if
  
  call RuntimeParameters_get('xmin', krn_xmin)
  call RuntimeParameters_get('xmax', krn_xmax)
  call RuntimeParameters_get('ymin', krn_ymin)
  call RuntimeParameters_get('ymax', krn_ymax)
  call RuntimeParameters_get('zmin', krn_zmin)
  call RuntimeParameters_get('zmax', krn_zmax)
  krn_dx = krn_xmax - krn_xmin
  krn_dy = krn_ymax - krn_ymin
  krn_dz = krn_zmax - krn_zmin
  
  call RuntimeParameters_get('smlrho', krn_smallrho)
  call RuntimeParameters_get('smallt', krn_smalltemp)
  
  krn_smoothType = SMOOTH_PART_TYPE
  krn_real_smoothType = real(SMOOTH_PART_TYPE)

  krn_partType = partType
  krn_real_partType = real(partType)
  
  krn_num_particles = 0
  krn_max_particles = 100000
  
end subroutine

