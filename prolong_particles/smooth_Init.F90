

subroutine smooth_Init

  use smooth_Data
  !use tree, ONLY : lrefine, nodetype, newchild, nchild, nfaces, lnblocks, maxblocks_tr
  use tree, ONLY : maxblocks_tr
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  !use Particles_data, ONLY : particles
  use Driver_interface, ONLY : Driver_getMype  

  implicit none

#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"
  
  integer :: ierr, c, i

  if(smooth_initialized) return
  smooth_initialized = .true.
  
  call RuntimeParameters_get("lrefine_smooth_max", lrefine_smooth_max)
  call RuntimeParameters_get("lrefine_smooth_min", lrefine_smooth_min)
  call RuntimeParameters_get("smooth_radius_factor", smooth_radius_factor)
  call RuntimeParameters_get("SmoothParticles", SmoothParticles)
  call RuntimeParameters_get("max_num_smooth_part", max_num_smooth_part)
  call Driver_getMype(GLOBAL_COMM, smooth_mype)
  
  if(smooth_mype .eq. MASTER_PE) then
    print*,  "Initializing particle smoothing."
    print*, "SmoothParticles = ", SmoothParticles
    print*, "lrefine_smooth_min = ", lrefine_smooth_min
    print*, "lrefine_smooth_max = ", lrefine_smooth_max
    print*, "smooth_radius_factor = ", smooth_radius_factor
  endif

  
  if(.not. SmoothParticles) return
  
  nbuf_restrict = maxblocks_tr
  nbuf_prolong = maxblocks_tr/8
  
  allocate(send_prolong_data(NXB,NYB,NZB,nchild,nbuf_prolong), &
       stat=ierr)
  allocate(recv_prolong_data(NXB,NYB,NZB), stat=ierr)
  allocate(send_prolong_req(nchild*nbuf_prolong), stat=ierr)

  hg_restrict_n1 = NXB/2
  hg_restrict_n2 = NYB/2
  hg_restrict_n3 = NZB/2

  allocate(send_restrict_data(hg_restrict_n1,hg_restrict_n2,&
           hg_restrict_n3,nbuf_restrict), stat=ierr)
  allocate(recv_restrict_data(hg_restrict_n1, &
           hg_restrict_n2,hg_restrict_n3), stat=ierr)
  allocate(send_restrict_req(nbuf_restrict), stat=ierr)
  


  n1off(1) = 0
  n2off(1) = 0
  n3off(1) = 0
  n1off(2) = NXB/2
  n2off(2) = 0
  n3off(2) = 0
  n1off(3) = 0
  n2off(3) = NYB/2
  n3off(3) = 0
  n1off(4) = NXB/2
  n2off(4) = NYB/2
  n3off(4) = 0
  n1off(5) = 0
  n2off(5) = 0
  n3off(5) = NZB/2
  n1off(6) = NXB/2
  n2off(6) = 0
  n3off(6) = NZB/2
  n1off(7) = 0
  n2off(7) = NYB/2
  n3off(7) = NZB/2
  n1off(8) = NXB/2
  n2off(8) = NYB/2
  n3off(8) = NZB/2

  allocate(Px(-2:2,NXB,nchild), stat=ierr)
  allocate(Py(-2:2,NYB,nchild), stat=ierr)
  allocate(Pz(-2:2,NZB,nchild), stat=ierr)


  do c = 1, nchild
     do i = 1, NXB-1, 2
        Px(:,i,  c) = (/ -3./128.,  11./64., 1., -11./64.,  3./128. /)
        Px(:,i+1,c) = (/  3./128., -11./64., 1.,  11./64., -3./128. /)
     enddo
     do i = 1, NYB-1, 2
        Py(:,i,  c) = (/ -3./128.,  11./64., 1., -11./64.,  3./128. /)
        Py(:,i+1,c) = (/  3./128., -11./64., 1.,  11./64., -3./128. /)
     enddo
     do i = 1, NZB-1, 2
        Pz(:,i,  c) = (/ -3./128.,  11./64., 1., -11./64.,  3./128. /)
        Pz(:,i+1,c) = (/  3./128., -11./64., 1.,  11./64., -3./128. /)
     enddo
  enddo
  

  hg_ili = 1 + NGUARD
  hg_iui = NXB + NGUARD
  hg_jli = 1 + NGUARD
  hg_jui = NYB + NGUARD
  hg_kli = 1 + NGUARD
  hg_kui = NZB + NGUARD

  hg_ile = 1
  hg_iue = NXB + 2*NGUARD
  hg_jle = 1
  hg_jue = NYB + 2*NGUARD
  hg_kle = 1
  hg_kue = NZB + 2*NGUARD
  
  allocate(smooth_saveNodeType(maxblocks_tr),stat=ierr)
  allocate(smooth_saveNewChild(maxblocks_tr),stat=ierr)
  
  nmax = NXB*NYB*NZB
  
  blk_og      = 1
  blk_coarse  = 2
  proc_og     = 3
  proc_coarse = 4
  posx        = 5
  posy        = 6
  posz        = 7
  accx        = 8
  accy        = 9
  accz        = 10
  mass        = 11
  pid         = 12
  lref_dest   = 13
  
  num_smooth_part_props = 13
  
  allocate(particles_smooth(num_smooth_part_props, max_num_smooth_part))
  
  smooth_posAttrib(IAXIS) = posx
  smooth_posAttrib(JAXIS) = posy
  smooth_posAttrib(KAXIS) = posz
  
  num_smooth = 0
  !particles(SMOOTHTAG_PART_PROP,:) = 0.0

  
  return

end subroutine smooth_Init
