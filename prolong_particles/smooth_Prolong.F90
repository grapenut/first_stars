
!! ARGUMENTS
!!
!!  level - the child level getting prolonged
!!  ito   - the grid variable to prolong into (on the parents)
!!  ifrom - the grid variable to prolong from (on the children)
!! mapMode - how to do the prolongation
!!           mapMode =  1 ---> umap3 prolongation
!!           mapMode =  2 ---> direct insertion


subroutine smooth_Prolong(level, ito, ifrom, mapMode)

  use Grid_data, ONLY : gr_meshComm, gr_meshMe, gr_intpol, gr_dirGeom, gr_smallx
  use Driver_interface, ONLY : Driver_abortFlash
  use tree, ONLY : lrefine, lnblocks,child,nchild,parent,nodetype
  use smooth_Data
  use physicaldata, ONLY: unk
  use workspace, ONLY : work
  use Grid_interface, ONLY : Grid_getCellCoords
 
 

  implicit none


#include "constants.h"
#include "Flash.h"  
#include "Flash_mpi.h"
  
  integer, intent(in)           :: ifrom, ito, level
  integer, intent(in)           :: mapMode

  integer                      :: b, c, p, h, i, j, k, ii, jj, kk, xx, yy, zz
  integer                      :: x, y, z, idest, ia, ja, ka, ib, jb, kb
  integer                      :: ierr1, ierr2, ierr3, blockID
  integer                      :: i1, i2, j1, j2, k1, k2, ichild
  integer                      :: ierr, nsent
  integer                      :: ioff,joff,koff
  real, pointer                :: solnData(:,:,:,:)
  integer                      :: status(MPI_STATUS_SIZE)
  integer                      :: send_status(MPI_STATUS_SIZE,nbuf_restrict)
  logical                      :: any_sent
  real                         :: prolongedSection(1:2,1:2,1:2)

  
  integer,parameter :: niver = 1
  integer,parameter :: lw = 4*(NXB+2*NGUARD+4) + K2D*(4*(NYB+2*NGUARD +4)+&
       5*(NXB+2*NGUARD+4)*NUNK_VARS) + K3D*(4*(NZB+2*NGUARD+4)+&
       5*(NXB+2*NGUARD+4)*(NYB+2*NGUARD+4)*NUNK_VARS)
  integer,parameter :: liw=  NXB+NYB+NZB+12+2*NGUARD
  integer,parameter :: ref_ratio = 2
  integer,parameter :: intpol_guard = 2
  integer,parameter :: twice_iguard = 2*intpol_guard
  integer,parameter :: lxiend = (NXB+2*NGUARD)/2 + twice_iguard
  integer,parameter :: linxu  = lxiend*niver
  integer,parameter :: lyiend = max(1,K2D*((NYB+2*NGUARD)/2 + twice_iguard))
  integer,parameter :: lziend = max(1,K3D*((NZB+2*NGUARD)/2 + twice_iguard))
  integer,parameter :: lxoend = NXB+2*NGUARD
  integer,parameter :: lonxu  = lxoend*niver
  integer,parameter :: lyoend = NYB+2*NGUARD
  integer,parameter :: lzoend = NZB+2*NGUARD
  integer :: offia,offib,offic,offoa,offob,offoc
  real, dimension(GRID_ILO_GC:GRID_IHI_GC) :: xo, xi, xb
  real, dimension(GRID_JLO_GC:GRID_JHI_GC) :: yo, yi, yb
  real, dimension(GRID_KLO_GC:GRID_KHI_GC) :: zo, zi, zb

  real, dimension(linxu,lyiend,lziend) :: pu
  real, dimension(lonxu,lyoend,lzoend) :: qu
 
  integer       :: inxu, onxu, imapcx
  real          :: pdx, pdy, pdz, cdx, cdy, cdz 
  integer, dimension(MDIM) :: geom
  real, dimension(lw)  :: wmap  
  integer, dimension(liw) :: iwmap  
  logical       :: conserved_var = .FALSE.
  integer       :: xoend, yoend, zoend, xiend, yiend, ziend

  integer       :: mvx_m, mvm_m, mvxu_m, mvxmu_m, mui_m
  integer       :: n_dens, nu_d, nud, iud, i_dens, nui, iu
  integer       :: ip_i_iiu, ip_ii_iiu, ip_jj_iiujj, ip_iiu_iiujj



#include "umap.h"
  COMMON  /amrmapi_p/                   &
       n_dens, nud(mui_m), iud(0:mui_m,mui_m),     &
       ip_i_iiu     (mvxu_m),          &
       ip_ii_iiu    (mvxu_m),          &
       ip_jj_iiujj  (mvxmu_m),         &
       ip_iiu_iiujj (mvxmu_m)

!==============================================================================


  !mapMode=1 ---> umap3 triquadratic mapping
  !mapMode=2 ---> direct insertion
  
  if(mapMode .ne. 1 .and. mapMode .ne. 2) then
     call Driver_abortFlash("smooth_Prolong: mapMode not supported")
  endif
  


  ia = 1+NGUARD
  ja = 1+NGUARD
  ka = 1+NGUARD
  ib = NXB+NGUARD
  jb = NYB+NGUARD
  kb = NZB+NGUARD  
  n_dens = 0

  nsent = 0
  h = 1
  
  do b = 1, lnblocks
     if((lrefine(b) .eq. level) .and. (smooth_saveNodeType(b) .ge. PARENT_BLK)) then
  



        
        ! direct insertion prolongation
        if(mapMode .eq. 2) then
           
           do c = 1, nchild
              
              send_prolong_data(:,:,:,c,h) = 0.0
              
              ! looping over children cells of block b
              do i = 1, NXB
                 do j = 1, NYB
                    do k = 1, NZB
                       ! prolongation by direct insertion:
                       
                       send_prolong_data(i,j,k,c,h) = unk(ifrom,NGUARD+n1off(c)+1+(i-1)/2,&
                            NGUARD+n2off(c)+1+(j-1)/2, & 
                            NGUARD+n3off(c)+1+(k-1)/2, b)
                       
                    enddo
                 enddo
              enddo
              
           enddo   ! block children
           
        endif
        






        if(mapMode .eq. 1) then 
           
           
           ! prolongation with umap3

        call Grid_getCellCoords(IAXIS,b,CENTER,.true.,xb,GRID_IHI_GC)
        call Grid_getCellCoords(JAXIS,b,CENTER,.true.,yb,GRID_JHI_GC)
        call Grid_getCellCoords(KAXIS,b,CENTER,.true.,zb,GRID_KHI_GC)
        
        pdx = (xb(2) - xb(1)) / 2.0
        cdx = pdx / 2.0
        ! assume uniform grid in all directions (pdx = pdy = pdz etc.):
        pdy = pdx
        pdz = pdx
        cdy = cdx
        cdz = cdx
        
        ! we also assume (but check) that xiend=yiend=ziend and xoend=yoend=zoend
        
        do c = 1, nchild
           
           send_prolong_data(:,:,:,c,h) = 0.0
           
           ! ths following is based on amr_prolong_gen_work1_fun
           ! and thus umap3
           
           ioff = n1off(c)
           joff = n2off(c)
           koff = n3off(c)
       
           offia  = (ia-1+NGUARD)/ref_ratio - intpol_guard
           offoa  = ia - 1
           xoend = ib-ia+1
           !xiend = (ib-ia+ref_ratio)/ref_ratio + min(1,mod(ib-ia,ref_ratio)*mod(ib+NGUARD,ref_ratio)) + twice_iguard
           xiend = NXB

           offib  = (ja-1+NGUARD)/ref_ratio - intpol_guard
           offob = ja - 1
           yoend = jb-ja+1
           !yiend = (jb-ja+ref_ratio)/ref_ratio + min(1,mod(jb-ja,ref_ratio)*mod(jb+NGUARD,ref_ratio)) + twice_iguard
           yiend = NYB
           
           offic  = (ka-1+NGUARD)/ref_ratio - intpol_guard
           offoc = ka - 1
           zoend = kb-ka+1
           !ziend = (kb-ka+ref_ratio)/ref_ratio + min(1,mod(kb-ka,ref_ratio)*mod(kb+NGUARD,ref_ratio)) + twice_iguard
           ziend = NZB
           
           ! (x,y,z)i = coords of parent starting two GC layers out
           ! (x,y,z)o = coords of interior child cells
           
           xi(1) = xb(1+intpol_guard+ioff)
           yi(1) = yb(1+intpol_guard+joff)
           zi(1) = zb(1+intpol_guard+koff)
           
           if(xiend .ne. yiend .or. xiend .ne. ziend) &
                call Driver_abortFlash("xiend, yiend, ziend aren't equal!")
           
           do i = 2, xiend
              xi(i) = xi(i-1) + 2.0 * pdx
              yi(i) = yi(i-1) + 2.0 * pdy
              zi(i) = zi(i-1) + 2.0 * pdz
           enddo
           

          
           xo(1) = xi(1) + twice_iguard * pdx - cdx
           yo(1) = yi(1) + twice_iguard * pdy - cdy
           zo(1) = zi(1) + twice_iguard * pdz - cdz
           
           if(xoend .ne. yoend .or. xoend .ne. zoend) &
                call Driver_abortFlash("xoend, yoend, zoend aren't equal!")
                      

           do i = 2, xoend
              xo(i) = xo(i-1) + 2.0 * cdx
              yo(i) = yo(i-1) + 2.0 * cdy
              zo(i) = zo(i-1) + 2.0 * cdz
           enddo

           
           inxu=niver*xiend
           onxu=niver*xoend  
           pu = 0.0
           geom=gr_dirGeom
           imapcx=1
           
           do k = 1,ziend
              k1 = (offic+koff)+k
              do j = 1,yiend
                 j1 = (offib+joff)+j
                 do i = 1,xiend
                    i1 = offia+ioff+i
                    pu(i,j,k) = unk(ifrom,i1,j1,k1,b)
                 end do
              end do
           end do
           
           !print*, 'parameters=', lxiend, lyiend, lziend, xiend, yiend, ziend
           !print*, 'more', linxu, inxu, pdx, pdy, pdz
           
           
           call umap3 (lxiend,lyiend,lziend,xiend,yiend,ziend,&
                &      linxu,inxu,xi,pdx,yi,pdy,zi,pdz,pu,&
                &      lxoend,lyoend,lzoend,xoend,yoend,zoend,&
                &      lonxu,onxu,xo,cdx,yo,cdy,zo,cdz,qu,&
                &      niver,gr_intpol,imapcx,&
                &      geom(1),geom(2),geom(3),ref_ratio,ref_ratio,ref_ratio,&
                &      wmap,lw,iwmap,liw,ierr,conserved_var,gr_smallx)
           


           do k = 1, zoend
              k1 = offoc + k
              do j = 1, yoend
                 j1 = offob + j
                 do i = 1, xoend
                    i1 = offoa + i
                    
                    send_prolong_data(i,j,k,c,h) = qu(i,j,k)
                    
     
                 enddo
              enddo
           enddo

           
        enddo           

           
     endif
        
        
      
           
        ! Now, loop over this block's children. If child is on processor, 
        ! set its values directly. If it is off-proc, send the data to the 
        ! owning processor
        
        
        
        any_sent = .false.
        
        do c = 1, nchild
           if(child(2,c,b) .eq. gr_meshMe) then ! local child
              
              blockID = child(1,c,b)
              
              unk(ito,NGUARD+1:NGUARD+NXB, & 
                   NGUARD+1:NGUARD+NYB, & 
                   NGUARD+1:NGUARD+NZB,blockID) = &
                   unk(ito,NGUARD+1:NGUARD+NXB, & 
                   NGUARD+1:NGUARD+NYB, & 
                   NGUARD+1:NGUARD+NZB,blockID) + & 
                   send_prolong_data(1:NXB,1:NYB,1:NZB,c,h)
!!$              
!!$              unk(ito,NGUARD+1:NGUARD+NXB, & 
!!$                   NGUARD+1:NGUARD+NYB, & 
!!$                   NGUARD+1:NGUARD+NZB,blockID) = &
!!$                   send_prolong_data(1:NXB,1:NYB,1:NZB,c,h)
              
           else  ! non-local child
              
              any_sent = .true.
              nsent = nsent + 1
              call mpi_issend(send_prolong_data(1,1,1,c,h), nmax, &
                   FLASH_REAL, child(2,c,b), child(1,c,b), &
                   gr_meshComm, send_prolong_req(nsent), ierr)
              
           endif
           
        enddo
        

        if(any_sent) h = h + 1
        if(h .gt. nbuf_prolong) call Driver_abortFlash("Buffer space exceeded in smooth_Prolong")
        
        
     end if ! lrefine = level?
  end do ! blocks
  
!  print*, "smooth_Prolong, h = ", h
  
  
  ! now child blocks receive messages from parents
  
  do b = 1, lnblocks
     if ((lrefine(b) .eq. level+1) .and. (parent(2,b) .ne. gr_meshMe)) then
        
        call mpi_recv(recv_prolong_data(1,1,1), nmax, &
             FLASH_REAL, parent(2,b), b, gr_meshComm, & 
             status, ierr)
        
        unk(ito, NGUARD+1:NGUARD+NXB, &
             NGUARD+1:NGUARD+NYB, & 
             NGUARD+1:NGUARD+NZB,b) = &
             unk(ito, NGUARD+1:NGUARD+NXB, &
             NGUARD+1:NGUARD+NYB, & 
             NGUARD+1:NGUARD+NZB,b) + & 
             recv_prolong_data(1:NXB,1:NYB,1:NZB)
!!$        
!!$        unk(ito, NGUARD+1:NGUARD+NXB, &
!!$             NGUARD+1:NGUARD+NYB, & 
!!$             NGUARD+1:NGUARD+NZB,b) = &
!!$             recv_prolong_data(1:NXB,1:NYB,1:NZB)
        
           
        endif  ! appropriate block
     enddo  ! child blocks
  
  
  call mpi_waitall(nsent, send_prolong_req, send_status, ierr)
  


  return

end subroutine smooth_Prolong





