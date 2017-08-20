!!****if* source/Grid/GridSolvers/Multigrid/gr_hgRestrict
!!
!! NAME
!!  gr_hgRestrict
!!
!! SYNOPSIS
!!
!!  gr_hgRestrict(integer, intent(in) :: level,
!!                integer, intent(in) :: ito,
!!                integer, intent(in) :: ifrom)
!!
!! DESCRIPTION
!!  
!!  Restrict the interior data from blocks on a given level to their parents.
!!  This is used for defining the gravity source and residual on all levels.
!!  It works by zeroing the parent blocks, and then averaging 2^DIM child
!!  cells into a single parent cell, handling on and off-processor blocks
!!  separately.
!!
!! ARGUMENTS
!!
!!  level - the child level being restricted from
!!  ito   - the grid variable to restrict into (on the parents)
!!  ifrom - the grid variable to restrict from (on the children)
!!
!! NOTES
!!  
!!  This routine is significantly faster than the alternatives like
!!  gr_restrictTree
!!
!! SEE ALSO
!!
!! gr_hgProlong, gr_restrictTree
!!
!!***

!==============================================================================

!!REORDER(5): unk
!!REORDER(4): solnData

subroutine smooth_Restrict(level, ito, ifrom)

  use Grid_data, ONLY : gr_meshComm, gr_meshMe
  use Driver_interface, ONLY : Driver_abortFlash
  use tree, ONLY : lrefine, lnblocks,child,nchild,parent,nodetype
  use physicaldata, ONLY: unk
  use paramesh_interfaces, only : amr_restrict_unk_genorder
  use paramesh_dimensions, ONLY : nvar
  use smooth_Data

  implicit none

#include "constants.h"
#include "Flash.h"  
#include "Flash_mpi.h"
  
  integer, intent(in)          :: ifrom, ito, level
  
  integer                      :: b, c, p, h, i, j, k, ii, jj, kk
  integer                      :: kc, jc, ic
  integer                      :: ierr1, ierr2, ierr3
  integer                      :: i1, i2, j1, j2, k1, k2, ichild
  integer                      :: ierr, nsent
  real, pointer                :: solnData(:,:,:,:)
  integer                      :: status(MPI_STATUS_SIZE)
  integer                      :: send_status(MPI_STATUS_SIZE,nbuf_restrict)
  integer                      :: order = 1
  real                         :: dataout(nvar,1:2*NGUARD+NXB,1:2*NGUARD+NYB,1:2*NGUARD+NZB)
  
  !=======================================================================


  nsent = 0
  h = 1
  do b = 1, lnblocks
     if (lrefine(b) == level) then

        ! First, compute restricted values for this block.
!!$       
           do k = 1, hg_restrict_n3
              k1 = hg_kli + 2*k - 2
              k2 = k1 + 1
              do j = 1, hg_restrict_n2
                 j1 = hg_jli + 2*j - 2
                 j2 = j1 + 1
                 do i = 1, hg_restrict_n1
                    i1 = hg_ili + 2*i - 2
                    i2 = i1 + 1
                    send_restrict_data(i,j,k,h) = 0.125*(unk(ifrom,i1,j1,k1,b) + &
                         unk(ifrom,i2,j1,k1,b)+&
                         unk(ifrom,i1,j2,k1,b)+&
                         unk(ifrom,i2,j2,k1,b)+&
                         unk(ifrom,i1,j1,k2,b)+&
                         unk(ifrom,i2,j1,k2,b)+&
                         unk(ifrom,i1,j2,k2,b)+&
                         unk(ifrom,i2,j2,k2,b))
                 enddo
              enddo
           enddo

!!$        ! different method:
!!$        
!!$        call amr_restrict_unk_genorder(unk(:,:,:,:,b), dataout, order, ifrom)
!!$        
!!$        ! put data in more compact form:
!!$        do k = 1+NGUARD, NZB+NGUARD, 2
!!$           kk = (k-NGUARD)/2 + 1 + NGUARD
!!$           do j = 1+NGUARD, NYB+NGUARD, 2
!!$              jj = (j-NGUARD)/2 + 1 + NGUARD
!!$              do i = 1+NGUARD, NXB+NGUARD, 2
!!$                 ii = (i-NGUARD)/2 + 1 + NGUARD
!!$               
!!$                 send_restrict_data(ii,jj,kk,h) = dataout(ifrom,i,j,k)
!!$                 
!!$                 if(dataout(ifrom,i,j,k) .lt. 0.0) then
!!$                       call Driver_abortFlash("negative value in restrict")
!!$                 end if
!!$                 
!!$                 if(dataout(ifrom,i,j,k) .gt. 1.0) then
!!$                    call Driver_abortFlash("Too high a value in restrict")
!!$                 end if
!!$
!!$              enddo
!!$           enddo
!!$        enddo
        
        
        ! Next, if parent is on this processor, copy the restricted values directly.
        ! If parent is off-processor, send the data to the owning processor
        ! (non-blocking).
        if (parent(2,b) == gr_meshMe) then ! local parent
           ! p = block ID of parent if on same proc
           p = parent(1,b)
           do ichild = 1, nchild
              if (child(1,ichild,p) == b) then
                 c = ichild
                 exit
              endif
              if (ichild == nchild) then
                call Driver_abortFlash("[smooth_Restrict] could not find child block");
              endif
           enddo
           ! first method:
           
           if ((c == 1) .or. (c == 3) .or. (c == 5) .or. (c == 7)) then
              i1 = hg_ili
              i2 = hg_ili + hg_restrict_n1 - 1
           else
              i1 = hg_iui - hg_restrict_n1 + 1
              i2 = hg_iui
           endif
           if ((c == 1) .or. (c == 2) .or. (c == 5) .or. (c == 6)) then
              j1 = hg_jli
              j2 = hg_jli + hg_restrict_n2 - 1
           else
              j1 = hg_jui - hg_restrict_n2 + 1
              j2 = hg_jui
           endif
           if ((c == 1) .or. (c == 2) .or. (c == 3) .or. (c == 4)) then
              k1 = hg_kli
              k2 = hg_kli + hg_restrict_n3 - 1
           else
              k1 = hg_kui - hg_restrict_n3 + 1
              k2 = hg_kui
           endif
           unk(ito,i1:i2,j1:j2,k1:k2,p) = send_restrict_data(:,:,:,h)
           
!!$           ! second method - adapted from mpi_amr_1blk_restrict.F90
!!$
!!$
!!$           ic = n1off(c) + NGUARD
!!$           jc = n2off(c) + NGUARD
!!$           kc = n3off(c) + NGUARD
!!$           do i = 1, hg_restrict_n1
!!$              ii = i + ic
!!$              do j = 1, hg_restrict_n2
!!$                 jj = j + jc
!!$                 do k = 1, hg_restrict_n3
!!$                    kk = k + kc
!!$                    
!!$                    unk(ito,ii,jj,kk,p) = send_restrict_data(i,j,k,h)
!!$                    
!!$                 enddo
!!$              enddo
!!$           enddo
       
           
        else                          ! remote parent
           nsent = nsent + 1
           
           call mpi_issend(send_restrict_data(1,1,1,h), &
                hg_restrict_n1*hg_restrict_n2*hg_restrict_n3, &
                FLASH_REAL, parent(2,b), b, &
                gr_meshComm, send_restrict_req(nsent), ierr)

           h = h + 1
        endif
        
        if (h > nbuf_restrict) &
             call Driver_abortFlash("Buffer space exceeded in smooth_Restrict")
        
     endif
  enddo
  
  do b = 1, lnblocks
     if ((lrefine(b) == level-1) .and. (nodetype(b) .gt. LEAF)) then
        do c = 1, nchild
           if (child(2,c,b) /= gr_meshMe) then
              
              ! If child is on another processor, receive the restricted data from
              ! the child (blocking).
              
              call mpi_recv(recv_restrict_data(1,1,1), &
                   hg_restrict_n1*hg_restrict_n2*hg_restrict_n3, &
                   FLASH_REAL, child(2,c,b), child(1,c,b), &
                   gr_meshComm, status, ierr)
              
              ! old method:
              if ((c == 1) .or. (c == 3) .or. (c == 5) .or. (c == 7)) then
                 i1 = hg_ili
                 i2 = hg_ili + hg_restrict_n1 - 1
              else
                 i1 = hg_iui - hg_restrict_n1 + 1
                 i2 = hg_iui
              endif
              if ((c == 1) .or. (c == 2) .or. (c == 5) .or. (c == 6)) then
                 j1 = hg_jli
                 j2 = hg_jli + hg_restrict_n2 - 1
              else
                 j1 = hg_jui - hg_restrict_n2 + 1
                 j2 = hg_jui
              endif
              if ((c == 1) .or. (c == 2) .or. (c == 3) .or. (c == 4)) then
                 k1 = hg_kli
                 k2 = hg_kli + hg_restrict_n3 - 1
              else
                 k1 = hg_kui - hg_restrict_n3 + 1
                 k2 = hg_kui
              endif
              
              unk(ito,i1:i2,j1:j2,k1:k2,b) = recv_restrict_data(:,:,:)
              
!!$              ! new method:
!!$              ic = n1off(c) + NGUARD
!!$              jc = n2off(c) + NGUARD
!!$              kc = n3off(c) + NGUARD
!!$              do i = 1, hg_restrict_n1
!!$                 ii = i + ic
!!$                 do j = 1, hg_restrict_n2
!!$                    jj = j + jc
!!$                    do k = 1, hg_restrict_n3
!!$                       kk = k + kc
!!$                       
!!$                       unk(ito,ii,jj,kk,b) = recv_restrict_data(i,j,k)
!!$                       
!!$                    enddo
!!$                 enddo
!!$              enddo
              
              
           endif
        enddo
     endif
  enddo
  
  call mpi_waitall(nsent, send_restrict_req, send_status, ierr)
  
  !call timer_stop("gr_hgRestrict")
  
  !==============================================================================
  
  return
end subroutine smooth_Restrict
