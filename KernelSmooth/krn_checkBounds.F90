
! roughly check if the particle smoothing volume overlaps with a block bounding box
! rough smoothing volume defined by center px,py,pz and sides 2*len

recursive subroutine krn_checkBounds(px, py, pz, len, bndBox, overlap, again)

  use Kernel_data

  implicit none

#include "constants.h"
  
  real, intent(in) :: px, py, pz, len
  real, intent(in), dimension(LOW:HIGH, MDIM) :: bndBox
  logical, intent(inout) :: overlap
  logical, intent(in) :: again

  real, dimension(LOW:HIGH, MDIM) :: partBox

  ! exit early if we already detected overlap
  if (overlap) return

  partBox(LOW, IAXIS) = px - len
  partBox(HIGH, IAXIS) = px + len
  partBox(LOW, JAXIS) = py - len
  partBox(HIGH, JAXIS) = py + len
  partBox(LOW, KAXIS) = pz - len
  partBox(HIGH, KAXIS) = pz + len
  
  overlap = overlap.or.&
            (( (bndBox(LOW,IAXIS).le.partBox(HIGH,IAXIS)) .and. (bndBox(HIGH,IAXIS).ge.partbox(LOW,IAXIS)) ).and. &
             ( (bndBox(LOW,JAXIS).le.partBox(HIGH,JAXIS)) .and. (bndBox(HIGH,JAXIS).ge.partbox(LOW,JAXIS)) ).and. &
             ( (bndBox(LOW,KAXIS).le.partBox(HIGH,KAXIS)) .and. (bndBox(HIGH,KAXIS).ge.partbox(LOW,KAXIS)) ))

  if (.not.overlap.and.krn_periodic.and.again) then
    if ((px - len).lt.krn_xmin) call krn_checkBounds(px+krn_dx, py, pz, len, bndBox, overlap, .false.)
    if ((px + len).gt.krn_xmax) call krn_checkBounds(px-krn_dx, py, pz, len, bndBox, overlap, .false.)

    if ((py - len).lt.krn_ymin) call krn_checkBounds(px, py+krn_dy, pz, len, bndBox, overlap, .false.)
    if ((py + len).gt.krn_ymax) call krn_checkBounds(px, py-krn_dy, pz, len, bndBox, overlap, .false.)
    
    if ((pz - len).lt.krn_zmin) call krn_checkBounds(px, py, pz+krn_dz, len, bndBox, overlap, .false.)
    if ((pz + len).gt.krn_zmax) call krn_checkBounds(px, py, pz-krn_dz, len, bndBox, overlap, .false.)
  end if

end subroutine
  
