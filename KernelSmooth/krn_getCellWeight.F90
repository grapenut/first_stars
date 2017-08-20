
! calculate weight for mapping particle at px,py,pz to cell at cx,cy,cz

recursive subroutine krn_getCellWeight(px, py, pz, cx, cy, cz, len, wgt, again)

  use Kernel_data
  use Kernel_interface, ONLY : krn_kernel

  implicit none
  
  real, intent(in) :: px, py, pz, cx, cy, cz, len
  real, intent(inout) :: wgt
  logical, intent(in) :: again  !check to prevent multiple levels of recursion
  
  real :: dx, dy, dz, dr

  if (wgt.gt.0.0) return
  
  dx = cx - px
  dy = cy - py
  dz = cz - pz
  dr = sqrt(dx*dx + dy*dy + dz*dz)

  if (dr.lt.len) then
    call krn_kernel(dr/len, wgt)
  end if
  
  if (.not.(wgt.gt.0.0).and.krn_periodic.and.again) then
    if ((px - len).lt.krn_xmin) call krn_getCellWeight(px+krn_dx, py, pz, cx, cy, cz, len, wgt, .false.)
    if ((px + len).gt.krn_xmax) call krn_getCellWeight(px-krn_dx, py, pz, cx, cy, cz, len, wgt, .false.)

    if ((py - len).lt.krn_ymin) call krn_getCellWeight(px, py+krn_dy, pz, cx, cy, cz, len, wgt, .false.)
    if ((py + len).gt.krn_ymax) call krn_getCellWeight(px, py-krn_dy, pz, cx, cy, cz, len, wgt, .false.)
    
    if ((pz - len).lt.krn_zmin) call krn_getCellWeight(px, py, pz+krn_dz, cx, cy, cz, len, wgt, .false.)
    if ((pz + len).gt.krn_zmax) call krn_getCellWeight(px, py, pz-krn_dz, cx, cy, cz, len, wgt, .false.)
  end if

end subroutine

