
!! build logarithmic indices
subroutine phch_makeLogIndices(list, inum, imin, imax)
  implicit none
  
  real, dimension(inum), intent(out) :: list
  integer, intent(in) :: inum
  real, intent(in) :: imin, imax
  integer :: i
  
  do i=1,inum
    list(i) = imin * (imax/imin)**((i-1.0)/(inum-1.0))
  end do
  
end subroutine


