
!! support function for finding the number of unique indices
subroutine phch_findNumIndices(list, rows, inum, imin, imax)
  
  implicit none
  
  integer, intent(in) :: rows
  real, dimension(rows), intent(in) :: list
  integer, intent(out) :: inum
  real, intent(out) :: imin
  real, intent(out) :: imax

  logical, dimension(rows) :: mask
  logical :: keep_going
  
  real :: last_value, new_value

  imin = minval(list)
  imax = maxval(list)
  last_value = -1.0
  inum = 0
  keep_going = .true.
  mask = .true.
  do while (keep_going)
    new_value = minval(list, mask)
    if (new_value .gt. last_value) then
      inum = inum + 1
      where(list.le.new_value)
        mask = .false.
      elsewhere
        mask = .true.
      end where
    else
      keep_going = .false.
    end if
    last_value = new_value
  end do
  inum = inum - 1

end subroutine


