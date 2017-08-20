
! kernel weighting function
! q = r / 2h

subroutine krn_kernel(q,w)
  implicit none
  
  real, intent(in) :: q
  real, intent(inout) :: w


  if (q.gt.1.0) then
    return
  else if (q.gt.0.5) then
    w = w + 2.0 * (1.0-q)**3
  else if (q.gt.0.0) then
    w = w + 1.0 - 6.0*q*q + 6.0*q*q*q
  end if
  
end subroutine

