
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine phch_calculateRate(temp, alpha)
  
  use PhotoChem_data, ONLY : phch_current_invscale3, phch_current_scale2

implicit none

  real, intent(in) :: temp
  real, intent(out) :: alpha
  
  real :: T
  
  T = temp * phch_current_scale2
  
  !alpha = 3e-13 * phch_current_invscale3
  
  !alpha = 4.881e-6 * T**(-1.5) * (1.0 + 114.8 * T**(-0.407) )**(-2.242) * phch_current_invscale3
  
  ! this is in m^3 / s
  !alpha = 2.0e-16 * temp**(-0.75)

  alpha = 2.90e-10 * temp**(-0.77) * phch_current_invscale3

  return

end subroutine phch_calculateRate

