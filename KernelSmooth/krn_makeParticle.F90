
! create a single smoothed particle with the given quantities

subroutine krn_makeParticle(x, y, z, mass, px, py, pz, temp, len)

  use Kernel_data

  implicit none
  
  real, intent(in) :: x, y, z, mass, px, py, pz, temp, len
  integer :: num
  
  krn_num_particles = krn_num_particles + 1

  if (krn_num_particles.gt.krn_max_particles) then
    print *,'KERNEL: krn_num_particles too big', krn_num_particles, krn_max_particles
  end if

  num = krn_num_particles

  krn_particles(krn_prop_x,num) = x
  krn_particles(krn_prop_y,num) = y
  krn_particles(krn_prop_z,num) = z
  krn_particles(krn_prop_mass,num) = mass
  krn_particles(krn_prop_px,num) = px
  krn_particles(krn_prop_py,num) = py
  krn_particles(krn_prop_pz,num) = pz
  krn_particles(krn_prop_temp,num) = temp
  krn_particles(krn_prop_len,num) = len
  krn_particles(krn_prop_wgt,num) = 0.0
  
end subroutine
