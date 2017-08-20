!-----------------------------------------------------------------------------
!
!  Copyright (C) 1997-2010 Krzysztof M. Gorski, Eric Hivon,
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke
!
!
!  This file is part of HEALPix.
!
!  HEALPix is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  HEALPix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with HEALPix; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!  For more information about HEALPix see http://healpix.jpl.nasa.gov
!
!-----------------------------------------------------------------------------
!
!  Modified APRIL 2011 by Jeremy Ritter for use with FLASH AMR to HealPIX mapping
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!=======================================================================
!     pix2ang_ring
!
!     renders theta and phi coordinates of the nominal pixel center
!     for the pixel number ipix (RING scheme)
!     given the map resolution parameter nside
!
!=======================================================================
  subroutine pix2ang(ipix, theta, phi)
  
    use Driver_interface, ONLY : Driver_abortFlash
    use PhotoChem_data, ONLY: phch_nside, phch_npix

implicit none

#include "constants.h"
    
    INTEGER, INTENT(IN)  :: ipix
    REAL,     INTENT(OUT) :: theta, phi

    INTEGER ::  nl2, nl4, iring, iphi
    INTEGER ::  npix, ncap, ip, nside
    REAL ::  fodd, dnside
    real, parameter :: half = 0.500000000000000
    real, parameter :: one  = 1.000000000000000
    real, parameter :: three = 3.00000000000000
    real, parameter :: threehalf = 1.50000000000000
    real, parameter :: HALFPI = 0.5 * PI
    !-----------------------------------------------------------------------
    
    nside = phch_nside
    
    if (nside < 1) call Driver_abortFlash("nside out of range")

    npix = phch_npix
    if (ipix <0 .or. ipix>npix-1) call Driver_abortFlash ("ipix out of range")

    nl2  = 2*nside
    ncap = nl2*(nside-1) ! points in each polar cap, =0 for nside =1
    dnside = real(nside)

    if (ipix < ncap) then ! North Polar cap -------------

       iring = nint( sqrt( (ipix+1) * half )) ! counted from North pole
       iphi  = ipix - 2*iring*(iring - 1)

       theta = ACOS( one - (iring/dnside)**2 / three )
       phi   = (real(iphi) + half) * HALFPI/iring

    elseif (ipix < npix-ncap) then ! Equatorial region ------

       ip    = ipix - ncap
       nl4   = 4*nside
       iring = INT( ip / nl4 ) + nside ! counted from North pole
       iphi  = iand(ip, nl4-1)

       fodd  = half * ( iand(iring+nside+1,1) )  ! 0 if iring+nside is odd, 1/2 otherwise
       theta = ACOS( (nl2 - iring) / (threehalf*dnside) )
       phi   = (real(iphi) + fodd) * HALFPI / dnside

    else ! South Polar cap -----------------------------------

       ip    = npix - ipix
       iring = nint( sqrt( ip * half ))     ! counted from South pole
       iphi  = 2*iring*(iring + 1) - ip

       theta = ACOS( (iring/dnside)**2 / three  - one)
       phi   = (real(iphi) + half) * HALFPI/iring

    endif

    return
  end subroutine pix2ang

!=======================================================================
!     ang2pix_ring
!
!     renders the pixel number ipix (RING scheme) for a pixel which contains
!     a point on a sphere at coordinates theta and phi, given the map
!     resolution parameter nside
!=======================================================================
  subroutine ang2pix(theta, phi, ipix)

    use Driver_interface, ONLY : Driver_abortFlash
    use PhotoChem_data, ONLY : phch_nside

implicit none
    
    REAL,     INTENT(IN)  :: theta, phi
    INTEGER, INTENT(OUT) :: ipix

#include "constants.h"

    INTEGER   ::  nl4, jp, jm, ir, ip, nside
    REAL     ::  z, za, tt, tp, tmp, temp1, temp2
    INTEGER ::  kshift

    REAL, PARAMETER :: pi = PI;
    REAL, PARAMETER :: twopi = 2.0 * PI;
    REAL, PARAMETER :: halfpi = 0.5 * PI;
    real, parameter :: twothird = 2.0/3.0

    !-----------------------------------------------------------------------
    nside = phch_nside
    
    if (nside.lt.1) call Driver_abortFlash("nside out of range")

    if (theta<0.0 .or. theta>pi)  then
       print*,"ANG2PIX_RING: theta : ",theta," is out of range [0, Pi]"
       call Driver_abortFlash("ANG2PIX: theta out of range")
    endif

    z = COS(theta)
    za = ABS(z)
    tt = MODULO( phi, twopi) / halfpi  ! in [0,4)

    nl4 = 4*nside
    if ( za <= twothird ) then ! Equatorial region ------------------
       temp1 = nside*(0.50000 + tt)
       temp2 = nside* 0.75000 * z
       jp = int(temp1-temp2) ! index of  ascending edge line
       jm = int(temp1+temp2) ! index of descending edge line

       ir = nside + 1 + jp - jm 
       
       kshift = 0
       if (MODULO(ir,2).eq.0.0) kshift = 1

       ip = INT( ( jp+jm - nside + kshift + 1 ) / 2 ) + 1
       if (ip > nl4) ip = ip - nl4

       ipix = 2*nside*(nside-1) + nl4*int(ir-1) + ip

    else ! North & South polar caps -----------------------------

       tp = tt - INT(tt)      !MODULO(tt,1.0)
       tmp = nside * SQRT( 3.0*(1.0 - za) )

       jp = INT(  tp          * tmp) ! increasing edge line index
       jm = INT((1.0 - tp) * tmp) ! decreasing edge line index

       ir = jp + jm + 1        ! ring number counted from the closest pole
       ip = INT( tt * ir) + 1     ! in {0,4*ir-1}
       if (ip > 4*ir) ip = ip - 4*ir

       if (z>0.0) then
          ipix = 2*ir*(ir-1) + ip
       else
          ipix = 3*nside*int(nl4) - 2*ir*(ir+1) + ip
       endif

    endif

    return
  end subroutine ang2pix

!=======================================================================
!     vec2pix_ring
!
!     renders the pixel number ipix (RING scheme) for a pixel which contains
!     a point on a sphere at coordinate vector (=x,y,z), given the map
!     resolution parameter nside
!=======================================================================
  subroutine vec2pix  (vector, ipix)

    use Driver_interface, ONLY: Driver_abortFlash
    use PhotoChem_data, ONLY: phch_nside

    REAL,     INTENT(IN), dimension(1:) :: vector
    INTEGER, INTENT(OUT)               :: ipix

    INTEGER   :: nl4, jp, jm, ir, ip
    REAL     :: z, za, tt, tp, tmp, dnorm, phi,temp1,temp2
    INTEGER :: kshift
    
    real, parameter :: halfpi = 0.5 * PI
    REAL, PARAMETER :: twopi = 2.0 * PI;
    real, parameter :: twothird = 2.0/3.0

    !-----------------------------------------------------------------------
    dnorm = SQRT(vector(1)**2+vector(2)**2+vector(3)**2)
    z = vector(3) / dnorm
    phi = 0.0
    if (vector(1) /= 0.0 .or. vector(2) /= 0.0) &
         &     phi = ATAN2(vector(2),vector(1)) ! phi in ]-pi,pi]

    za = ABS(z)
    if (phi < 0.0)     phi = phi + twopi ! phi in [0,2pi[
    tt = phi / halfpi   ! in [0,4)

    nl4 = 4*phch_nside
    if ( za <= twothird ) then ! Equatorial region ------------------
       temp1 = phch_nside*(0.50000 + tt)
       temp2 = phch_nside* 0.75000 * z
       jp = int(temp1-temp2) ! index of  ascending edge line
       jm = int(temp1+temp2) ! index of descending edge line
!        jp = INT(phch_nside*(0.5 + tt - z*0.75)) index of  ascending edge line
!        jm = INT(phch_nside*(0.5 + tt + z*0.75)) index of descending edge line

       ir = phch_nside + jp - jm ! in {0,2n} (ring number counted from z=2/3)
       kshift = iand(ir, 1) ! kshift=1 if ir is odd, 0 otherwise

       ip = INT( ( jp+jm - phch_nside + kshift + 1 ) / 2) ! in {0,4n-1}
       if (ip >= nl4) ip = ip - nl4

       ipix = 2*phch_nside*(phch_nside-1) + nl4*int(ir) + ip

    else ! North & South polar caps -----------------------------

       tp = tt - INT(tt)      !MODULO(tt,1.0)
       tmp = phch_nside * SQRT( 3.0*(1.0 - za) )

       jp = INT(  tp          * tmp) ! increasing edge line index
       jm = INT((1.0 - tp) * tmp) ! decreasing edge line index

       ir = jp + jm + 1        ! ring number counted from the closest pole
       ip = INT( tt * ir)     ! in {0,4*ir-1}
       if (ip >= 4*ir) ip = ip - 4*ir

       if (z>0.) then
          ipix = 2*ir*(ir-1) + ip
       else
          ipix = 3*phch_nside*int(nl4) - 2*ir*(ir+1) + ip
       endif

    endif

    return
  end subroutine vec2pix

!=======================================================================
!     pix2vec_ring
!
!     renders vector (x,y,z) coordinates of the nominal pixel center
!     for the pixel number ipix (RING scheme)
!     given the map resolution parameter phch_nside
!     also returns the (x,y,z) position of the 4 pixel vertices (=corners)
!     in the order N,W,S,E
!=======================================================================

  subroutine pix2vec  (ipix, vector, vertex)

    use Driver_interface, ONLY: Driver_abortFlash
    use PhotoChem_data, ONLY: phch_nside, phch_npix

    INTEGER, INTENT(IN)                             :: ipix
    REAL,     INTENT(OUT),dimension(1:)              :: vector
    REAL,     INTENT(OUT),dimension(1:,1:), optional :: vertex

    INTEGER :: nl2, nl4, iring, iphi
    INTEGER ::  npix, ncap, ip
    REAL ::  fact1, fact2, fodd, z, sth, phi
    real, parameter :: half = 0.500000000000000
    real, parameter :: halfpi = 0.5 * PI

    real :: phi_nv, phi_wv, phi_sv, phi_ev, sin_phi, cos_phi
    real :: z_nv, z_sv, sth_nv, sth_sv
    real :: hdelta_phi
    integer :: iphi_mod, iphi_rat
    logical :: do_vertex
    integer :: diff_phi
    !-----------------------------------------------------------------------

    npix = phch_npix
    if (ipix <0 .or. ipix>npix-1) call Driver_abortFlash("ipix out of range")

    nl2   = 2*phch_nside
    ncap  = nl2*(phch_nside-1) ! points in each polar cap, =0 for phch_nside =1

    do_vertex = .false.
    if (present(vertex)) then
       if (size(vertex,dim=1) >= 3 .and. size(vertex,dim=2) >= 4) then
          do_vertex = .true.
       else
          call Driver_abortFlash(" pix2vec_ring : vertex array has wrong size ")
       endif
    endif

    if (ipix < ncap) then ! North Polar cap -------------

       iring = nint( sqrt( (ipix+1) * half )) ! counted from North pole
       iphi  = ipix - 2*iring*(iring - 1)

       fact2 = (3.00000*phch_nside)*phch_nside
       z =  1.0 - (iring / fact2) * iring
       phi   = (real(iphi) + half) * HALFPI/iring

       if (do_vertex) then
          hdelta_phi = PI/(4.0*iring)   ! half pixel width
          z_nv = 1.0 - (iring-1) / fact2 * (iring-1)
          z_sv = 1.0 - (iring+1) / fact2 * (iring+1)
          iphi_mod = MODULO(iphi, iring) ! in {0,1,... iring-1}
          iphi_rat = (iphi) / iring      ! in {0,1,2,3}
          phi_nv = 0.0
          if (iring > 1) phi_nv = HALFPI * (iphi_rat +  iphi_mod   /real(iring-1))
          phi_sv                = HALFPI * (iphi_rat + (iphi_mod+1)/real(iring+1))
          diff_phi = 3 ! both phi_nv and phi_sv different from phi
       endif


    elseif (ipix < npix - ncap) then ! Equatorial region ------

       nl4 = 4*phch_nside
       ip    = ipix - ncap
       iring = INT( ip / nl4 ) + phch_nside ! counted from North pole
       iphi  = iand(ip, nl4-1)

       fact1 =  1.50000*phch_nside
       fodd  = half * ( iand(iring+phch_nside+1,1) )  ! 0 if iring+phch_nside is odd, 1/2 otherwise
       z = (nl2 - iring) / fact1
       phi   = (real(iphi) + fodd) * HALFPI / phch_nside

       if (do_vertex) then
          fact2 = (3.00000*phch_nside)*phch_nside
          hdelta_phi = PI/(4.0*phch_nside)   ! half pixel width
          phi_nv = phi
          phi_sv = phi
          z_nv = (nl2 - iring +1) / fact1
          z_sv = (nl2 - iring -1) / fact1
          diff_phi = 0 ! phi_nv = phi_sv = phi
          if (iring == phch_nside) then ! northern transition
             z_nv = 1.0 - (phch_nside-1) / fact2 * (phch_nside-1)
             iphi_mod = iand(iphi, phch_nside-1) ! in {0,1,... phch_nside-1}
             iphi_rat = (iphi) / phch_nside      ! in {0,1,2,3}
             if (phch_nside > 1) then
                phi_nv = HALFPI * (iphi_rat +  iphi_mod /real(phch_nside-1))
                diff_phi = 1
             endif
          elseif (iring == 3*phch_nside) then ! southern transition
             z_sv = -1.0 + (phch_nside-1) / fact2 * (phch_nside-1)
             iphi_mod = iand(iphi, phch_nside-1) ! in {0,1,... iring-1}
             iphi_rat = (iphi) / phch_nside      ! in {0,1,2,3}
             if (phch_nside > 1) then
                phi_sv = HALFPI * (iphi_rat +  iphi_mod   /real(phch_nside-1))
                diff_phi = 2
             endif
          endif
       endif

    else ! South Polar cap -----------------------------------

       ip    = npix - ipix
       iring = nint( sqrt( ip * half ))     ! counted from South pole
       iphi  = 2*iring*(iring + 1) - ip

       fact2 = (3.00000*phch_nside)*phch_nside
       z = -1.0 + (iring / fact2) * iring
       phi   = (real(iphi) + half) * HALFPI/iring

       if (do_vertex) then
          hdelta_phi = PI/(4.0*iring)   ! half pixel width
          z_nv = -1.0 + (iring+1)/ fact2 * (iring+1)
          z_sv = -1.0 + (iring-1)/ fact2 * (iring-1)
          iphi_mod = MODULO(iphi, iring) ! in {0,1,... iring-1}
          iphi_rat = (iphi) / iring      ! in {0,1,2,3}
          phi_nv                = HALFPI * (iphi_rat + (iphi_mod+1)/real(iring+1))
          phi_sv = 0.0
          if (iring > 1) phi_sv = HALFPI * (iphi_rat +  iphi_mod   /real(iring-1))
          diff_phi = 3
       endif

    endif

    ! pixel center
    sth = SQRT((1.0-z)*(1.0+z))
    cos_phi = cos(phi)
    sin_phi = sin(phi)
    vector(1) = sth * cos_phi
    vector(2) = sth * sin_phi
    vector(3) = z

    if (do_vertex) then
       ! west vertex
       phi_wv      = phi - hdelta_phi
       vertex(1,2) = sth * COS(phi_wv)
       vertex(2,2) = sth * SIN(phi_wv)
       vertex(3,2) = z

       ! east vertex
       phi_ev      = phi + hdelta_phi
       vertex(1,4) = sth * COS(phi_ev)
       vertex(2,4) = sth * SIN(phi_ev)
       vertex(3,4) = z

       ! north and south vertices
       sth_nv = SQRT((1.0-z_nv)*(1.0+z_nv))
       sth_sv = SQRT((1.0-z_sv)*(1.0+z_sv))
       if (diff_phi == 0) then
          vertex(1,1) = sth_nv * cos_phi
          vertex(2,1) = sth_nv * sin_phi
          vertex(1,3) = sth_sv * cos_phi
          vertex(2,3) = sth_sv * sin_phi
       else
          vertex(1,1) = sth_nv * COS(phi_nv)
          vertex(2,1) = sth_nv * SIN(phi_nv)
          vertex(1,3) = sth_sv * COS(phi_sv)
          vertex(2,3) = sth_sv * SIN(phi_sv)
       endif
       vertex(3,1) = z_nv
       vertex(3,3) = z_sv
    endif

    return
  end subroutine pix2vec


  !=======================================================================
  subroutine ang2vec(theta, phi, vector)
    !=======================================================================
    !     renders the vector (x,y,z) corresponding to angles
    !     theta (co-latitude measured from North pole, in [0,Pi] radians)
    !     and phi (longitude measured eastward, in radians)
    !     North pole is (x,y,z)=(0,0,1)
    !     added by EH, Feb 2000
    !=======================================================================
    REAL, INTENT(IN) :: theta, phi
    REAL, INTENT(OUT), dimension(1:) :: vector

    REAL :: sintheta
    !=======================================================================

    if (theta<0.0 .or. theta>PI)  then
       print*,"ANG2VEC: theta : ",theta," is out of range [0, Pi]"
    endif
    sintheta = SIN(theta)

    vector(1) = sintheta * COS(phi)
    vector(2) = sintheta * SIN(phi)
    vector(3) = COS(theta)

    return
  end subroutine ang2vec

!=======================================================================
  subroutine vec2ang(x,y,z, r, theta, phi)
    !=======================================================================
    !     renders the angles theta, phi corresponding to vector (x,y,z)
    !     theta (co-latitude measured from North pole, in [0,Pi] radians)
    !     and phi (longitude measured eastward, in [0,2Pi[ radians)
    !     North pole is (x,y,z)=(0,0,1)
    !     added by EH, Feb 2000
    !=======================================================================

#include "constants.h"

    REAL, INTENT(IN) :: x, y, z
    REAL, INTENT(OUT) :: r, theta, phi

    real, parameter :: twopi = 2.0 * PI

    !=======================================================================

    r = SQRT(x**2+y**2+z**2)
    
    if (r.eq.0.0) then
      theta = ACOS( 0.0 )
    else
      theta = ACOS(z / r)
    end if

    phi = 0.0
    if (x /= 0.0 .or. y /= 0.0) &
         &     phi = ATAN2(y,x) ! phi in ]-pi,pi]
    if (phi < 0.0)     phi = phi + twopi ! phi in [0,2pi[

    return
  end subroutine vec2ang


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine xyz2pix(x,y,z,pix)
  use phch_interface, ONLY : vec2ang, ang2pix

  implicit none
  
  real, intent(IN) :: x, y, z
  integer, intent(OUT) :: pix
  
  real :: r, theta, phi

  call vec2ang(x, y, z, r, theta, phi)
  
  call ang2pix(theta, phi, pix)

end subroutine xyz2pix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine xyz2bin(x, y, z, pix, bin)
  use PhotoChem_data, ONLY: phch_bin_min, phch_nbin, phch_bin_log_max_over_min
  use phch_interface, ONLY : vec2ang, ang2pix

implicit none

  real, intent(IN) :: x, y, z
  integer, intent(OUT) :: pix, bin
  
  real :: r, theta, phi
  
  call vec2ang(x, y, z, r, theta, phi)
  
  call ang2pix(theta, phi, pix)

  if (r.lt.phch_bin_min) then
    bin = 0
    return
  end if  

  !!bin = int(phch_nbin * (log10(r / phch_bin_max) + 1))
  
  bin = int(floor(real(phch_nbin) * log10(r/phch_bin_min)/phch_bin_log_max_over_min)) + 1

  if (bin.gt.phch_nbin) bin = -1

  return
end subroutine xyz2bin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! given a point in space, x y z, and the half-dimensions
!! of a cube, dx dy dz, determine if the cube lies entirely
!! within a pixel. if it does not, split the cell into 8ths
!! and try again. if it is in a pixel, add up volumetric weights

recursive subroutine cube2pix(x,y,z,dx,dy,dz,level)
  use PhotoChem_data, ONLY: phch_max_octree_level, phch_pixel_volume, phch_pixel_number, phch_current_dens, &
                             phch_octree_level, phch_pixel_lambda, phch_current_alpha, &
                             phch_binzero_volume, phch_binzero_number, phch_binzero_lambda, &
                             phch_binzero_hplus, phch_pixel_hplus, phch_total_h_fraction, phch_hplus_fraction, phch_nbin

  use phch_interface, ONLY : xyz2bin

implicit none

  real, intent(IN) :: x, y, z, dx, dy, dz
  integer, intent(IN) :: level
  
  integer :: pixA, pixB, binA, binB
  real :: hdx, hdy, hdz, volume, mass
  integer :: nextlevel, pix, bin, idx, min_bin, max_level
  logical :: splitcell
  
  real, dimension(8,3) :: pos
  
  real :: ntotal
  
  if (level.gt.phch_octree_level) phch_octree_level = level

  splitcell = .false.
  
  nextlevel = level+1
  
  call xyz2bin(x,y,z,pix,bin)
  
  ! we can exit early if the cell center is outside the last bin
  ! the max_level would be 0 and it won't split or be added to anything
  if (bin.lt.0) return
  
  max_level = phch_max_octree_level - floor(real(phch_max_octree_level) * real(bin) / real(1+phch_nbin))

  if (nextlevel.le.max_level) then
    ! label opposite pair of corners so we can check 1-4 with 5-8
    ! can exit early if we find a non-matching pair
    pos(1,:) = (/ x-dx, y-dy, z-dz /)
    pos(5,:) = (/ x+dx, y+dy, z+dz /)

    pos(2,:) = (/ x+dx, y-dy, z-dz /)
    pos(6,:) = (/ x-dx, y+dy, z+dz /)

    pos(3,:) = (/ x+dx, y+dy, z-dz /)
    pos(7,:) = (/ x-dx, y-dy, z+dz /)

    pos(4,:) = (/ x-dx, y+dy, z-dz /)
    pos(8,:) = (/ x+dx, y-dy, z+dz /)

    do idx=1,4
      call xyz2bin(pos(idx,1), pos(idx,2), pos(idx,3), pixA, binA)
      call xyz2bin(pos(idx+4,1), pos(idx+4,2), pos(idx+4,3), pixB, binB)

      if ((pixA.ne.pixB).or.(binA.ne.binB)) then
        splitcell = .true.
        exit
      endif
    enddo
  endif

  if (splitcell) then
  
    hdx = 0.5*dx
    hdy = 0.5*dy
    hdz = 0.5*dz
  
    call cube2pix(x-hdx,y-hdy,z-hdz,hdx,hdy,hdz,nextlevel)
    call cube2pix(x+hdx,y-hdy,z-hdz,hdx,hdy,hdz,nextlevel)
    call cube2pix(x+hdx,y+hdy,z-hdz,hdx,hdy,hdz,nextlevel)
    call cube2pix(x-hdx,y+hdy,z-hdz,hdx,hdy,hdz,nextlevel)
    call cube2pix(x-hdx,y-hdy,z+hdz,hdx,hdy,hdz,nextlevel)
    call cube2pix(x+hdx,y-hdy,z+hdz,hdx,hdy,hdz,nextlevel)
    call cube2pix(x+hdx,y+hdy,z+hdz,hdx,hdy,hdz,nextlevel)
    call cube2pix(x-hdx,y+hdy,z+hdz,hdx,hdy,hdz,nextlevel)
  
  else             !! add current cube to its center's pixel
    
    volume = 8.0*dx*dy*dz
    mass = phch_current_dens * volume
    ntotal = phch_total_h_fraction * phch_current_dens

    if (bin.eq.0) then
      phch_binzero_volume = phch_binzero_volume + volume
      phch_binzero_number = phch_binzero_number + phch_total_h_fraction * mass
      phch_binzero_lambda = phch_binzero_lambda + phch_current_alpha * ntotal**2.0 * volume
      phch_binzero_hplus = phch_binzero_hplus + mass * phch_hplus_fraction
    else if (bin.gt.0) then
      phch_pixel_volume(pix,bin) = phch_pixel_volume(pix,bin) + volume
      phch_pixel_number(pix,bin) = phch_pixel_number(pix,bin) + phch_total_h_fraction * mass
      phch_pixel_lambda(pix,bin) = phch_pixel_lambda(pix,bin) + phch_current_alpha * ntotal**2.0 * volume
      phch_pixel_hplus(pix,bin) = phch_pixel_hplus(pix,bin) + mass*phch_hplus_fraction
    endif
    
  endif
  
end subroutine cube2pix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cell2pix(x,y,z,dx,dy,dz, dens, temp)

  use PhotoChem_data, ONLY: phch_current_dens, phch_octree_level, phch_current_alpha
  use phch_interface, ONLY : cube2pix, phch_calculateRate

implicit none

  real, intent(IN) :: x,y,z,dx,dy,dz, dens, temp
  
  phch_current_dens = dens
  call phch_calculateRate(temp, phch_current_alpha)
  phch_octree_level = 0
  
  call cube2pix(x,y,z,dx,dy,dz,0)

end subroutine cell2pix

