!!****if* source/Grid/GridMain/paramesh/gr_unmarkRefineByLogRadius
!!
!! NAME
!!  gr_unmarkRefineByLogRadius
!!
!!  
!! SYNOPSIS 
!!  call gr_unmarkRefineByLogRadius(real(in) :: xc, 
!!                                  real(in) :: yc,
!!                                  real(in) :: zc)
!!  
!! DESCRIPTION
!!  Cancel refinement flags for all blocks that are too far away from a
!!  center given by (xc,yc,zc).
!!  The determination whether a block is 'too far away' depends on the
!!  current block size as well as the distance, and is made using the
!!  runtime parameter gr_lrefineMaxRedRadiusFact.
!!  
!! ARGUMENTS 
!!  xc -   Center of the interval/circle/sphere : IAXIS
!!  yc -                                          JAXIS
!!  zc -                                          KAXIS
!!               (Coordinates for nonexistent dimensions are ignored.)
!!  
!! SIDE EFFECTS
!!
!!  Elements in the PARAMESH logical array refine(:) may be modified.
!!
!! NOTES
!! 
!!  This routine has not been tested well. It probably should be viewed only as a
!!  guideline for a user's implementation.
!!  
!!  If the geometry is SPHERICAL or POLAR, the distance is measured in the radial
!!  direction (X-direction) alone, and is taken as the radial distance from a
!!  sphere of radius given by xc.  In particular, the distance is the distance
!!  from the coordinate center if xc = 0.0.
!!
!!  3D cylindrical geometry is not supported.
!!***

subroutine gr_unmarkRefineByLogRadius(xc, yc, zc)

!-------------------------------------------------------------------------------
  use tree, ONLY : refine, derefine, lrefine, bsize, coord, lnblocks, nodetype, lrefine_min
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_geometry, gr_lrefineMaxRedRadiusSq
  use Cosmology_data, ONLY : csm_scaleFactor
#include "constants.h"
#include "Flash.h"
  implicit none

! Arguments

  real, intent(IN)      :: xc, yc, zc

! Local data

  real, dimension(MDIM) :: blockCenter
  real                  :: RadiusSq
  real                  :: pc, rvir2, rext2, scale2
  real, save		:: radius_virial, radius_exterior
  integer               :: b
  integer, save		:: lref, lref_outside_virial, lref_exterior
  logical, save		:: first_call = .true.

  if (first_call) then
    first_call = .false.
    pc = 3.08e18

    radius_virial = 1.0e4 * pc
    radius_virial = radius_virial * radius_virial
    lref_outside_virial = lrefine_min + 2
    
    radius_exterior = 2.5e4 * pc
    radius_exterior = radius_exterior * radius_exterior
    lref_exterior = lrefine_min
  end if
  scale2 = 1.0 / (csm_scalefactor * csm_scalefactor)
  rvir2 = radius_virial * scale2
  rext2 = radius_exterior * scale2

  do b = 1, lnblocks
     if (nodetype(b) == LEAF) then
        blockCenter(:) = coord(:,b)
        radiusSq = (BlockCenter(1) - xc)**2 + (BlockCenter(2) - yc)**2 + (BlockCenter(3) - zc)**2
        lref = lrefine(b)
        
        if (radiusSq.gt.rext2) then
          if (lref.ge.lref_exterior) then
            refine(b) = .false.
          end if
        else if (radiusSq.gt.rvir2) then
          if (lref.ge.lref_outside_virial) then
            refine(b) = .false.
          !else if (lref.gt.lref_virial) then
          !  refine(b) = .false.
          !  derefine(b) = .true.
          end if
        end if
     end if
  end do

  return
end subroutine gr_unmarkRefineByLogRadius
