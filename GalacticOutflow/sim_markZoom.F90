!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! mark grid for refinement where needed to resolve new supernovae before insertion

! restrict resolution by stepping down max grid refinement by radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sim_markZoom()

  !-------------------------------------------------------------------------------

  use Simulation_data, ONLY : sim_num_stars, sim_zoom_status, sim_galaxy_refine_level, &
                              sim_mype, sim_ismaster, sim_nova_r2, sim_galaxy_r2, &
                              sim_zoom_time, sim_zoom_level, sim_numPEs, sim_starburst_refine_level, &
                              sim_nova_started, sim_nova_time, sim_star_coords, sim_starburst_r2
  use tree, ONLY: refine, derefine, lrefine, bsize, coord, nodetype, lnblocks
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getSimTime, Driver_getDt
  use Grid_data, ONLY : gr_maxRefine
  
  use PhotoChem_data, ONLY : phch_source_x, phch_source_y, phch_source_z
  
  use Cosmology_interface, ONLY : Cosmology_getRedshift
  
  implicit none

#include "constants.h"
#include "Flash.h"
#include "zoom.h"
#include "Flash_mpi.h"

  ! Local data

  real, dimension(MDIM) :: blockCenter, blockSize, bleft, bright, coords, dist2source, dist2star
  integer               :: b, n, ierr
  logical               :: ready(sim_num_stars)
  real			:: time, r2
  
  integer :: zoom_level

  logical :: refine_sn
  
  integer :: global_num_blocks
  
  real :: redshift, invscale, invscale2, galaxy_r2, starburst_r2
  
  call Driver_getSimTime(time)
  call Cosmology_getRedshift(redshift)
  invscale = 1.0 + redshift
  invscale2 = invscale * invscale
  galaxy_r2 = sim_galaxy_r2 * invscale2
  starburst_r2 = sim_starburst_r2 * invscale2 ! starburst_radius
  
  ! are we nearing our global block limit?
  ! we can allow a step down in maximum refinement level to compensate
  !allow_step = global_num_blocks .gt. sim_block_usage_threshold*MAXBLOCKS*sim_numPEs

  ! source refinement location/level
  coords(1) = phch_source_x
  coords(2) = phch_source_y
  coords(3) = phch_source_z
  
  ! reset local_ready from the previous step in case blocks moved processors
  where (sim_zoom_status.eq.ZOOM_LOCAL_READY)
    sim_zoom_status = ZOOM_REFINING
  end where

  ! SN refinement only happens if we are zooming
  refine_sn = any(sim_zoom_status .eq. ZOOM_REFINING)
  ready = .true.

  ! setup some variables used by all blocks
  !lmax = floor(1.0 + log(sim_box_size / (invscale * maxval(sim_nova_radius)))/log(2.0)) + 2
  !delta_lref = lmax - sim_derefine_min

  !rmin = sim_box_size / 2.0**(lmax - 1 - sim_derefine_radius_base)
  !rmax = sim_box_size / 2.0**(sim_derefine_min - 1 - sim_derefine_radius_base)
  !log_rmax_over_rmin = log(rmax / rmin)

  !new_step = sim_derefine_step

  do b = 1, lnblocks
    if (nodetype(b) == LEAF) then
      blockCenter = coord(:,b)
      blockSize  = 0.5 * bsize(:,b)
      
      ! first find the distance from the closest block corner to the photo source
      bleft = blockCenter - blockSize - coords
      bright = blockCenter + blockSize - coords
      
      where (bleft*bright .gt. 0)
        dist2source = min(bleft*bleft, bright*bright)
      elsewhere
        dist2source = 0.0
      end where

      !blockWidth = bsize(:,b)
      !blockWidth = blockWidth * blockWidth  !!! squared
      
      ! make sure the galaxy containing the source is kept at galaxy_refine_level
      if (all(dist2source.lt.galaxy_r2)) then
        if (all(dist2source.lt.starburst_r2)) then
          if (lrefine(b).lt.sim_starburst_refine_level) then
            refine(b) = .true.
            derefine(b) = .false.
            print *,'ZOOM: starburst ', lrefine(b), sim_starburst_refine_level, dist2source
          else if (lrefine(b).eq.sim_starburst_refine_level) then
            derefine(b) = .false.
          end if
        else
          if (lrefine(b).lt.sim_galaxy_refine_level) then
            refine(b) = .true.
            derefine(b) = .false.
            print *,'ZOOM: galaxy ', lrefine(b), sim_galaxy_refine_level, dist2source
          else if (lrefine(b).eq.sim_galaxy_refine_level) then
            derefine(b) = .false.
          end if
        end if
      end if
      
      if (refine_sn) then
        !!!!!!!! now do supernova zoom
        do n = 1,sim_num_stars

          if (sim_zoom_status(n) .eq. ZOOM_REFINING) then
            zoom_level = sim_zoom_level(n)
            r2 = sim_nova_r2(n)

            bleft = blockCenter - blockSize - sim_star_coords(n,:)
            bright = blockCenter + blockSize - sim_star_coords(n,:)
            
            where (bleft*bright .gt. 0)
              dist2star = min(bleft*bleft, bright*bright)
            elsewhere
              dist2star = 0.0
            end where

            if (all(dist2star.lt.r2)) then
              if (lrefine(b).lt.zoom_level) then
                refine(b)   = .true.
                derefine(b) = .false.
                ready(n) = .false.
                print '(A,I3,A,I3,A,I3,A,I3,A,I3,A)','ZOOM: proc#', sim_mype, ' block#', b, ' zoom_level', lrefine(b), ' of', zoom_level, '(', gr_maxRefine, ')'
              else if (lrefine(b).ge.zoom_level) then
                derefine(b) = .false.
              end if
            end if
          end if

        end do
      end if ! refine_sn
      
      ! restrict max refinement by radius
      ! lmax = upper refinement scale
      ! lmin = lower refinement scale
      ! (lmax - lmin) sets the width of the radial bins
      ! rmin = radius scale of the first jump = N * actual_blockWidth
      ! rmax = sim_box_size / 2
      ! 1) when a block at max refinement passes through the boundary radius for that level
      !    then we can reduce gr_maxRefine
      ! 2) any block whose closest side is past the boundary will be restricted
      
      !radius = max(sqrt(sum(dist2source),1.0)
      !level_ref = floor(sim_derefine_radius_base - log(radius / sim_box_size) / sim_log2) + 1
      
      !level_ref = lmax - 1 - floor(delta_lref * log(radius/rmin) / log_rmax_over_rmin)
      !level_radius = rmin * exp(log_rmax_over_rmin*real(lmax-lrefine(b)-1)/real(delta_lref))
      
      !nref = lrefine(b)+1
      !if (lrefine(b).ge.level_ref .and. nref.ne.sim_derefine_min) then
        !if (refine(b).and.lrefine(b).ge.(gr_maxRefine - 1)) print *,'DEBUG3: Restricting ', lrefine(b), ' >= ', level_ref, ' at ', radius
      !  refine(b) = .false.
      !end if
        
      !if (lrefine(b).gt.level_ref .and. lrefine(b).ge.gr_maxRefine) then
      !  if (allow_step) then
      !    !print *,'DEBUG: Stepping ', lrefine(b), ' > ', level_ref, ' at ', radius
      !    new_step = sim_derefine_step + 1
      !  !else
      !  !  print *,'DEBUG: Suppressing step ', lrefine(b), ' > ', level_ref, ' at ', radius
      !  end if
      !end if
      
    end if ! LEAF
  end do

  ! check if the local processor has finished zooming for supernovae
  where ((sim_zoom_status .eq. ZOOM_REFINING).and. ready)
    sim_zoom_status = ZOOM_LOCAL_READY
  end where

  ! check if we've changed the maximum refinement due to a radius step
  !old_max = sim_derefine_step
  !if (refine_sn) then
  !  sim_derefine_step = 0
  !else
  !  call MPI_ALLREDUCE(new_step, sim_derefine_step, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
  !end if
  
  !if (sim_ismaster) then
  !  if (old_max.ne.sim_derefine_step) then
  !    print *,'REFINE: Making refinement step ', sim_derefine_step
  !  else
  !    print *,'REFINE: Keeping gr_maxRefine ', gr_maxRefine
  !  end if
  !end if

  !old_max = gr_maxRefine
  !max_zoom_level = maxval(sim_zoom_level)
  !gr_maxRefine = max(max_zoom_level - sim_derefine_step, sim_derefine_min)

  call MPI_ALLREDUCE(lnblocks, global_num_blocks, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (sim_ismaster) then
    print *,'REFINE: Block usage = ', global_num_blocks / real(MAXBLOCKS*sim_numPEs) * 100.0
    print *,'REFINE: time until next SN = ', abs(minval(sim_zoom_time,sim_nova_started.eq.0) - time)/3.1557e7
    print *,'REFINE: time since last SN = ', abs(maxval(sim_nova_time,sim_nova_started.eq.1) - time)/3.1557e7
    print *,'REFINE: time from first SN = ', abs(sim_nova_time(1) - time)/3.1557e7

    !if (gr_maxRefine.lt.old_max) print *,'REFINE: Lowering maximum refinement'
    !if (gr_maxRefine.gt.old_max) print *,'REFINE: Raising maximum refinement'
    !print *,'REFINE: max = zoom - step :', gr_maxRefine, max_zoom_level, sim_derefine_step
  end if

  !-------------------------------------------------------------------------------
  
  return
end subroutine sim_markZoom

