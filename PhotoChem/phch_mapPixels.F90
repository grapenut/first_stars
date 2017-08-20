
subroutine phch_mapPixels(numblocks, blocklist)
  use PhotoChem_data, ONLY : phch_npix, phch_nbin, phch_ismaster, phch_proton_mass, phch_bin_min_cubed, &
                             phch_pixel_number, phch_pixel_volume, phch_stromgren_radius, phch_stromgren_r2, &
                             phch_pixel_density, phch_bin_min, phch_bin_max, phch_bin_radius, phch_current_invscale3, phch_source_x, &
                             phch_source_y, phch_source_z, phch_current_scale2, phch_xi_max, &
                             phch_photon_number, phch_pixel_lambda, phch_current_scale3, &
                             phch_binzero_volume, phch_binzero_number, phch_binzero_lambda, phch_binzero_density, &
                             phch_binzero_photons, phch_pixel_photons, phch_pixel_hplus, &
                             phch_binzero_hplus, phch_total_h_fraction, phch_hplus_fraction, &
                             phch_bin_min_physical, phch_bin_max_physical, phch_bin_min_squared, &
                             phch_bin_radius_physical, phch_bin_r3, phch_bin_r3_physical

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_putPointData, Grid_getPointData
    
  use Cosmology_interface, ONLY : Cosmology_getRedshift
  
  use phch_interface, ONLY : cell2pix, xyz2bin

  implicit none      

#include "constants.h"
#include "Flash.h"

  include "Flash_mpi.h"

  integer, INTENT(in) :: numblocks
  integer, intent(in), dimension(numblocks) :: blocklist

  real :: center_x, center_y, center_z

  integer :: i, j, k, b, blockId, pix, bin
  
  real :: xx, yy, zz, dx, dy, dz
  
  real, allocatable,dimension(:) ::xCenter,yCoord,zCoord
  
  real,  allocatable, dimension(:,:) :: total_pixel_number, total_pixeL_volume, total_pixel_lambda
  
  real :: zone_dens, zone_temp, current_sum, last_sum, absorbed_photons
  real :: r, r2, inside, cell_volume
  
  real :: total_binzero_number, total_binzero_lambda, total_binzero_volume, total_binzero_hplus

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX,sizeY,sizeZ, ierr
  integer, dimension(MDIM) :: axis

  logical :: gcell = .true.
  character(20) :: num, fmt
  character(40) :: time_string
  
  logical :: stromgren_breakout
  
  real :: redshift, scale, invscale, invscale2
  
  real, parameter :: fourpi = 4.0 * PI
  
  real :: local_flux, hydrogen_fraction, hplus_fraction
  
  ! precompute and store scale factors for all calculations so we don't do this in every function
  call Cosmology_getRedshift(redshift)
  invscale = (1.0 + redshift)
  invscale2 = invscale * invscale
  scale = 1.0 / invscale
  phch_current_scale2 = scale * scale
  phch_current_scale3 = phch_current_scale2 * scale
  phch_current_invscale3 = invscale2 * invscale
  
  ! scale bin min/max/widths into comoving coordinates
  phch_bin_min = phch_bin_min_physical * invscale
  phch_bin_max = phch_bin_max_physical * invscale

  phch_bin_min_squared = phch_bin_min * phch_bin_min
  phch_bin_min_cubed = phch_bin_min_squared * phch_bin_min

  phch_bin_radius(:) = phch_bin_radius_physical(:) * invscale
  phch_bin_r3(:) = phch_bin_r3_physical(:) * phch_current_invscale3
  
  ! clear out the aggregate variables
  phch_pixel_number = 0.0
  phch_pixel_volume = 0.0
  phch_pixel_density = 0.0
  phch_pixel_lambda = 0.0
  phch_pixel_photons = 0.0
  phch_pixel_hplus = 0.0
  
  phch_binzero_number = 0.0
  phch_binzero_volume = 0.0
  phch_binzero_density = 0.0
  phch_binzero_lambda = 0.0
  phch_binzero_photons = 0.0
  phch_binzero_hplus = 0.0
  
  phch_stromgren_radius = 0.0
  phch_stromgren_r2 = 0.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Map grid density/temperature to pixel/bin

  call current_date_time(time_string)
  if (phch_ismaster) print '(3A)', "PHCH: ", trim(time_string), " Mapping grid to pixels"

  center_x = phch_source_x
  center_y = phch_source_y
  center_z = phch_source_z

  do b=1,numblocks
    blockId = blocklist(b)

    ! Setup block limits with/without guard cells
    call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

    ! Retrieve the overall size of the block with guard cells
    sizeX = blkLimitsGC(HIGH,IAXIS)
    sizeY = blkLimitsGC(HIGH,JAXIS)
    sizeZ = blkLimitsGC(HIGH,KAXIS)

    ! Allocate coordinate arrays
    allocate(xCenter(sizeX))
    allocate(yCoord(sizeY))
    allocate(zCoord(sizeZ))

    ! Initialize coordinate arrays
    xCenter = 0.0
    yCoord = 0.0
    zCoord = 0.0

    ! Setup coordinate arrays
    if (NDIM == 3) call Grid_getCellCoords(KAXIS, blockId, CENTER,gcell, zCoord, sizeZ)
    if (NDIM >= 2) call Grid_getCellCoords(JAXIS, blockId, CENTER,gcell, yCoord, sizeY)
    call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, sizeX)
    
    ! cell half widths
    dx = 0.5*(xCenter(2) - xCenter(1))
    dy = 0.5*(yCoord(2) - yCoord(1))
    dz = 0.5*(zCoord(2) - zCoord(1))

    do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)

       zz = zCoord(k) - center_z ! get cell center coordinates in  z-direction
       axis(KAXIS) = k
       
       do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)

          axis(JAXIS) = j
          yy = yCoord(j)- center_y ! get cell center coordinates in the y-direction

          do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

             axis(IAXIS) = i
             xx  = xCenter(i) - center_x ! get cell center, left, and right positions in x
             
             call Grid_getPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, zone_dens)
             call Grid_getPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, zone_temp)

             call Grid_getPointData(blockId, CENTER, H_SPEC, EXTERIOR, axis, hydrogen_fraction)
             call Grid_getPointData(blockId, CENTER, HPLU_SPEC, EXTERIOR, axis, hplus_fraction)
             !call Grid_getPointData(blockId, CENTER, HMIN_SPEC, EXTERIOR, axis, hminus_fraction)
             !call Grid_getPointData(blockId, CENTER, HTWO_SPEC, EXTERIOR, axis, h2_fraction)
             !call Grid_getPointData(blockId, CENTER, HTWP_SPEC, EXTERIOR, axis, h2plus_fraction)
             !call Grid_getPointData(blockId, CENTER, HD_SPEC, EXTERIOR, axis, hd_fraction)
             
             phch_total_h_fraction = (hydrogen_fraction + hplus_fraction) / phch_proton_mass
             phch_hplus_fraction = hplus_fraction / phch_proton_mass
             
             call cell2pix(xx, yy, zz, dx, dy, dz, zone_dens, zone_temp)

          enddo
       enddo
    enddo

    ! Clean up coordinate arrays
    deallocate(xCenter)
    deallocate(yCoord)
    deallocate(zCoord)
  
  end do

  ! sum the aggregate variables across all processors
  allocate(total_pixel_number(phch_npix, phch_nbin))
  allocate(total_pixel_volume(phch_npix, phch_nbin))
  allocate(total_pixel_lambda(phch_npix, phch_nbin))

  call MPI_ALLREDUCE(phch_pixel_number, total_pixel_number, phch_npix*phch_nbin, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  phch_pixel_number = total_pixel_number

  call MPI_ALLREDUCE(phch_pixel_hplus, total_pixel_number, phch_npix*phch_nbin, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  phch_pixel_hplus = total_pixel_number

  call MPI_ALLREDUCE(phch_pixel_volume, total_pixel_volume, phch_npix*phch_nbin, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(phch_pixel_lambda, total_pixel_lambda, phch_npix*phch_nbin, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

  ! and the stuff interior to the bins (binzero)
  call MPI_ALLREDUCE(phch_binzero_number, total_binzero_number, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(phch_binzero_hplus, total_binzero_hplus, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(phch_binzero_volume, total_binzero_volume, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(phch_binzero_lambda, total_binzero_lambda, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  
  ! replace the individual local values with the global summed values
  phch_pixel_volume = total_pixel_volume
  phch_pixel_lambda = total_pixel_lambda
  
  phch_binzero_number = total_binzero_number
  phch_binzero_volume = total_binzero_volume
  phch_binzero_lambda = total_binzero_lambda
  phch_binzero_hplus = total_binzero_hplus

  ! divide by total volume for the volume-weighted sums in binzero
  if (phch_binzero_number.gt.0.0 .and. phch_binzero_volume.gt.0.0) then
    phch_binzero_density = phch_binzero_number / phch_binzero_volume
    phch_binzero_lambda = phch_binzero_lambda / phch_binzero_volume
  else
    phch_binzero_density = 0.0
    phch_binzero_lambda = 0.0
  end if

  ! do the same for the rest of the bins
  ! oh yeah! we can just save pixel_lambda as pixel_photons right here!
  where (phch_pixel_volume.gt.0.0)
    phch_pixel_photons = phch_pixel_lambda 
    phch_pixel_lambda = phch_pixel_lambda / phch_pixel_volume
    phch_pixel_density = phch_pixel_number / phch_pixel_volume
  elsewhere
    phch_pixel_photons = 0.0
    phch_pixel_lambda = 0.0
    phch_pixel_density = 0.0
  end where
  
  !do pix=1,phch_npix
  !  do bin=1,phch_nbin
  !    if (phch_pixel_volume(pix,bin).gt.0.0) then
  !      phch_pixel_density(pix,bin) = phch_pixel_number(pix,bin) / phch_pixel_volume(pix,bin)
  !      phch_pixel_lambda(pix,bin) = phch_pixel_lambda(pix,bin) / phch_pixel_volume(pix,bin)
  !    else
  !      phch_pixel_density(pix,bin) = 0.0
  !      phch_pixel_lambda(pix,bin) = 0.0
  !    end if
  !  end do
  !end do
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Calculate Stromgren radius for each pixel

  call current_date_time(time_string)
  if (phch_ismaster) print '(3A)', "PHCH: ", trim(time_string), " Calculating Stromgren radius"

  ! calculate what would be the number of photons required to ionize all of binzero
  phch_binzero_photons = fourpi / 3.0 * phch_binzero_lambda * phch_bin_min_cubed
  
  ! skip the stromgren calculation if there are not even enough photons to ionize binzero
  if (phch_binzero_photons.gt.phch_photon_number) then
    ! set the minimum stromgren radius and reduce phch_binzero_photons to the number emitted
    phch_stromgren_radius = phch_bin_min
    phch_stromgren_r2 = phch_bin_min_squared

    if (phch_ismaster) print *,'PHCH: stromgren radius is less than phch_bin_min!'
    
    stromgren_breakout = .false.

  else
  
    stromgren_breakout = .true.

    ! integrate outward along each pixel to count the photons absorbed
    do pix=1,phch_npix
      current_sum = phch_binzero_photons

      ! sum each bin's photon absorption count until the sum equals the total photon number emitted
      do bin=1,phch_nbin
        ! skip empty bins
        if (phch_pixel_density(pix,bin).le.0.0) then
          ! unless it's the last one
          if (bin.eq.phch_nbin) then
            phch_stromgren_radius(pix) = phch_bin_max
            if (phch_ismaster) print *,'PHCH: stromgren radius exceeds phch_bin_max on pixel ', pix
            exit
          end if
          cycle
        end if

        last_sum = current_sum
        current_sum = phch_binzero_photons + sum(phch_pixel_photons(pix,1:bin))

        ! check if the stromgren radius is within this bin
        if (current_sum.gt.phch_photon_number) then
          phch_stromgren_radius(pix) = ( (phch_photon_number - last_sum) * 3.0 / (fourpi * phch_pixel_lambda(pix,bin)) + phch_bin_r3(bin) )**(1.0/3.0)
          exit
        end if
        
        ! check if we failed to find the stromgren radius
        if (bin.eq.phch_nbin) then
          phch_stromgren_radius(pix) = phch_bin_max
          if (phch_ismaster) print *,'PHCH: stromgren radius exceeds phch_bin_max on pixel ', pix
          exit
        end if

      enddo
      
      ! setup Rs^2 for faster comparison to block positions later
      phch_stromgren_r2(pix) = phch_stromgren_radius(pix) * phch_stromgren_radius(pix)
    enddo
  end if

  deallocate(total_pixel_number)
  deallocate(total_pixel_volume)
  deallocate(total_pixel_lambda)
  
  if (phch_ismaster) then
    !! write out to a file
    open(9, file="stromgren_radius.txt", form='formatted', action='write', position='append')
    write(num,*) phch_npix
    write(fmt,*) '(', trim(num), 'e17.8)'
    write(9,fmt) (phch_stromgren_radius(i), i=1,phch_npix)
    close(9)
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Calculate grid flux from pixel/bin data

  call current_date_time(time_string)
  if (phch_ismaster) print '(3A)', "PHCH: ", trim(time_string), " Mapping pixels to grid"

  !! Loop through all cells again and add up the value of FLUX_VAR if inside stromgren radius
  do b=1,numblocks
    blockId = blocklist(b)

    call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

    ! Retrieve the overall size of the block with guard cells
    sizeX = blkLimitsGC(HIGH,IAXIS)
    sizeY = blkLimitsGC(HIGH,JAXIS)
    sizeZ = blkLimitsGC(HIGH,KAXIS)

    ! Allocate coordinate arrays
    allocate(xCenter(sizeX))
    allocate(yCoord(sizeY))
    allocate(zCoord(sizeZ))

    ! Initialize coordinate arrays
    xCenter = 0.0
    yCoord = 0.0
    zCoord = 0.0

    ! Setup coordinate arrays
    if (NDIM == 3) call Grid_getCellCoords(KAXIS, blockId, CENTER,gcell, zCoord, sizeZ)
    if (NDIM >= 2) call Grid_getCellCoords(JAXIS, blockId, CENTER,gcell, yCoord, sizeY)
    call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, sizeX)
    
    dx = (xCenter(2) - xCenter(1))
    dy = (yCoord(2) - yCoord(1))
    dz = (zCoord(2) - zCoord(1))

    do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)

       zz = zCoord(k) - center_z ! get cell center coordinates in  z-direction
       axis(KAXIS) = k
       
       do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)

          axis(JAXIS) = j
          yy = yCoord(j) - center_y ! get cell center coordinates in the y-direction

          do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

             axis(IAXIS) = i
             xx  = xCenter(i) - center_x ! get cell center, left, and right positions in x
             
             cell_volume = dx*dy*dz
             
             call xyz2bin(xx, yy, zz, pix, bin)
             r2 = xx*xx + yy*yy + zz*zz
             
             local_flux = 0.0
             
             if (r2.le.phch_bin_min_squared) then
               ! we're inside binzero
               
               ! check to make sure the cell isn't on the point
               r = 0.5 * dx
               if (r2.le.r*r) then
                 ! set radius to half the cell width and square it
                 ! this just keeps it from blowing up and makes it smooth since
                 ! the next cell over will be approximately 1.5 * dx and 2.5 * dx, etc
                 r2 = r * r
               end if
               
               ! let's just set the flux as if all the radiation was hitting each layer within binzero
               ! this will ensure that it is always heated enough to drive out the gas
               ! and should be smooth in the tranisition from binzero to bin 1
               if (stromgren_breakout) then
                 local_flux = phch_photon_number * invscale2 / (fourpi * r2)
               else
                 local_flux = phch_xi_max
               end if

             else if (r2.le.phch_stromgren_r2(pix)) then
               ! within the stromgren sphere proper, calculate flux
               ! total photon number minus the photons absorbed in interior bins
               ! should precompute the attenuation from interior bins and just calculate
               ! additional attention from the edge of the bin to the cell

               ! number of photons absorbed in this particular bin, between in the interior edge and r
               ! r3 = r2**1.5
               inside = fourpi / 3.0 * phch_pixel_lambda(pix,bin) * (r2**1.5 - phch_bin_r3(bin))
               
               ! sum the photons absorbed interior to this cell
               if (bin.eq.1) then
                 absorbed_photons = phch_binzero_photons + inside
               else
                 absorbed_photons = phch_binzero_photons + sum(phch_pixel_photons(pix,1:bin-1)) + inside
               end if

               ! the local flux is the total photon number minus absorbed photon number divided by 4pi r^2
               local_flux = (phch_photon_number - absorbed_photons) * invscale2 / (fourpi * r2)
             end if
             
             ! put the flux back on the grid
             call Grid_putPointData(blockId, CENTER, FLUX_VAR, EXTERIOR, axis, local_flux)
          enddo
       enddo
    enddo

    ! Clean up coordinate arrays
    deallocate(xCenter)
    deallocate(yCoord)
    deallocate(zCoord)
  
  end do

  return

end subroutine phch_mapPixels

