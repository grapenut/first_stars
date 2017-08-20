!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sim_initNova(n)

  use Simulation_data, ONLY : sim_ismaster, sim_mype, sim_star_coords, &
                              sim_nova_radius, sim_star_mass, sim_nova_yield, sim_nova_temp, sim_nova_velocity, M_Sun,&
                              sim_nova_particles, sim_nova_started, sim_nova_time, sim_derefine_step
  use tree, ONLY: nodetype, lnblocks
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
      Grid_getCellCoords, Grid_putPointData, Grid_getBlkPtr
  use Eos_interface, ONLY : Eos_wrapped
  use Particles_data, ONLY: pt_numLocal, particles, pt_maxPerProc, useParticles
  use Driver_interface, ONLY : Driver_getSimTime, Driver_abortFlash
  !use Cool_data, ONLY : cool_species_map, cool_metal_list, cool_nonmetal_list, cool_num_metals, cool_num_nonmetals
  use Cosmology_interface, ONLY : Cosmology_getRedshift
  
  implicit none

#include "constants.h"
#include "Flash.h"  

  integer, intent(IN) :: n

  integer :: i, j, k, b, s, px, py, pz, np
  
  real, allocatable, dimension(:) :: cell_x, cell_y, cell_z
  real				:: nova_dens, nova_velx, nova_vely, nova_velz
  real				:: radius2, cx, cy, cz, ic, jc, kc, cdx, cdy, cdz, xpos, ypos, zpos
  real				:: cdist2, nova_eint, nova_ener, cell_volume, time
  integer, dimension(2,MDIM) 	:: blkLimits, blkLimitsGC
  integer 			:: sizeX, sizeY, sizeZ, particles_per_axis
  integer, dimension(MDIM) 	:: axis
  !real, DIMENSION(:,:,:,:), POINTER :: solnData
  
  real :: redshift, invscale, invscale2, nova_radius, nova_velocity, nova_temp

  logical :: gcell = .true., doprint = .true., doprint_cpu = .false.
  
  integer :: new_particles = 0
  real :: fraction
  
  sim_nova_started(n) = 1
  sim_derefine_step = 0
  
  call Driver_getSimTime(time)
  sim_nova_time(n) = time
  
  !! put nova quantities in comoving coordinates
  call Cosmology_getRedshift(redshift)
  invscale = (1.0 + redshift)
  invscale2 = invscale*invscale
  
  nova_radius = sim_nova_radius(n)
  nova_dens = 3.0 * sim_star_mass(n) * M_Sun / (4.0 * PI * nova_radius**3.0)
  nova_velocity = sim_nova_velocity(n)
  nova_temp = invscale2 * sim_nova_temp
  radius2 = nova_radius*nova_radius
  
  ic = sim_star_coords(n, 1)
  jc = sim_star_coords(n, 2)
  kc = sim_star_coords(n, 3)

  if (sim_ismaster) print *,'INITNOVA: (r,x,y,z)=', nova_radius, ic, jc, kc

  np = pt_numLocal

  do b = 1,lnblocks
  
    if (nodetype(b).eq.LEAF) then

      !! Loop over the cells in this block and check their distance
          
      ! Setup block limits with/without guard cells
      call Grid_getBlkIndexLimits(b,blkLimits,blkLimitsGC)
      !call Grid_getBlkPtr(b, solnData)

      ! Number of cells w/guardcells
      sizeX = blkLimitsGC(HIGH,IAXIS)
      sizeY = blkLimitsGC(HIGH,JAXIS)
      sizeZ = blkLimitsGC(HIGH,KAXIS)

      ! Allocate coordinate arrays
      allocate(cell_x(sizeX))
      allocate(cell_y(sizeY))
      allocate(cell_z(sizeZ))

      ! Initialize coordinate arrays
      cell_x = 0.0
      cell_y = 0.0
      cell_z = 0.0

      ! Setup coordinate arrays
      call Grid_getCellCoords(IAXIS, b, CENTER, gcell, cell_x, sizeX)
      call Grid_getCellCoords(JAXIS, b, CENTER, gcell, cell_y, sizeY)
      call Grid_getCellCoords(KAXIS, b, CENTER, gcell, cell_z, sizeZ)
      
      cdx = (cell_x(blkLimits(LOW,IAXIS)+1) - cell_x(blkLimits(LOW,IAXIS)))
      cdy = (cell_y(blkLimits(LOW,JAXIS)+1) - cell_y(blkLimits(LOW,JAXIS)))
      cdz = (cell_z(blkLimits(LOW,KAXIS)+1) - cell_z(blkLimits(LOW,KAXIS)))

      cell_volume = cdx * cdy * cdz
      
      !! Can cause integer overflow on large blocks
      !! Calculate particles per axis directly after distance checking
      !!  to make sure we're not bothering to calculate on large blocks
      !!  far outside the supernova's radius
      !particles_per_cell = max(1,floor(sim_nova_particles * (3.0 * cell_volume) / (4.0 * PI * sim_nova_radius**3.0)))
      !particles_per_axis = int(real(particles_per_cell)**(1.0/3.0))
      !particles_per_cell = int(real(particles_per_axis)**3.0)
      
      doprint = .true.
      new_particles = 0
      
      ! Loop over cells in the block.  For each, compute the physical position of 
      ! its left and right edge and its center as well as its physical width.  
      ! Then decide which side of the initial discontinuity it is on and initialize 
      ! the hydro variables appropriately.

!      if (doprint_cpu) then
!        call current_date_time(time_string)
!        print '(A6,I4,I4,A,A)','START ', sim_mype, b, '       - ', time_string
!      end if

      do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        axis(KAXIS) = k
        cz = cell_z(k) - kc
         
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
          axis(JAXIS) = j
          cy = cell_y(j) - jc

          do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
            axis(IAXIS) = i
            cx = cell_x(i) - ic
                 
            ! Find minimum and maximum distance from (ic,jc,kc) for each dimension.
            ! For each coordinate, if both "left" and "right" distances have the same sign,
            ! then the smaller magnitude is the minimum.  Otherwise (ic,jc,kc) is
            ! contained within the interval for that dimension, so the minimum is 0.
            ! The maximum distance is always the larger of the two magnitudes.
            ! Nonexistent dimensions have had all distances set to zero, so they are
            ! ignored.

            cdist2 = cx**2 + cy**2 + cz**2

            if (cdist2 <= radius2) then

              particles_per_axis = max(1,floor(( sim_nova_particles * (3.0 * cell_volume) / (4.0 * PI) )**(1.0/3.0) / nova_radius ))
            
              if (doprint) then
                !call current_date_time(time_string)
                !print '(A6,I4,I4,I6,A,A)','NOVA  ', sim_mype, b, particles_per_axis, ' - ', time_string
                doprint = .false.
                doprint_cpu = .true.
              endif
    
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !!  Supernova initial values  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              nova_velx = cx / nova_radius * nova_velocity
              nova_vely = cy / nova_radius * nova_velocity
              nova_velz = cz / nova_radius * nova_velocity

              call Grid_putPointData(b, CENTER, DENS_VAR, EXTERIOR, axis, nova_dens)
              call Grid_putPointData(b, CENTER, TEMP_VAR, EXTERIOR, axis, nova_temp)   
              call Grid_putPointData(b, CENTER, VELX_VAR, EXTERIOR, axis, nova_velx)
              call Grid_putPointData(b, CENTER, VELY_VAR, EXTERIOR, axis, nova_vely)
              call Grid_putPointData(b, CENTER, VELZ_VAR, EXTERIOR, axis, nova_velz)
              call Grid_putPointData(b, CENTER, ENER_VAR, EXTERIOR, axis, 0.5 * (nova_velx**2.0 + nova_vely**2.0 + nova_velz**2.0))
              
              do s=1,NSPECIES-1
                call Grid_getPointData(b, CENTER, SPECIES_BEGIN - 1 + s, EXTERIOR, axis, fraction)
                call Grid_putPointData(b, CENTER, SPECIES_BEGIN - 1 + s, EXTERIOR, axis, (1.0-sim_nova_yield) * fraction)
              enddo

              call Grid_putPointData(b, CENTER, Z_SPEC, EXTERIOR, axis, sim_nova_yield)
              
              
              if (useParticles .and. np.lt.pt_maxPerProc) then
              
                !! Add new particles evenly in the cell
                do px = 1,particles_per_axis
                  do py = 1,particles_per_axis
                    do pz = 1,particles_per_axis
                    
                      if (np + 1 .lt. pt_maxPerProc) then
                        np = np + 1
                        
                        new_particles = new_particles + 1
                        
                        xpos = cell_x(i) + cdx*(-0.5 + real(px)/real(particles_per_axis))
                        ypos = cell_y(j) + cdy*(-0.5 + real(py)/real(particles_per_axis))
                        zpos = cell_z(k) + cdz*(-0.5 + real(pz)/real(particles_per_axis))
                        
                        !particles(:,2:np) = particles(:,1:np-1)
                        
                        particles(:,np) = 0.0
                        
                        particles(TYPE_PART_PROP,np) = PASSIVE_PART_TYPE

                        particles(POSX_PART_PROP,np) = xpos
                        particles(POSY_PART_PROP,np) = ypos
                        particles(POSZ_PART_PROP,np) = zpos
                        
                        particles( BLK_PART_PROP,np) = real(b)
                        particles( TAG_PART_PROP,np) = real(sim_mype*100000000 + new_particles*10 + n)
                      else
                        print *,'ERROR: exceeded pt_maxPerProc'
                        call Driver_abortFlash('INITNOVA: can not add any more particles');
                      end if !pt_maxPerProc

                    end do
                  end do
                end do

              else
                print *,'ERROR! can not add more particles in initNova!'
                call Driver_abortFlash('INITNOVA: can not add any more particles');
              end if  !! useParticles

            end if  ! dist < radius
          enddo
        enddo
      enddo
      
      pt_numLocal = np

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Set pressures and energy
      call Eos_wrapped(MODE_DENS_TEMP, blkLimits, b)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Set ENER to EINT + 1/2 MV^2
      do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)

       axis(KAXIS) = k
        
       do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)

         axis(JAXIS) = j

          do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!  Supernova initial values  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            call Grid_getPointData(b, CENTER, EINT_VAR, EXTERIOR, axis, nova_eint)   
            call Grid_getPointData(b, CENTER, VELX_VAR, EXTERIOR, axis, nova_velx)
            call Grid_getPointData(b, CENTER, VELY_VAR, EXTERIOR, axis, nova_vely)
            call Grid_getPointData(b, CENTER, VELZ_VAR, EXTERIOR, axis, nova_velz)

            nova_ener = nova_eint + 0.5 * (nova_velx**2.0 + nova_vely**2.0 + nova_velz**2.0)

            call Grid_putPointData(b, CENTER, ENER_VAR, EXTERIOR, axis, nova_ener)
          enddo
        enddo
      enddo

      ! Clean up coordinate arrays
      deallocate(cell_x)
      deallocate(cell_y)
      deallocate(cell_z)

!      if (doprint_cpu) then
!        call current_date_time(time_string)
!        print '(A6,I4,I4,I6,A,A)','DONE  ', sim_mype, b, new_particles, ' - ', time_string
!      end if
    
    end if ! LEAF
  end do !lnblocks
  
!  call current_date_time(time_string)  
!  print '(A,I4,A,A)','---- CPU DONE ---- ', sim_mype, ' - ', time_string
  
end subroutine sim_initNova
