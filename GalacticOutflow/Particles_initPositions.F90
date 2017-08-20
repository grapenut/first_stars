!!****if* source/Particles/ParticlesInitialization/Particles_initPositions
!!
!! NAME
!!    Particles_initPositions
!!
!! SYNOPSIS
!!
!!   call Particles_initPositions( logical(inout) :: partPosInitialized,
!!                                 logical(out)   :: updateRefine)
!!
!!
!! DESCRIPTION
!!
!!    Initialize particle locations. This routine calls pt_initPositions
!!    which a Particles unit's local routine to initialize the positions
!!    on leaf blocks. The routine also creates tags for all the particles
!!    This routine will initialize based on Lattice or with Density 
!!    distribution depending upon which of the two is selected. 
!!
!! ARGUMENTS
!!
!!  partPosInitialized : boolean indicating whether particles positions were 
!!            successfully initialized.
!!            On entry, a value of TRUE is taken to mean that position
!!            initialization has already been completed previously,
!!            so the routine returns immediately (leaving partPosInitialized
!!            TRUE).
!!            If particles are disabled (as per runtime parameter useParticles),
!!            this implementation also returns immediately with
!!            partPosInitialized set to TRUE.
!!            Otherwise, partPosInitialized will be to TRUE if all particles
!!            have been placed in the domain successfully. A return value
!!            of FALSE may indicate that only some particles have been placed
!!            in the domain, perhaps because of space limitations in the
!!            particles array in some MPI tasks; partially initialized
!!            particles data of this kind may still be useful during FLASH
!!            initialization, in particular for the the purpose of providing
!!            refinement criteria for the initial Grid construction if the
!!            runtime parameter refine_on_particle_count is TRUE, but a
!!            fully successful Particles_initPositions invocation is still
!!            required before the simulation is allowed to proceed with
!!            its main evolution loop.
!!
!!  updateRefine : is set to TRUE if the routine wishes to indicate that during
!!                 the initial iterative construction of an AMR Grid (see
!!                 Grid_initDomain and gr_expandDomain), the initialization of
!!                 particle positions need not be repeated for each iteration if
!!                 all particles have already been placed in the domain in a
!!                 previous iteration.  Under that condition, subsequent calls
!!                 to Particles_initPositions from the Grid construction loop
!!                 will have partPosInitialized=.TRUE.  so will return
!!                 immediately, and Particles_updateRefinement will be called in
!!                 each iteration when the number of blocks or the block
!!                 distribution may have changed, to make sure that the retained
!!                 particles get moved to the correct block (hence the name of
!!                 the dummy argument).
!!
!!                 This implementation always returns FALSE.
!!                 Alternative implementations may wish to return TRUE
!!                 instead if initialization is very expensive.
!!                 
!!
!! NOTES
!!
!!  
!!
!! SIDE EFFECTS
!!
!!  Updates particles data that is private to the Particles unit.
!!
!!  May modify an internal flag (pt_posInitialized) that keeps track of
!!  whether initialization of particle positions is complete.
!!
!! SEE ALSO
!!
!!  Driver_initFlash
!!  Grid_initDomain
!!  Particles_initData
!!***

!!#define DEBUG_PARTICLES

subroutine Particles_initPositions (partPosInitialized,updateRefine)


  use Grid_data, ONLY : gr_lrefineMinInit
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_sortParticles
  use Driver_interface, ONLY : Driver_abortFlash
  use pt_interface, ONLY : pt_initPositions, pt_createTag, pt_initLocal, pt_updateTypeDS
  use Particles_data, ONLY : pt_posInitialized, pt_numLocal, useParticles, pt_velInitialized, &
       pt_typeInfo, particles, pt_meshNumProcs, pt_meshMe, pt_maxPerProc, pt_indexList, pt_indexCount
  
  use Simulation_data, ONLY : sim_baryon_init_file, sim_darkmatter_init_file, sim_baryons_loaded, &
       sim_darkmatter_loaded, sim_init_refine_done, sim_smallrho, sim_particle_max_limit
  use tree, ONLY : lrefine_min
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
 
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "Particles.h"

  logical, INTENT(INOUT) :: partPosInitialized
  logical, INTENT(OUT) :: updateRefine

  real, dimension(2,MDIM):: boundBox
  integer :: blockID
  integer :: blkCount
  integer,dimension(MAXBLOCKS) :: blkList

!----------------------------------------------------------------------

  integer, parameter :: num_baryon_cols = 10
  integer, parameter :: num_dm_cols = 7
  
  integer :: max_particles, total_particles
  integer, save :: baryon_loop_count, darkmatter_loop_count, fileBaryon, fileDM, status, call_counter, load_limit
  integer :: ierr
  integer, dimension(MAXBLOCKS, NPART_TYPES) :: particlesPerBlk

  real :: baryon_array(num_baryon_cols), darkmatter_array(num_dm_cols)
  
  logical :: finishedBaryon, finishedDM
  logical,save :: firstcall = .true.
  logical, parameter:: regrid = .true.
  
  real :: px, py, pz
  
  real, save :: scale, scale2, scale3, redshift
  real, save :: xmin, xmax, ymin, ymax, zmin, zmax
  
  integer, dimension(NPART_TYPES) :: counts
  integer, save :: dest_proc

  character (len=40) :: time_string
  
!----------------------------------------------------------------------

  if(.not.useParticles) then
     partPosInitialized = .true.
  end if
  if(partPosInitialized) return

  ! keep particles between iterations that were previously loaded
  updateRefine = .true.

  !    use master to load X particles and distribute them
  !    keep loading particles and distributing until we reach 90% pt_max
  !    maybe we should retain particles between loops? I think we can just
  !    continue reading using a saved file descriptor

  !open particle file and loop through, loading particles until finished
  if (firstcall) then
    if (pt_meshMe .eq. MASTER_PE) then
      fileBaryon = 143
      open(unit=fileBaryon, file=sim_baryon_init_file, form='formatted', action='read', iostat=status)
      rewind(fileBaryon)
      
      if (status /= 0) print '(2A)', "Particles_initPositions: error opening baryon file: ", sim_baryon_init_file
      
      fileDM = 147
      open(unit=fileDM, file=sim_darkmatter_init_file, form='formatted', action='read', iostat=status)
      rewind(fileDM)

      if (status /= 0) print '(2A)', "Particles_initPositions: error opening dark matter file: ", sim_darkmatter_init_file
      
      dest_proc = MASTER_PE
    end if
    
    baryon_loop_count = 0
    darkmatter_loop_count = 0
    pt_numLocal = 0
    
    call_counter = 0

    call RuntimeParameters_get("xmin", xmin)
    call RuntimeParameters_get("xmax", xmax)
    call RuntimeParameters_get("ymin", ymin)
    call RuntimeParameters_get("ymax", ymax)
    call RuntimeParameters_get("zmin", zmin)
    call RuntimeParameters_get("zmax", zmax)

    call RuntimeParameters_get('zInitial', redshift)
    
    scale = 1.0 / (1.0 + redshift)
    scale2 = scale * scale
    scale3 = scale2 * scale
    
    load_limit = 0.5*sim_particle_max_limit

    if (pt_meshMe.eq.MASTER_PE) print *, "Particles_initPositions: START = ", xmax, ymax, zmax, scale
    
    firstcall = .false.
  end if

  finishedBaryon = sim_baryons_loaded
  finishedDM = sim_darkmatter_loaded
  
  call_counter = call_counter + 1

  ! allow the grid to refine to the minimum first
  if (call_counter.lt.lrefine_min .or. .not.sim_init_refine_done) return
  
  call current_date_time(time_string)
  if (pt_meshMe.eq.MASTER_PE) print *,'Particles_initPositions: call counter == ', call_counter, trim(time_string)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! BARYON INITIALIZATION
  
  ! keep going until the end of the file
  do while (.not.finishedBaryon)
    baryon_loop_count = baryon_loop_count + 1

    ! load particles on master only, until it is full
    do while (pt_meshMe.eq.MASTER_PE .and. pt_numLocal.lt.load_limit)
    
      ! fill particle array from particle file
      read(fileBaryon, *, iostat=status) baryon_array
    
      ! check to see if we reached end of file
      if (status .gt. 0) then
        print *, 'PARTICLES: Gas particle file read error! ', status
        exit
      else if (status .lt. 0) then
        print *,'PARTICLES: Gas particle file done'
        finishedBaryon = .true.
        exit
      end if
      
      px = baryon_array(1) / scale
      if (px < xmin .or. px > xmax) then
        print *, 'SKIPPED: ', baryon_array(1:3) / scale
        cycle
      end if
      
      py = baryon_array(2) / scale
      if (py < ymin .or. py > ymax) then
        print *, 'SKIPPED: ', baryon_array(1:3) / scale
        cycle
      end if

      pz = baryon_array(3) / scale
      if (pz < zmin .or. pz > zmax) then
        print *, 'SKIPPED: ', baryon_array(1:3) / scale
        cycle
      end if
    
      pt_numLocal = pt_numLocal + 1
    
      particles(:, pt_numLocal) = 0.0
      particles(BLK_PART_PROP, pt_numLocal) = 1.0
      particles(TYPE_PART_PROP, pt_numLocal) = PASSIVE_PART_TYPE
      particles(POSX_PART_PROP, pt_numLocal) = baryon_array(1) / scale			! x
      particles(POSY_PART_PROP, pt_numLocal) = baryon_array(2) / scale			! y
      particles(POSZ_PART_PROP, pt_numLocal) = baryon_array(3) / scale			! z
      particles(VELX_PART_PROP, pt_numLocal) = baryon_array(7) * baryon_array(4) / scale	! m*vx
      particles(VELY_PART_PROP, pt_numLocal) = baryon_array(7) * baryon_array(5) / scale	! m*vy
      particles(VELZ_PART_PROP, pt_numLocal) = baryon_array(7) * baryon_array(6) / scale	! m*vz
      particles(MASS_PART_PROP, pt_numLocal) = baryon_array(7)				! mass
      particles(GRID_DENS_PART_PROP, pt_numLocal) = baryon_array(8) * scale3			! density
      particles(SMOOTHTAG_PART_PROP, pt_numLocal) = baryon_array(9) / scale			! smoothing length
      particles(MASS_SAVE_PART_PROP, pt_numLocal) = baryon_array(7) * baryon_array(10) / scale2	! mass * temperature
    end do


    ! if we loaded all particles from the file, this means we are done
    call MPI_BCAST(finishedBaryon, 1, FLASH_LOGICAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    if (finishedBaryon) then
    
      sim_baryons_loaded = .true.
    
      if (pt_meshMe.eq.MASTER_PE) then
        close(fileBaryon)
      end if
    end if

    ! Just copy particles in bulk to the next processor in line and keep chugging away
    ! Call Grid_moveParticles and Grid_sortParticles only after they have all been loaded

    ! increment the destination processor and broadcast it so everybody knows
    dest_proc = dest_proc + 1
    call MPI_BCAST(dest_proc, 1, FLASH_INTEGER, MASTER_PE, MPI_COMM_WORLD, ierr)

    ! send particles from master to next proc in line
    if (pt_meshMe.eq.MASTER_PE) then
      ! send num particles and particles array from master
      call MPI_SEND(pt_numLocal, 1, FLASH_INTEGER, dest_proc, 111, MPI_COMM_WORLD, ierr)
      if (pt_numLocal.gt.0) call MPI_SEND(particles(:,1:pt_numLocal), pt_numLocal*NPART_PROPS, FLASH_REAL, dest_proc, 222, MPI_COMM_WORLD, ierr)
      pt_numLocal = 0
    else if (pt_meshMe.eq.dest_proc) then
      ! receive num particles and particles array on child
      call MPI_RECV(pt_numLocal, 1, FLASH_INTEGER, MASTER_PE, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if (pt_numLocal.gt.0) call MPI_RECV(particles(:,1:pt_numLocal), pt_numLocal*NPART_PROPS, FLASH_REAL, MASTER_PE, 222, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    call MPI_ALLREDUCE(pt_numLocal, max_particles, 1, FLASH_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(pt_numLocal, total_particles, 1, FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call current_date_time(time_string)
    if (pt_meshMe.eq.MASTER_PE) print '(A,A,A,I,A,I,A,I)', "PART: ", trim(time_string), " gas loop ", baryon_loop_count, " loaded ", total_particles, " max ", max_particles

    ! exit early if any procs are already full, pass it on to the smoothing and refinement before coming back here
    ! we will return here as long as partPosInitialized is false.
    !if (max_particles .gt. safety_limit) then
    !  if (pt_meshMe.eq.MASTER_PE) print *,'WARNING: Exceeded safety limit in Particles_initPositions:baryons'
    !  return
    !end if
  end do

  if (.not.finishedBaryon) then
    ! still loading baryons, don't even worry about dark matter yet
    return
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! DARK MATTER INITIALIZATION
  
  ! keep going until the end of the file
  do while (.not.finishedDM)
    darkmatter_loop_count = darkmatter_loop_count + 1

    ! load particles on master only, until it is full
    do while (pt_meshMe.eq.MASTER_PE .and. pt_numLocal.lt.load_limit)
    
      ! fill particle array from particle file
      read(fileDM, *, iostat=status) darkmatter_array
    
      ! check to see if we reached end of file
      if (status .gt. 0) then
        print *, 'PARTICLES: DM particle file read error! ', status
        exit
      else if (status .lt. 0) then
        print *,'PARTICLES: DM particle file done'
        finishedDM = .true.
        exit
      end if
    
      px = darkmatter_array(1) / scale
      if (px < xmin .or. px > xmax) then
        print *, 'SKIPPED: ', darkmatter_array(1:3) / scale
        cycle
      end if
      
      py = darkmatter_array(2) / scale
      if (py < ymin .or. py > ymax) then
        print *, 'SKIPPED: ', darkmatter_array(1:3) / scale
        cycle
      end if

      pz = darkmatter_array(3) / scale
      if (pz < zmin .or. pz > zmax) then
        print *, 'SKIPPED: ', darkmatter_array(1:3) / scale
        cycle
      end if

      pt_numLocal = pt_numLocal + 1
    
      particles(:, pt_numLocal) = 0.0
      particles(BLK_PART_PROP, pt_numLocal) = 1.0
      particles(TYPE_PART_PROP, pt_numLocal) = DARK_PART_TYPE
      particles(POSX_PART_PROP, pt_numLocal) = darkmatter_array(1) / scale
      particles(POSY_PART_PROP, pt_numLocal) = darkmatter_array(2) / scale
      particles(POSZ_PART_PROP, pt_numLocal) = darkmatter_array(3) / scale
      particles(VELX_PART_PROP, pt_numLocal) = darkmatter_array(4) / scale
      particles(VELY_PART_PROP, pt_numLocal) = darkmatter_array(5) / scale
      particles(VELZ_PART_PROP, pt_numLocal) = darkmatter_array(6) / scale
      particles(MASS_PART_PROP, pt_numLocal) = darkmatter_array(7)	!mass
      particles(GRID_DENS_PART_PROP, pt_numLocal) = sim_smallrho	!density
    end do

    ! if we loaded all particles from the file, this means we are done
    call MPI_BCAST(finishedDM, 1, FLASH_LOGICAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    if (finishedDM) then
    
      sim_darkmatter_loaded = .true.
    
      if (pt_meshMe.eq.MASTER_PE) then
        close(fileDM)
      end if
    end if
    
    ! Just copy particles in bulk to the next processor in line and keep chugging away
    ! Call Grid_moveParticles and Grid_sortParticles only after they have all been loaded

    ! increment the destination processor and broadcast it so everybody knows
    dest_proc = dest_proc + 1
    call MPI_BCAST(dest_proc, 1, FLASH_INTEGER, MASTER_PE, MPI_COMM_WORLD, ierr)

    ! send particles from master to next proc in line
    if (pt_meshMe.eq.MASTER_PE) then
      ! send num particles and particles array from master
      call MPI_SEND(pt_numLocal, 1, FLASH_INTEGER, dest_proc, 111, MPI_COMM_WORLD, ierr)
      if (pt_numLocal.gt.0) call MPI_SEND(particles(:,1:pt_numLocal), pt_numLocal*NPART_PROPS, FLASH_REAL, dest_proc, 222, MPI_COMM_WORLD, ierr)
      pt_numLocal = 0
    else if (pt_meshMe.eq.dest_proc) then
      ! receive num particles and particles array on child
      call MPI_RECV(pt_numLocal, 1, FLASH_INTEGER, MASTER_PE, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if (pt_numLocal.gt.0) call MPI_RECV(particles(:,1:pt_numLocal), pt_numLocal*NPART_PROPS, FLASH_REAL, MASTER_PE, 222, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    call MPI_ALLREDUCE(pt_numLocal, max_particles, 1, FLASH_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(pt_numLocal, total_particles, 1, FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call current_date_time(time_string)
    if (pt_meshMe.eq.MASTER_PE) print '(A,A,A,I,A,I,A,I)', "PART: ", trim(time_string), " DM loop ", darkmatter_loop_count, " loaded ", total_particles, " max ", max_particles

    ! exit early if any procs are already full, pass it on to the smoothing and refinement before coming back here
    ! we will return here as long as partPosInitialized is false.
    !if (max_particles .gt. safety_limit) then
    !  if (pt_meshMe.eq.MASTER_PE) print *,'WARNING: Exceeded safety limit in Particles_initPositions:darkmatter'
    !  return
    !end if
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if (finishedBaryon .and. finishedDM) then
    call current_date_time(time_string)
    if (pt_meshMe.eq.MASTER_PE) print "(A,A,A)","PART: ", trim(time_string), " moving particles"

    ! use grid_moveparticles to redistribute particles to their proper positions
    call Grid_moveParticles(particles, NPART_PROPS, pt_maxPerProc, pt_numLocal, pt_indexList, pt_indexCount, regrid) 

    call Grid_sortParticles(particles, NPART_PROPS, pt_numLocal, NPART_TYPES, pt_maxPerProc, particlesPerBlk, BLK_PART_PROP, TYPE_PART_PROP)
    call pt_updateTypeDS(particlesPerBlk)

    partPosInitialized = .true.
    pt_posInitialized = partPosInitialized
    pt_velInitialized = .true.
    
    call current_date_time(time_string)
    if (pt_meshMe.eq.MASTER_PE) print '(A)', "Particles_initPositions: all done ", trim(time_string)
    
  end if

  !call pt_createTag()
  return
  
end subroutine Particles_initPositions

