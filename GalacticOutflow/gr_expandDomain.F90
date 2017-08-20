!!****if* source/Grid/GridMain/paramesh/gr_expandDomain
!!
!!  NAME
!!     gr_expandDomain
!!
!!  SYNOPSIS
!!     call gr_expandDomain(logical(OUT) :: particlesInitialized)
!!
!!  DESCRIPTION
!!
!!    The grid is initialized in gr_createDomain, with a specified
!!    number of root blocks, typically one single block. This routine
!!    refines appropriate portions of the initialized physical 
!!    domain according to the given refinement criteria, and applies
!!    initial conditions to the AMR domain.
!!
!!    In simulations with particles, under certain conditions particle
!!    positions will also be initialized.  Currently this is the case
!!    if and only if the runtime parameter refine_on_particle_count is
!!    true.
!!
!!  ARGUMENTS
!!    particlesInitialized : is true if this routine initialized particles positions
!!
!!  SIDE EFFECTS
!!
!!    Particle positions may be initialized, see DESCRIPTION above.
!!***

#define DEBUG_PARTICLES

subroutine gr_expandDomain (particlesInitialized)

  use Grid_data, ONLY : gr_domainBC,gr_eosModeInit,gr_refineOnParticleCount, gr_meshMe,&
       gr_meshNumProcs, gr_lrefineMinInit, gr_gcellsUpToDate, gr_meshComm
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY : Logfile_stamp, Logfile_stampVarMask
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_updateRefinement, &
       Grid_getLocalNumBlks, Grid_getListOfBlocks, &
       Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_moveParticles, Grid_sortParticles
  use gr_interface, ONLY : gr_updateRefinement
  use tree, ONLY : lrefine, lrefine_min, lrefine_max, grid_changed, refine, stay, derefine, newchild
  use paramesh_interfaces, ONLY : amr_refine_derefine, amr_restrict
  use Eos_interface, ONLY : Eos_wrapped
  use pt_interface, ONLY : pt_updateTypeDS
  use gr_ptInterface, ONLY : gr_ptMarkRefineDerefine
  use Particles_data, ONLY : pt_numLocal, particles, pt_maxPerProc, pt_indexList, pt_indexCount, pt_typeInfo
  use Kernel_interface, ONLY : Kernel_mapParticles
  !use IO_interface, ONLY : IO_outputInitial
  
  use Simulation_data, ONLY : m_H
  use sim_interface, ONLY : sim_markInitial, sim_markOverdensity, sim_markZoom

#include "Flash.h"

#ifndef FLASH_GRID_PARAMESH2
  use physicaldata, ONLY: no_permanent_guardcells, mpi_pattern_id
#endif
  use Simulation_interface, ONLY : Simulation_initBlock
  use Simulation_data, ONLY : sim_init_refine_done, sim_particle_refine_done
  use Particles_interface, ONLY : Particles_accumCount, &
    Particles_initPositions, Particles_updateGridVar, &
    Particles_updateRefinement
  use Driver_interface, ONLY : Driver_abortFlash
  use RadTrans_interface, ONLY: RadTrans_sumEnergy
  
  use PhotoChem_data, ONLY : phch_default_source_x, phch_default_source_y, phch_default_source_z
  implicit none

#include "constants.h"
#include 'Flash_mpi.h'
#include "Particles.h"

  real, pointer:: solnData(:,:,:,:)
  logical, intent(out) :: particlesInitialized
  integer :: lnblocks, lrefineMinSave


  !!          Local variables and functions

  integer :: ntimes, i, j, local_num

  integer, dimension(2,MDIM) :: blkLimits
  integer, dimension(2,MDIM) :: blkLimitsGC
  integer :: bcount, cur_treedepth, grid_changed_anytime
  logical :: restart = .false.
  logical :: particlesPosnsDone, retainParticles
  integer :: level = FINEST  !not yet implemented, 1 is a dummy value
  integer ,dimension(MAXBLOCKS) :: blkList
  character(len=32), dimension(2,2) :: block_buff
  character(len=32)                 :: int_to_str
  integer :: gridDataStruct, whichBlocks
  logical, parameter :: regrid = .true.
  integer, dimension(MAXBLOCKS, NPART_TYPES) :: particlesPerBlk
  integer :: ierr, num_smoothed, global_smoothed, max_smoothed, refcount, global_refcount
  logical :: keep_going
  integer :: iopt, iempty
  integer, dimension(NPART_TYPES) :: counts

  !!============================================================================



  !!============================================================================

  !! The beginning timestep number, time, and timestep.
  !! If the initial redshift (zinitial) is physical (>= 0),
  !! we use it to initialize the time; otherwise we set the
  !! redshift to zero and get the initial time from tinitial.
  !! The latter case (no cosmology) is the default.

  ! initialize the step counter and the simulation time
  ! the timestep initialization is moved to after the initialization,
  ! so we can check whether it is > t_cfl

  particlesInitialized=.false.

  call gr_initParameshArrays(restart,        &
       gr_domainBC(LOW,IAXIS),gr_domainBC(HIGH,IAXIS), &
       gr_domainBC(LOW,JAXIS),gr_domainBC(HIGH,JAXIS), &
       gr_domainBC(LOW,KAXIS),gr_domainBC(HIGH,KAXIS))

  ! The Paramesh call above may have already resulted in some block refining,
  ! so try get the current max level from the lrefine array. This is only used for
  ! diagnostic output. Note that this assumes that lrefine on the current 
  ! processor is representative of the grid as a whole.
  cur_treedepth = maxval(lrefine)

  gridDataStruct = CENTER
#if NFACE_VARS > 0
  gridDataStruct = CENTER_FACES
#endif

  lrefineMinSave = lrefine_min
  lrefine_min = min(gr_lrefineMinInit,lrefine_max)

  grid_changed_anytime = grid_changed ! save value established by previous Paramesh initialization

  retainParticles=.false.
  particlesPosnsDone = .false.
  
  call gr_updateData()

  !do ntimes = 1, lrefine_max+2
  ntimes = 0
  do while ((.not.particlesPosnsDone).and.ntimes.lt.1000)
     ntimes = ntimes + 1
     if (ntimes .EQ. gr_lrefineMinInit) then
        lrefine_min = lrefineMinSave
     end if
     write (block_buff(1,1), '(a)') 'iteration'
     write (int_to_str, '(i7,a1)') ntimes, ','
     write (block_buff(1,2), '(a,1x,a)') trim(adjustl(int_to_str))
     
     write (block_buff(2,1), '(a)') 'create level'
     write (int_to_str, '(i7)') min(cur_treedepth+1,lrefine_max)
     write (block_buff(2,2), '(a)') trim(adjustl(int_to_str))

     call Logfile_stamp( block_buff, 2, 2, '[GRID gr_expandDomain]')

     call Grid_getLocalNumBlks(lnblocks)
     

#ifndef FLASH_GRID_PARAMESH2
     if (no_permanent_guardcells) then
        call gr_commSetup(gridDataStruct)
     end if
#endif

     if(.not.particlesPosnsDone) then

        if(.not.retainParticles) then
           particlesPosnsDone=.false.
        end if
        sim_particle_refine_done = .false.
        call Particles_initPositions(particlesPosnsDone,retainParticles)
     end if

     call Grid_moveParticles(particles, NPART_PROPS, pt_maxPerProc, pt_numLocal, pt_indexList, pt_indexCount, regrid) 
     call Grid_sortParticles(particles, NPART_PROPS, pt_numLocal, NPART_TYPES, pt_maxPerProc, particlesPerBlk, BLK_PART_PROP, TYPE_PART_PROP)
     call pt_updateTypeDS(particlesPerBlk)

     ! initial zoom into center of grid, not sure why i put it in this function
     if (.not.sim_init_refine_done) call sim_markInitial()

     ! refine on particle smoothing length up to some point, just to get particles on the grid
     !call Particles_updateRefinement(lnblocks)

     grid_changed_anytime = max(grid_changed, grid_changed_anytime)
     grid_changed = 0              ! will be 1 after amr_refine_derefine if the grid actually changed  
     call amr_refine_derefine()
     call gr_updateData()
#ifndef FLASH_GRID_PARAMESH2
     if (grid_changed .NE. 0) mpi_pattern_id = -abs(mpi_pattern_id) !make it different from recognized values
#endif           
     cur_treedepth = max(maxval(lrefine),min(cur_treedepth+1,lrefine_max))

     gr_gcellsUpToDate = .false.
     !end if

  end do !ntimes

  !!!!! START GAS SMOOTHING
  
#ifndef FLASH_GRID_PARAMESH2
     if (no_permanent_guardcells) then
        call gr_commSetup(gridDataStruct)
     end if
#endif

  ! update particle blockIDs and block info, etc
  call Grid_moveParticles(particles, NPART_PROPS, pt_maxPerProc, pt_numLocal, pt_indexList, pt_indexCount, regrid) 
  call Grid_sortParticles(particles, NPART_PROPS, pt_numLocal, NPART_TYPES, pt_maxPerProc, particlesPerBlk, BLK_PART_PROP, TYPE_PART_PROP)
  call pt_updateTypeDS(particlesPerBlk)

  call MPI_ALLREDUCE(pt_typeInfo(PART_LOCAL,:), counts, NPART_TYPES, FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (gr_meshMe.eq.MASTER_PE) print *,'GRIDINIT: counts', sum(counts), counts(:)
  
  ! initial mapping of gas to density/temperature/velocity
  call Kernel_mapParticles(PASSIVE_PART_TYPE)
  
  ! let's dump it to disk here for checking, won't be needed always
  !call IO_outputInitial(0, 0.0)
  
  ! refine on gas density until it is resolved, then map again
  ! repeat until no more refining is done
  keep_going = .true.
  ntimes = 0
  do while (keep_going)
    ntimes = ntimes + 1
    keep_going = .false.  ! only do this again if the grid refined last time
    
    ! keep refining on the mapped density until no more refinements are made
    i = 0
    grid_changed = 1
    global_refcount = -1
    do while (global_refcount.ne.0)
      i = i + 1
      
      call gr_ptFillBlkParticleInfo()
      
      iopt=1; iempty=1
      call amr_restrict(gr_meshMe,iopt,iempty)
      
      refine(:) = .false.
      stay(:) = .false.
      derefine(:) = .false.
      newchild(:) = .false.
      
      call gr_ptMarkRefineDerefine()
      call sim_markOverdensity()
      call sim_markZoom()
      !call gr_markVarThreshold(DENS_VAR, -1.0,  1, -1, 1)
      refcount = count(refine(1:lnblocks))
      call MPI_ALLREDUCE(refcount, global_refcount, 1, FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

      grid_changed = 0
      !call amr_refine_derefine()
      !call gr_updateData()

      call gr_updateRefinement()
      
      grid_changed_anytime = max(grid_changed, grid_changed_anytime) !leave global flag true if grid changed in ANY iteration
      if (global_refcount.ne.0) keep_going = .true.
      if (gr_meshMe.eq.MASTER_PE) print '(A,I3,A,I3,A,I,A,I)','GASREFINE:', ntimes, ' x', i, ' changed=', grid_changed, ' ref=', global_refcount
    end do

#ifndef FLASH_GRID_PARAMESH2
     if (no_permanent_guardcells) then
        call gr_commSetup(gridDataStruct)
     end if
#endif
    
    ! relocate particles and resort particles after refining
    call Grid_moveParticles(particles, NPART_PROPS, pt_maxPerProc, pt_numLocal, pt_indexList, pt_indexCount, regrid) 
    call Grid_sortParticles(particles, NPART_PROPS, pt_numLocal, NPART_TYPES, pt_maxPerProc, particlesPerBlk, BLK_PART_PROP, TYPE_PART_PROP)
    call pt_updateTypeDS(particlesPerBlk)
    
    ! post-refinement remapping
    call Kernel_mapParticles(PASSIVE_PART_TYPE)
    
    ! make another checkpoint for comparison
    !call IO_outputInitial(ntimes, 0.0)
  end do

  ! remove gas particles, code taken from Particles_clean()
  local_num=pt_numLocal
  j=1
  do i=1,local_num
     if(particles(TYPE_PART_PROP,j).eq.real(PASSIVE_PART_TYPE)) then
        particles(:,j)=particles(:,pt_numLocal)
        pt_numLocal=pt_numLocal-1
     else
        j=j+1
     end if
  end do
  
  ! assign a unique id tag to each of the DM particles left
  call pt_createTag()
  
  if (gr_meshMe.eq.MASTER_PE) then
    ! add another particle that we can use for a photoeq source
    ! has TAG_PART_PROP == 0
    pt_numLocal = pt_numLocal + 1
    particles(:, pt_numLocal) = 0.0
    particles(BLK_PART_PROP, pt_numLocal) = 1.0
    particles(TAG_PART_PROP, pt_numLocal) = 0.0
    particles(TYPE_PART_PROP, pt_numLocal) = real(DARK_PART_TYPE)
    particles(POSX_PART_PROP, pt_numLocal) = phch_default_source_x
    particles(POSY_PART_PROP, pt_numLocal) = phch_default_source_y
    particles(POSZ_PART_PROP, pt_numLocal) = phch_default_source_z
    particles(MASS_PART_PROP, pt_numLocal) = m_H
  end if
  
  ! resort after deleting to make sure things are in good order
  call Grid_moveParticles(particles, NPART_PROPS, pt_maxPerProc, pt_numLocal, pt_indexList, pt_indexCount, regrid) 
  call Grid_sortParticles(particles, NPART_PROPS, pt_numLocal, NPART_TYPES, pt_maxPerProc, particlesPerBlk, BLK_PART_PROP, TYPE_PART_PROP)
  call pt_updateTypeDS(particlesPerBlk)
  
  ! initialization for dark matter smoothing
  call smooth_Init()
  
  ! This is here for safety, in case the user did not take care to make things
  ! thermodynamically consistent in the initial state.- KW
  !whichBlocks = LEAF
  !if (ntimes == 1) whichBlocks = ALL_BLKS
  call Grid_getListOfBlocks(LEAF, blkList, bcount)

  call Timers_start("eos")
  do i = 1, bcount
    call Grid_getBlkIndexLimits(blkList(i), blkLimits, blkLimitsGC)
    call Eos_wrapped(MODE_DENS_TEMP,blkLimits,blkList(i))
  end do
  call Timers_stop("eos")

  grid_changed = max(grid_changed, grid_changed_anytime) !leave global flag true if grid changed in ANY iteration

  if(gr_refineOnParticleCount) then
     if(.not.particlesPosnsDone) call Driver_abortFlash(&
       "This distribution of particles will not fit on the grid. Increase pt_maxPerProc, or decrease the particle count.")
     particlesInitialized=.true.
  end if

  lrefine_min = lrefineMinSave

  call gr_ensureValidNeighborInfo(10)

  call MPI_BARRIER(gr_meshComm, ierr)

  return
end subroutine gr_expandDomain
