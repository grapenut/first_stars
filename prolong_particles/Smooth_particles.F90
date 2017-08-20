module Smooth_particles


  
  contains
    
    
    subroutine SmoothParticleDensity(var)


      use smooth_Data
      use Particles_data, ONLY : particles, pt_numLocal, pt_posAttrib
      use Grid_data, ONLY : gr_meshMe, gr_meshComm
      use tree, ONLY : lnblocks, grid_changed
      use physicaldata, ONLY: unk
      use Particles_interface, ONLY: Particles_updateGridVar
      use Grid_data, ONLY : gr_meshComm
      use Grid_interface, ONLY : Grid_getLocalNumBlks, &
           Grid_getListOfBlocks, Grid_getBlkPhysicalSize,  & 
           Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr, & 
           Grid_getBlkPhysicalSize, Grid_fillGuardCells, Grid_mapMeshToParticles
      use gr_interface, ONLY : gr_checkGridState
      
      implicit none

      
#include "constants.h"  
#include "Flash.h"
#include "Flash_mpi.h"
#include "Particles.h"
      
      integer, intent(IN) :: var
      
      logical, save :: first_call = .true.
      
      integer :: globalCount, localCount
      integer :: ierr, i, m
      character (len=40) :: time_string
      
      integer :: numAttrib
      integer,dimension(2,1) :: attrib
      integer :: mapType = WEIGHTED
      
      real :: total_mass
      logical :: debug = .false.
      
      integer :: upper_lrefine_global, lower_lrefine_global

      real :: part_mass_smooth, part_mass_smooth_global
      integer, dimension(MDIM), save :: posProperties
      logical :: gridChanged


      if(first_call) then
         call smooth_Init()
         
         posProperties(1) = posx
         posProperties(2) = posy
         posProperties(3) = posz
         
         first_call = .false.
      endif
      
      if(.not. SmoothParticles) return

      
      ! record node types and calculate maxMeshRefineLevel
      call smooth_recordNodeTypeStates()

    
      ! particles need to know density of their host cell:
      numAttrib = 1
      attrib(GRID_DS_IND,numAttrib)=DENS_VAR
      attrib(PART_DS_IND,numAttrib)=GRID_DENS_PART_PROP
      call Grid_mapMeshToParticles(particles(:,1:pt_numlocal),NPART_PROPS,BLK_PART_PROP,& 
           pt_numLocal,pt_posAttrib,numAttrib,attrib,mapType,CENTER)
      
      call CreateSmoothParticles()

      localCount = num_smooth
      call MPI_ALLREDUCE(localCount,globalCount,1,FLASH_INTEGER,MPI_SUM,gr_meshComm,ierr)
      
      if(globalCount .eq. 0) then
         if(MASTER_PE .eq. gr_meshMe) then 
            print*, "No particle smoothing required!"
            print*, "max refine on grid, lrefine_smooth min/max=", maxMeshRefineLevel, lrefine_smooth_min, lrefine_smooth_max
         endif
         do i = 1, lnblocks
            unk(var,:,:,:,i) = unk(PDEN_VAR,:,:,:,i) 
         end do
         return
      endif
      
      if(MASTER_PE .eq. gr_meshMe) then 
         call current_date_time(time_string)
         if(gr_meshMe .eq. MASTER_PE) print '(3A)', "SMOOTH: ", trim(time_string), "  Preparing to smooth particle density..."
      endif
      
      
      ! zero out incoming variable:
      do i = 1, lnblocks
         unk(var,:,:,:,i) = 0.0
         !unk(VARI_VAR,:,:,:,i) = 0.0
      end do

      !! Particles which are being smoothed have had their
      !! mass zeroed out, so the following only maps
      !! non-smoothed particles to AUX_VAR density field
      call Particles_updateGridVar(MASS_PART_PROP,AUX_VAR)
      
      ! get smoothing particles to lrefine=lrefine_smooth blocks
      call MoveSmoothParticles(particles_smooth,num_smooth_part_props, &
           max_num_smooth_part,num_smooth)
      
      ! smoothing particles can exist on any block with refinemenet level
      ! between lrefine_smooth_min and lrefine_smooth_max inclusive

      ! map the smoothing particles' mass to the coarser levels
      call mpi_allreduce(lower_lrefine, lower_lrefine_global, 1, MPI_INTEGER, &
           MPI_MIN, gr_meshComm, ierr)
      call mpi_allreduce(upper_lrefine, upper_lrefine_global, 1, MPI_INTEGER, & 
           MPI_MAX, gr_meshComm, ierr)
      if (gr_meshMe .eq. MASTER_PE ) print*, "lower and upper lrefine global = ", lower_lrefine_global, upper_lrefine_global

      do m = lower_lrefine_global, upper_lrefine_global
        if (gr_meshMe.eq.MASTER_PE) print *, "SMOOTH: mapping level ", m
        call smooth_prepareNodes(m, .true.)
         
        !! amr_mg_init needed if in non-cautious paramesh mode
        !! ctss: also, need to save the value of grid_changed since amr_mg_init will reset it.
        !! this is only needed (as far as I know) for proper functioning of the
        !! pfft enhanced multigrid solver.
        !gridChanged = grid_changed
        !call amr_mg_init()
        !grid_changed = gridChanged
        
        call gr_commSetup(-1)
        call gr_checkGridState()
        
        call smooth_mapParticlesToMesh(particles_smooth,num_smooth_part_props, &
             num_smooth,max_num_smooth_part,mass,var,1,m,blk_coarse, posProperties)
      enddo
      

      !! At this point, all smooth particles should be mapped to blocks with
      !! refinement level between lrefine_smooth_min, and lrefine_smooth_max
   

      ! let's print out global mass of smooth particles
      ! and what was mapped to the grid
      if(debug) then
         call smooth_findTotalMass(var,total_mass)
         part_mass_smooth = 0.0
         if (num_smooth .gt. 0) then
            do i = 1, num_smooth
               part_mass_smooth  = part_mass_smooth + particles_smooth(mass,i)
            enddo
         endif
         call MPI_ALLREDUCE(part_mass_smooth,part_mass_smooth_global,1,FLASH_REAL,MPI_SUM,MPI_COMM_WORLD,ierr) 
         call MPI_ALLREDUCE(num_smooth,globalCount,1,FLASH_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr) 
         if(gr_meshMe .eq. MASTER_PE ) then
              print*, "global number smoothing particles after moving and mapping=", globalCount
              print*, "total mass of smoothing particles = ", part_mass_smooth_global
              print*, "total mass mapped onto grid = ", total_mass
           endif
           do m=lower_lrefine_global, upper_lrefine_global
              part_mass_smooth = 0.0
              do i=1, num_smooth
                 if(particles_smooth(lref_dest,i) .eq. real(m)) then
                    part_mass_smooth = part_mass_smooth + particles_smooth(mass,i) 
                 endif
              enddo
              call MPI_ALLREDUCE(part_mass_smooth,part_mass_smooth_global,1,FLASH_REAL,MPI_SUM,MPI_COMM_WORLD,ierr) 
                  if(gr_meshMe .eq. MASTER_PE ) print*, "lref_dest, tot mass smoothing part=", m, part_mass_smooth_global
           enddo
      endif
 
      
      call gr_freeCommRecvBuffer()
      
      !! Prolong this data down to leaf blocks:
      !! Passes .true. to prepareNodes and restoreNodeTypes for
      !! cautious mode, although this may only be needed for the
      !! smooth_mapParticlesToMesh call
      if(gr_meshMe .eq. MASTER_PE) then
         print*, "Prolonging data from ", lrefine_smooth_min, " to ", maxMeshRefineLevel
      endif
      do m = lrefine_smooth_min, maxMeshRefineLevel-1
         call smooth_prepareNodes(m, .true.)
         call smooth_fillGuardCells(m, var)
         call smooth_Prolong(m, var, var, 2)
      enddo

      call gr_freeCommRecvBuffer()
      
      ! restore node type states
      call smooth_restoreNodeTypes(.true.)
      
      
      ! add together the two contributions
      do i = 1, lnblocks
         unk(var,:,:,:,i) = unk(var,:,:,:,i) + unk(AUX_VAR,:,:,:,i)
         !unk(VARI_VAR,:,:,:,i) = unk(VARI_VAR,:,:,:,i) + unk(AUX_VAR,:,:,:,i)
      end do

      ! restore particle masses which were zeroed out
      call smooth_restoreParticleMasses()
      
      ! make sure mass is conserved between PDEN_VAR and PDE_VAR
      call smooth_checkMassConservation(LEAF)


      if(MASTER_PE .eq. gr_meshMe) then 
         call current_date_time(time_string)
         if(gr_meshMe .eq. MASTER_PE) print '(3A)', "SMOOTH: ", trim(time_string), "  Particle density smoothed"
      endif
      
  
      
      return
    end subroutine SmoothParticleDensity
    
  
  

    
    subroutine GetSmoothParticleAcceleration()
      
      use smooth_Data
      use Grid_data, ONLY : gr_meshMe, gr_meshComm
      use Gravity_interface, ONLY : Gravity_accelListOfBlocks
      use Grid_interface, ONLY : Grid_getListOfBlocks, &
           Grid_fillGuardCells, Grid_mapMeshToParticles
      use Particles_data, ONLY : particles, pt_numLocal
      use tree, ONLY : grid_changed
      
      implicit none
      
#include "constants.h"  
#include "Flash.h"
#include "Flash_mpi.h"
#include "Particles.h"
      
      integer :: localCount, globalCount, ierr
      integer :: blockCount, m
      integer :: numAttrib
      integer,dimension(2,1) :: attrib
      integer,dimension(MAXBLOCKS) :: blockList
      integer :: mapType = WEIGHTED
      character (len=40) :: time_string
      integer :: gridChanged
   
      integer :: lower_lrefine_global
      
      logical, save :: first_call = .true.
      
      if(first_call) then
         call smooth_Init()
         first_call = .false.
      endif
      
      if(.not. SmoothParticles) return

      call smooth_recordNodeTypeStates()

      
      call CreateSmoothParticles()

   
      localCount = num_smooth
      call MPI_ALLREDUCE(localCount,globalCount,1,FLASH_INTEGER,MPI_SUM,gr_meshComm,ierr)
      

      if(globalCount .eq. 0) then
         if(MASTER_PE .eq. gr_meshMe) print*, "No grid acceleration needed"
         return
      endif

      
      if(MASTER_PE .eq. gr_meshMe) then 
         call current_date_time(time_string)
         if(gr_meshMe .eq. MASTER_PE) print '(3A)', "SMOOTH: ", trim(time_string), "  Preparing to get smoothed particle acceleration"
      endif
      
      !! amr_mg_init needed if in non-cautious paramesh mode
      !! ctss: also, need to save the value of grid_changed since amr_mg_init will reset it.
      !! this is only needed (as far as I know) for proper functioning of the
      !! pfft enhanced multigrid solver.
      gridChanged = grid_changed
      call amr_mg_init()
      grid_changed = gridChanged
      
      
      call MoveSmoothParticles(particles_smooth, num_smooth_part_props, &
           max_num_smooth_part, num_smooth)
      
      
      !!!!
      ! Need to make sure higher levels have proper GPOT_VAR data
      ! to this end, here we restrict GPOT from leaf blocks and fill guard cell data
      
      !! See gr_hgSolve for basic strategy
      !! Passes .false. to prepareNodes and restoreNodeTypes for
      !! non-cautious mode
      
      call mpi_allreduce(lower_lrefine, lower_lrefine_global, 1, MPI_INTEGER, &
           MPI_MIN, gr_meshComm, ierr)

      call gr_freeCommRecvBuffer()
      
      ! I need to restrict GPOT
      if(gr_meshMe .eq. MASTER_PE) then
         print*, "Restricting GPOT from ", maxMeshRefineLevel, " to ", lrefine_smooth_min
      endif
      do m=maxMeshRefineLevel, lrefine_smooth_min+1, -1
         call smooth_Restrict(m,GPOT_VAR,GPOT_VAR)
      enddo

      !! After restriction, fill the guard cells, in low -> high order
      !! needed for particle acceleration
      if(gr_meshMe .eq. MASTER_PE) then
         print*, "Filling GPOT guardcells on level ", lrefine_smooth_min, " to ", maxMeshRefineLevel-1
      endif
      do m = lrefine_smooth_min, maxMeshRefineLevel-1
         call smooth_prepareNodes(m, .false.)
         call smooth_fillGuardCells(m, GPOT_VAR)
      enddo
      call gr_freeCommRecvBuffer()

      ! restore node type states
      call smooth_restoreNodeTypes(.false.)
  
      !!!!

      call Grid_getListOfBlocks(ALL_BLKS,blockList,blockCount)
      
      dont_skip_smoothed_particles = .true.
      numAttrib = 1
      attrib(GRID_DS_IND,numAttrib)=GRAC_VAR

      attrib(PART_DS_IND,numAttrib)=accx
      call Gravity_accelListOfBlocks(blockCount,blockList,IAXIS,GRAC_VAR)
      call Grid_mapMeshToParticles(particles_smooth(:,1:num_smooth),num_smooth_part_props,blk_coarse,&
           num_smooth,smooth_posAttrib,numAttrib,attrib,mapType)

      attrib(PART_DS_IND,numAttrib)=accy
      call Gravity_accelListOfBlocks(blockCount,blockList,JAXIS,GRAC_VAR)
      call Grid_mapMeshToParticles(particles_smooth(:,1:num_smooth),num_smooth_part_props,blk_coarse,&
           num_smooth,smooth_posAttrib,numAttrib,attrib,mapType)
      
      attrib(PART_DS_IND,numAttrib)=accz
      call Gravity_accelListOfBlocks(blockCount,blockList,KAXIS,GRAC_VAR)
      call Grid_mapMeshToParticles(particles_smooth(:,1:num_smooth),num_smooth_part_props,blk_coarse,&
           num_smooth,smooth_posAttrib,numAttrib,attrib,mapType)
   
      dont_skip_smoothed_particles = .false.
      
      call TransferSmoothedAccelerations(particles(:,1:pt_numLocal),NPART_PROPS,pt_numLocal,&
           particles_smooth,num_smooth_part_props, num_smooth,max_num_smooth_part)


      ! restore particle masses
      call smooth_restoreParticleMasses()

      
      if(MASTER_PE .eq. gr_meshMe) then 
         call current_date_time(time_string)
         if(gr_meshMe .eq. MASTER_PE) print '(3A)', "SMOOTH: ", trim(time_string), "  Got smoothed particle accelerations"
      endif


      
      return
      
    end subroutine GetSmoothParticleAcceleration
    
    





    ! It's a bit silly I'm writing an (essentially) identical routine
    ! to make this happen. Oh well.
    subroutine smooth_mapParticlesToMesh(particles,part_props,numParticles,&
         maxParticlesPerProc,propPart,varGrid,mode,level,blkProp,posProp)
      
      use smooth_Data, ONLY : posx, posz
      use gr_ptData, ONLY : gr_ptBlkList, gr_ptBlkCount, gr_ptBuf
      use gr_ptMapData, ONLY : gr_ptDomain, NUMBGUARDREGIONS, gr_ptSmearLen
      use Timers_interface, ONLY : Timers_start, Timers_stop
      use Logfile_interface, ONLY: Logfile_stampMessage
      use Driver_interface, ONLY : Driver_abortFlash
      use Grid_data, ONLY : gr_meshMe
      use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
           Grid_getBlkIndexLimits, Grid_getListOfBlocks, Grid_sortParticles, &
           Grid_getBlkCornerID, Grid_getDeltas, Grid_getBlkBoundBox
      use gr_ptInterface, ONLY :  gr_ptStoreOffBlockCells, gr_ptSameProcMap, &
           gr_ptOffProcMap, gr_ptMoveMappedData,gr_ptApplyBCsOneBlk
      use gr_interface, ONLY : gr_checkGridState
      use Particles_interface, ONLY : Particles_mapToMeshOneBlk 
      use tree, ONLY : lrefine
      
      implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "gr_ptMapToMesh.h"

  integer,intent(IN) :: maxParticlesPerProc,numParticles, part_props
  real,dimension(part_props,maxParticlesPerProc),intent(INOUT) :: particles
  integer, INTENT(in) :: propPart, varGrid
  integer, INTENT(in) :: mode, level, blkProp
  integer, dimension(MDIM), intent(IN) :: posProp
  

  integer,parameter :: particleTypes=1
  integer,dimension(MAXBLOCKS,particleTypes) :: particlesPerBlk, particlesPerBlk_level
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer :: blockID,pEnd, vmode, maplevel, blkProperty
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC,srcCoords,destCoords
  integer,dimension(MDIM) :: blkSize,blkSizeGC,guard
  real, allocatable, dimension(:) :: sendBuf, recvBuf
  integer,dimension(BLKID:REFLEVELDIF):: negh
  integer,dimension(MDIM) :: neghCornerID
  integer :: sendBufPtr, blkNo, numNegh, n, sendCount
  integer :: BufferSize, regionIter, ib, jb, kb, ie, je, ke, error
  integer :: localNumParticles
  logical :: mapByLevel
  integer, dimension(MDIM) :: posAttrib
  
  integer, dimension(MDIM) :: intPos
  real,dimension(MDIM) :: deltas,delInv,hFunc
  real,dimension(LOW:HIGH,MDIM) :: bndBox
  real :: dVolInv, particleAttributeValue
  logical, parameter :: fcBdry = .false.

  localNumParticles=numParticles

  vmode = mode

  maplevel = level
  mapByLevel = .true.
  
  blkProperty = blkProp
  
  posAttrib(1) = posProp(1)
  posAttrib(2) = posProp(2)
  posAttrib(3) = posProp(3)
    
  !print *, "Processor", gr_meshMe, "in Grid_mapParticlesToMesh"

  !! We are doing this because in paramesh all blocks have the
  !! same dimension, and this is paramesh specific implementation
  !! so this calculation can be kept out of the loop.
  blockID=1
  call Grid_getBlkIndexLimits(blockID,blkLimits, blkLimitsGC)

  guard=blkLimits(LOW,:)-blkLimitsGC(LOW,:)
  blkSizeGC=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
  blkSize=blkLimits(HIGH,:)-blkLimits(LOW,:)+1


  !! All the blocks are zeroed initially, rather than one at 
  !! a time within the blkNo loop. This is to prevent solnVec being
  !! zeroed when it contains valid values from an earlier iteration 
  !! of blkNo.
  if(mapByLevel) then
     call Grid_getListOfBlocks(REFINEMENT,gr_ptBlkList,gr_ptBlkCount,refinementLevel=mapLevel)
  else
     call Grid_getListOfBlocks(LEAF,gr_ptBlkList,gr_ptBlkCount)
  endif

  !if (gr_ptBlkCount.lt.1) then
  !  print *,'DEBUG: No block at level ', mapLevel, ' on processor ', gr_meshMe
  !  return
  !end if

#ifdef DEBUG_GRIDMAPPARTICLES_VERBOSE
  print *, "SMOOTHDEBUG Processor", gr_meshMe, "LEAF min refinement:", &
       minval(lrefine(gr_ptBlkList(1:gr_ptBlkCount))), "LEAF max refinement:", &
       maxval(lrefine(gr_ptBlkList(1:gr_ptBlkCount))), "no.particles:", numParticles
#endif

  do blkNo = 1, gr_ptblkCount, 1
     call Grid_getBlkPtr(gr_ptblkList(blkNo),solnVec,CENTER)
     if(vmode /= 0) then
        solnVec(varGrid,blkLimitsGC(LOW,IAXIS):blkLimitsGC(LOW,IAXIS)+guard(IAXIS)-1,:,:) = 0.0
        solnVec(varGrid,blkLimitsGC(HIGH,IAXIS)-guard(IAXIS)+1:blkLimitsGC(HIGH,IAXIS),:,:) = 0.0

        if(NDIM >= 2) then
           solnVec(varGrid,:,blkLimitsGC(LOW,JAXIS):blkLimitsGC(LOW,JAXIS)+guard(JAXIS)-1,:) = 0.0
           solnVec(varGrid,:,blkLimitsGC(HIGH,JAXIS)-guard(JAXIS)+1:blkLimitsGC(HIGH,JAXIS),:) = 0.0
        end if

        if(NDIM == 3) then
           solnVec(varGrid,:,:,blkLimitsGC(LOW,KAXIS):blkLimitsGC(LOW,KAXIS)+guard(KAXIS)-1) = 0.0
           solnVec(varGrid,:,:,blkLimitsGC(HIGH,KAXIS)-guard(KAXIS)+1:blkLimitsGC(HIGH,KAXIS)) = 0.0
        end if
     else
        solnVec(varGrid,:,:,:) = 0.0
     end if
     call Grid_releaseBlkPtr(gr_ptblkList(blkNo),solnVec,CENTER)
  end do



  ! The Grid_mapParticlesToMesh routine requires that
  ! the particles are sorted in block order
  call Grid_sortParticles(particles, part_props, localNumParticles, &
       particleTypes, maxParticlesPerProc, particlesPerBlk, blkProperty)
       
  

  call gr_ensureValidNeighborInfo(10) ! We want valid grid information as after a guardcell-filling operation
  call gr_checkGridState()


  !---------------------------------------------------------------------------------------
  ! In order to determine the size of the send / receive buffer, we will count the number 
  ! of cells that exist off-processor.  
  ! Whilst doing this, we will keep track of all the neighbors in a temporary data structure
  ! named gr_ptDomain.
  !---------------------------------------------------------------------------------------
  if (gr_ptSmearLen > 0 ) then
     allocate(gr_ptDomain(gr_ptBlkCount), STAT=error)
     if (error /= 0) then
        call Driver_abortFlash("Severe error. Memory cannot be allocated!")
     end if
     
     particlesPerBlk_level = 0
     particlesPerBlk_level(gr_ptBlkList(1:gr_ptBlkCount),:) = particlesPerBlk(gr_ptBlkList(1:gr_ptBlkCount),:)
     !print '(A,5I)','DEBUGOFF: ', gr_meshMe, gr_ptBlkCount, count(particlesPerBlk(:,:).gt.0), count(particlesPerBlk_level(:,:).gt.0), numParticles
     
     !if (gr_meshMe.eq.10) then
     !  print *,'DEBUGOFF: bcount=', gr_ptBlkCount
     !  print *,'DEBUGOFF: blists=', gr_ptBlkList(1:gr_ptBlkCount)
     !  print *,'DEBUGOFF: perblk=', particlesPerBlk(gr_ptBlkList(1:gr_ptBlkCount),:)
     !end if

     !Subroutine modifies module data structure named gr_ptDomain.
     !! CHECK particlesPerBlk matches gr_ptBlkList ordering?
     call gr_ptStoreOffBlockCells(particlesPerBlk, gr_ptBlkList, gr_ptBlkCount, &
          blkLimitsGC, blkSize, guard, BufferSize)

     !Each process has the same value for BufferSize (global maximum buffer 
     !size).  If it is zero then no data needs communicating.  Set it to 1 
     !so there is buffer space for gr_ptMoveMappedData to send / recv a single 
     !NULL element.
     if (BufferSize <= 0) then
        BufferSize = 1
        if (gr_meshMe == 0) then
           print *, "[Grid_mapParticlesToMesh]: No off-processor mapping required"
        endif
     end if  


     allocate(sendBuf(BufferSize), recvBuf(BufferSize), STAT=error)
     if (error /= 0) then
        call Driver_abortFlash("Severe error. Memory cannot be allocated!")
     end if
  end if



  ! ---------------------------------------------------------------------------------
  ! Now we can actually map the particles to the mesh.
  ! ---------------------------------------------------------------------------------
  allocate(gr_ptBuf(blkSizeGC(IAXIS),blkSizeGC(JAXIS),blkSizeGC(KAXIS)), STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("Severe error. Memory cannot be allocated!")
  end if


  sendBufPtr=1  !Pointer to the first element of sendBuf.
  sendCount=0
  pEnd=0


  do blkNo = 1, gr_ptBlkCount

     blockID = gr_ptBlkList(blkNo)
     if (particlesPerBlk(blockID,particleTypes) > 0) then


        call Grid_getBlkPtr(blockID,solnVec,CENTER)
        call Grid_getDeltas(blockID,deltas)
        call Grid_getBlkBoundBox(blockID,bndBox)
        delInv = 1.0
        delInv(1:NDIM)=1.0/deltas(1:NDIM) 
        intPos=1
        hfunc = 0.0
        dvolInv = delInv(IAXIS)*delInv(JAXIS)*delInv(KAXIS)

        gr_ptBuf = 0.0

        
        do n=1, numParticles
        
           if(int(particles(blkProperty,n)) .eq. blockID) then
              ! this particle should be mapped to this block 
              ! ala Particles_mapToMeshOneBlk:
              
              hfunc(1:NDIM) = (particles(posx:posz,n)-&
                   bndBox(LOW,1:NDIM))*delInv(1:NDIM)
              intPos(1:NDIM)= floor(hfunc(1:NDIM)) + 1 + guard(1:NDIM)
              hfunc(1:NDIM) = modulo(hfunc(1:NDIM),1.0) - 0.5
              particleAttributeValue=particles(propPart,n)*dvolInv
              call pt_mapOneParticle(blkLimitsGC,intPos,particleAttributeValue, &
                   fcBdry,hfunc,gr_ptBuf)
              
           endif
           
        enddo

        
    
        !Apply boundary conditions (BCs) to the guard cells of this block.
        !If for example, reflecting BCs are in place, a mass accumulated 
        !guard cell may map onto an internal cell of the SAME block.  This 
        !is NOT concerned with the global block configuration.  That is 
        !handled in gr_ptSameProcMap & gr_ptOffProcMap, where by periodic 
        !BCs will map a guard cell onto an internal cell of a DIFFERENT block.

        call gr_ptApplyBCsOneBlk(blkLimits,blkLimitsGC,blockID)

        !! All smearing within the block copied in 
        !! We have to add the gr_ptBuf contribution to what is already 
        !! there because an earlier iteration of blkNo could have written
        !! to this block.
        !! solnVec(varGrid,:,:,:) = solnVec(varGrid,:,:,:) + gr_ptBuf(:,:,:)
        !! Actually, only need internal region...
        solnVec(varGrid,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &
             solnVec(varGrid,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) + &
             gr_ptBuf(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))

        call Grid_releaseBlkPtr(blockID,solnVec,CENTER)



        if (gr_ptSmearLen > 0) then

           !! For each section of the gcells find the negh. If the negh is
           !! on the local processor at the same reflevel, copy the 
           !! values into the block. If the negh is there, but at different
           !! reflevel restrict or prolong the values
           !! If the negh is not on the processor, put it in the sendBuf
           !! Even in sendbuf, there can be two ways of putting info.
           !! The destination might be known before sending, or it may 
           !! not be known, both situations have to be handled.

           ! Initialise certain data for the case when NDIM < MDIM.
           ! Required because gr_ptDomain data structure only contains data for NDIM.
           srcCoords(LOW:HIGH, NDIM:MDIM) = 1
           destCoords(LOW:HIGH, NDIM:MDIM) = 1
           neghCornerID(NDIM:MDIM) = 0

           ! Loop over each guard cell region in a block.
           do regionIter = 1, NUMBGUARDREGIONS

              numNegh = gr_ptDomain(blkNo) % haloRegion(regionIter) % numNegh

              if (numNegh > 0) then
                 do n = 1, numNegh

                    negh(:) = gr_ptDomain(blkNo) % haloRegion(regionIter) % neighbor(n) % negh(:)
                    neghCornerID(1:NDIM) = &
                         gr_ptDomain(blkNo) % haloRegion(regionIter) % neighbor(n) % cornerID(1:NDIM)
                    srcCoords(LOW:HIGH,1:NDIM) = &
                         gr_ptDomain(blkNo) % haloRegion(regionIter) % neighbor(n) % srcCoords(LOW:HIGH,1:NDIM)
                    destCoords(LOW:HIGH,1:NDIM) = & 
                         gr_ptDomain(blkNo) % haloRegion(regionIter) % neighbor(n) % destCoords(LOW:HIGH,1:NDIM)

                    ib=srcCoords(LOW,IAXIS); ie=srcCoords(HIGH,IAXIS)
                    jb=srcCoords(LOW,JAXIS); je=srcCoords(HIGH,JAXIS)
                    kb=srcCoords(LOW,KAXIS); ke=srcCoords(HIGH,KAXIS)

                    !Move onto the next guard cell section if the current sections contains no data.
                    !Assuming each guard cell region contains at least some data, then
                    !at least one processor will use all of the space in sendBuf.
                    if(maxVal(abs(gr_ptBuf(ib:ie,jb:je,kb:ke)))/=0.0) then

                       if(negh(BLKPROC)==gr_meshMe) then
                          call gr_ptSameProcMap(srcCoords,destCoords,negh,varGrid)
                       else
                          !print *, "Processor", gr_meshMe, "calling gr_ptOffProcMap.  CornerID=", &
                          !neghCornerID(:,n)
                          call gr_ptOffProcMap(srcCoords,destCoords,BufferSize,&
                               sendBuf,sendCount,sendBufPtr,negh,neghCornerID)
                          !print *, "Processor", gr_meshMe, "packed", sendCount, & 
                          !"reals out of", BufferSize, "space, which will be sent to processor", negh(BLKPROC,n)
                       end if
                    end if

                 end do
              end if
           end do

        end if  !End -> if (gr_ptSmearLen > 0).
        
     end if !Test value in particlesPerBlk
  end do !End loop over blocks.


       
  deallocate(gr_ptBuf)

  if (gr_ptSmearLen > 0) then

     if(sendCount > BufferSize) then
        call Driver_abortFlash("Severe error. Communication buffer too small!!!!")
     end if
     
     call gr_ptMoveMappedData(varGrid,BufferSize,sendBuf,sendCount,recvBuf)


     deallocate(gr_ptDomain)  !Used in gr_ptMoveMappedData when it calls gr_ptDumpState.
     deallocate(sendBuf)
     deallocate(recvBuf)
  end if


      
      return
      
    end subroutine smooth_mapParticlesToMesh
    


  subroutine smooth_mapMeshToParticles (particles, part_props,part_blkID,&
                                        numParticles,posAttrib,&
                                        numAttrib, attrib,&
                                        mapType,gridDataStruct)

    use Driver_interface, ONLY : Driver_abortFlash
    use Grid_interface, ONLY : Grid_getBlkPtr,Grid_releaseBlkPtr,&
         Grid_getLocalNumBlks, Grid_getBlkBoundBox, Grid_getDeltas
    use Particles_interface, ONLY : Particles_mapFromMesh

    implicit none

#include "Flash.h"
#include "constants.h"
#include "GridParticles.h"
#include "Particles.h"


    integer, INTENT(in) :: part_props, numParticles, numAttrib, part_blkID
    real, INTENT(inout),dimension(part_props,numParticles) :: particles
    integer,dimension(MDIM), intent(IN) :: posAttrib
    integer, dimension(PART_ATTR_DS_SIZE,numAttrib),INTENT(in) :: attrib
    integer, INTENT(IN) :: mapType
    integer, optional, intent(IN) :: gridDataStruct

    integer :: i,j,currentBlk,blkCount, prevBlk

    real, pointer, dimension(:,:,:,:) :: solnVec
    real, dimension(LOW:HIGH,MDIM) :: bndBox
    real,dimension(MDIM) :: delta, pos
    integer :: gDataStruct
    real, dimension(numAttrib) :: partAttribVec

    logical :: smooth
    integer :: smooth_part_prop

    if(present(gridDataStruct)) then
       gDataStruct=gridDataStruct
    else
       gDataStruct=CENTER
    end if

    smooth = .true.
    
    if (smooth) then
      smooth_part_prop = 1
    else
      smooth_part_prop = SMOOTHTAG_PART_PROP
    endif
    
    if(numParticles>0) then

       call Grid_getLocalNumBlks(blkCount)

       currentBlk=int(particles(part_blkID,1))
       prevBlk=currentBlk
       call Grid_getBlkPtr(currentBlk,solnVec,gDataStruct)
       call Grid_getBlkBoundBox(currentBlk,bndBox)
       call Grid_getDeltas(currentBlk,delta)
       do i = 1, numParticles
#ifdef DEBUG_GRIDPARTICLES
          if((particles(part_blkID, i) < 0) .or. (particles(part_blkID, i) > blkCount)) then
             call Driver_abortFlash("BLK_PART_PROP out of bounds")
          end if
#endif
          
          if((particles(smooth_part_prop,i) .eq. 1.0) .and. (.not. smooth))  then
            cycle
          end if
          
          currentBlk=int(particles(part_blkID,i))
          if(currentBlk /= prevBlk)then
             call Grid_releaseBlkPtr(prevBlk,solnVec,gridDataStruct)
             call Grid_getBlkPtr(currentBlk,solnVec,gridDataStruct)
             call Grid_getBlkBoundBox(currentBlk,bndBox)
             call Grid_getDeltas(currentBlk,delta)
          end if
          do j = 1,MDIM
             pos(j)=particles(posAttrib(j),i)
          end do

          call Particles_mapFromMesh (mapType, numAttrib, attrib,&
               pos, bndBox,delta,solnVec, partAttribVec)
          do j = 1,numAttrib
             !if((particles(smooth_part_prop,i) .eq. 1.0) .and. (.not. smooth)) then
             !   ! do nothing
             !else
             !   particles(attrib(PART_DS_IND,j),i)=partAttribVec(j)
             !endif
             !if (smooth .or. (particles(smooth_part_prop,i).ne.1.0)) then
               particles(attrib(PART_DS_IND,j),i)=partAttribVec(j)
             !end if
          end do
          prevBlk=currentBlk
       enddo
       call Grid_releaseBlkPtr(currentBlk,solnVec,gridDataStruct)
    end if

  end subroutine smooth_mapMeshToParticles




  subroutine CreateSingleSmoothParticle(particle_data, part_props, lrefine_smooth, particle_id)
    
    use smooth_Data
    use Grid_data, ONLY : gr_meshMe
 
    implicit none

#include "constants.h"  
#include "Flash.h"
#include "Flash_mpi.h"

   
    integer, intent(IN) :: part_props, particle_id, lrefine_smooth
    real, intent(IN), dimension(part_props) :: particle_data
   

    num_smooth = num_smooth + 1
    
    !! Get max and min refinement levels where smooth particles will exist
    if(lrefine_smooth .gt. upper_lrefine) upper_lrefine = lrefine_smooth
    if(lrefine_smooth .lt. lower_lrefine) lower_lrefine = lrefine_smooth

    !! From original particle
    particles_smooth(blk_og,num_smooth)     = particle_data(BLK_PART_PROP)
    particles_smooth(proc_og,num_smooth)    = real(gr_meshMe)
    particles_smooth(posx,num_smooth)       = particle_data(POSX_PART_PROP)
    particles_smooth(posy,num_smooth)       = particle_data(POSY_PART_PROP)
    particles_smooth(posz,num_smooth)       = particle_data(POSZ_PART_PROP)
    particles_smooth(mass,num_smooth)       = particle_data(MASS_PART_PROP)
    particles_smooth(pid,num_smooth)        = real(particle_id)
    particles_smooth(lref_dest,num_smooth)  = real(lrefine_smooth)
    
    !! Unknown for now:
    particles_smooth(accx,num_smooth) = 0.0
    particles_smooth(accx,num_smooth) = 0.0
    particles_smooth(accx,num_smooth) = 0.0
    particles_smooth(blk_coarse,num_smooth) = 0.0
    particles_smooth(proc_coarse,num_smooth) = 0.0
    
    
    return
    
  end subroutine CreateSingleSmoothParticle
  
  



  subroutine FindRefinementLevelDest(mass_dm,grid_dens,lrefine_smooth)

    
    use smooth_Data
    use Simulation_data, ONLY : sim_box_size
    
    implicit none
    
#include "constants.h"  
#include "Flash.h"
    
    real, intent(IN) :: mass_dm, grid_dens
    integer, intent(OUT) :: lrefine_smooth
    
    real :: rs
    
    real, save :: onethird, log2, logfactor
    
    logical, save :: first_call = .true.
    
    
    if(first_call) then
       logfactor = sim_box_size / real(NXB)
       onethird = 1.0/3.0
       log2 = log10(2.0)
       first_call = .false.
    endif


    ! rs = comoving smoothing radius
    rs = smooth_radius_factor * (mass_dm / grid_dens)**onethird
    
    ! assumes box is uniform in all directions and NXB=NYB=NZB
    lrefine_smooth = floor(log10(logfactor / rs) / log2 + 1.0)
    
    if(lrefine_smooth .lt. lrefine_smooth_min) lrefine_smooth = lrefine_smooth_min
    if(lrefine_smooth .gt. lrefine_smooth_max) lrefine_smooth = lrefine_smooth_max
    
    return
    
  end subroutine FindRefinementLevelDest





  subroutine CreateSmoothParticles()
    
    use smooth_Data
    use tree, ONLY : lrefine
    use Grid_data, ONLY : gr_meshMe
    use Particles_data, ONLY : particles, pt_numLocal
    use Driver_interface, ONLY : Driver_abortFlash
 
    implicit none

#include "constants.h"  
#include "Flash.h"
#include "Flash_mpi.h"
    
    integer :: i, blockID, num_smooth_global, ierr

    integer :: lrefine_smooth
    
    logical :: vocal = .true.
    
    real :: darkType
    

    particles_smooth(:,:) = 0.0
    num_smooth = 0
    
    upper_lrefine = 0
    lower_lrefine = 100
    
    darkType = real(DARK_PART_TYPE)
    
    if(pt_numLocal .gt. 0) then
  
       particles(MASS_SAVE_PART_PROP,:) = particles(MASS_PART_PROP,:)
       particles(SMOOTHTAG_PART_PROP,:) = 0.0
       
       do i = 1, pt_numLocal
          ! only do this for dark matter! maybe better to use the sorted typeDS? -jr
          if (particles(TYPE_PART_PROP,i).ne.darkType) cycle
          
          blockID = int(particles(BLK_PART_PROP,i))
          particles(PROC_PART_PROP,i) = real(gr_meshMe)
          
          if(lrefine(blockID) .gt. lrefine_smooth_min) then
             
             ! if block is within range lrefine_smooth_min+1 and lrefine_smooth_max, 
             ! we may create a smoothing particle 
             
             call FindRefinementLevelDest(particles(MASS_PART_PROP,i), & 
                  particles(GRID_DENS_PART_PROP,i),lrefine_smooth)
             
             if(lrefine_smooth .lt. lrefine(blockID)) then
                call CreateSingleSmoothParticle(particles(:,i),NPART_PROPS,lrefine_smooth,i)
                particles(SMOOTHTAG_PART_PROP,i) = 1.0
                particles(MASS_PART_PROP,i) = 0.0
             endif
             
          endif
          
          
          !if(lrefine(blockID) .gt. lrefine_smooth_max) then            
          !  call CreateSingleSmoothParticle(particles(:,i),NPART_PROPS,lrefine_smooth_max,i)
          !   particles(SMOOTHTAG_PART_PROP,i) = 1.0
          !   particles(MASS_PART_PROP,i) = 0.0
          !endif
          
          
       enddo
    endif
    
    particles_smooth(:,num_smooth+1:max_num_smooth_part) = NONEXISTENT
    
    if(vocal) then
       call MPI_ALLREDUCE(num_smooth,num_smooth_global,1,FLASH_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr) 
       if(gr_meshMe .eq. MASTER_PE .and. num_smooth_global .gt. 0) &
            print*, "Global number of smoothing particles created = ", num_smooth_global
    endif
 
 
    return
    
  end subroutine CreateSmoothParticles






  subroutine smooth_restoreParticleMasses()

    use Particles_data, ONLY : particles, pt_numLocal
 
    implicit none

#include "constants.h"  
#include "Flash.h"
    
    if(pt_numLocal .gt. 0) &
    particles(MASS_PART_PROP,:) = particles(MASS_SAVE_PART_PROP,:)
 
    return
    
  end subroutine smooth_restoreParticleMasses







  subroutine MoveSmoothParticles(dataBuf, propCount, maxCount, localCount)
    
    use smooth_Data
    use Grid_interface, ONLY : Grid_getListOfBlocks 
    use Grid_data, ONLY : gr_meshMe
    use Driver_interface, ONLY : Driver_abortFlash
 
    implicit none

#include "constants.h"  
#include "Flash.h"
#include "Flash_mpi.h"
    
    integer, intent(IN) :: maxCount, propCount
    integer, intent(INOUT) :: localCount
    real, dimension(propCount, maxCount), intent(INOUT) :: dataBuf
    !real, allocatable, dimension(:,:) :: destBuf, sourceBuf

    
    integer :: m, numDest, i

    ! Want to move all particles in particles_smooth
    ! array to the block (and processor) with 
    ! lrefine = lsmooth, and which (obviously)
    ! also contains the particle
    
    ! see Grid_moveParticles and gr_ptMoveSieve

    
    allocate(destBuf(propCount,maxCount))
    allocate(sourceBuf(propCount,maxCount))
    
    
    call Grid_getListOfBlocks(ALL_BLKS,smooth_blkList,smooth_blkCount)
    
    dataBuf(blk_coarse,localCount+1:maxCount) = NONEXISTENT
    dataBuf(blk_coarse,1:localCount) = UNKNOWN
    sourceBuf(:,1:localCount) = dataBuf(:,1:localCount)
    
    m = 0
    numDest = 0
    
    do i = 1, localCount
       if(sourceBuf(posx,i) .lt. 1.0) &
            print*, "2. probably bad smooth particle position x", sourceBuf(:,i)
       if(sourceBuf(posy,i) .lt. 1.0) &
            print*, "2. probably bad smooth particle position y", sourceBuf(:,i)
       if(sourceBuf(posz,i) .lt. 1.0) &
            print*, "2. probably bad smooth particle position z", sourceBuf(:,i)
    enddo
    
    call smooth_localMatch(dataBuf,m,propCount,maxCount,sourceBuf,localCount,&
         destBuf,numDest)
    
    localCount = m
    
    dataBuf(blk_coarse,localCount+1:maxCount) = NONEXISTENT
    
    call smooth_moveSieve(dataBuf,localCount,propCount,&
         maxCount,numDest)

    ! set correct processor (probably not needed)
    dataBuf(proc_coarse,1:localCount) = real(gr_meshMe)
    
    ! particles_smooth should now be in the correct place. Let's verify this: 
    call smooth_verifySmoothParticleLocations(dataBuf,propCount,localCount,maxCount)


    deallocate(destBuf)
    deallocate(sourceBuf)
    

    return
    
  end subroutine MoveSmoothParticles
  




  
  subroutine smooth_localMatch(dataBuf,localCount,propCount,maxCount,&
       sourceBuf,numSource,destBuf,numDest)
    
    ! just like gr_ptLocalMatch, but instead of matching a leaf block
    ! this searches ALL blocks for lrefine = lrefine_smooth
    
    use smooth_Data, ONLY : smooth_blkList, smooth_blkCount, blk_coarse, &
         posx, posy, posz, lref_dest, smooth_saveNodeType
    use tree, ONLY : lrefine, nodetype
    use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_outsideBoundBox
    use Driver_interface, ONLY : Driver_abortFlash
    
    implicit none
    
#include "constants.h"  
#include "Flash.h"
#include "Flash_mpi.h"
    
    integer, intent(IN) :: propCount
    integer, intent(IN) :: numSource
    integer, intent(OUT) :: numDest
    integer, intent(INOUT) :: localCount
    integer, intent(IN) :: maxCount
    real,dimension(propCount,maxCount),intent(INOUT) :: dataBuf
    real,dimension(propCount,numSource),intent(INOUT) :: sourceBuf,destBuf
    
    integer ::   blockID, blockID_save
    integer :: i,m,j, lrefine_smooth
    real,dimension(MDIM)::pos
    
    logical :: found, outside
    real,dimension(LOW:HIGH,MDIM) :: bndBox
    integer,dimension(MDIM) :: Negh

    m = localCount
    numDest = 0
    
    do i = 1, numSource
 
       blockID = int(sourceBuf(blk_coarse,i))
       blockID_save = blockID
       pos(IAXIS) = sourceBuf(posx,i)
       pos(JAXIS) = sourceBuf(posy,i)
       pos(KAXIS) = sourceBuf(posz,i)
       lrefine_smooth = int(sourceBuf(lref_dest,i))
       
       ! loop over all local blocks with lrefine = lrefine_smooth
       
       !!!!!!!!!!!!!!!!!
       ! following segment is essentially gr_findBlock:
       found = .false.
       j = 0
       do while((.not. found) .and. (j .lt. smooth_blkCount))
          j=j+1
          blockID=smooth_blkList(j)
          if(lrefine(blockID) .eq. lrefine_smooth) then
             call Grid_getBlkBoundBox(blockID,bndBox)
             call Grid_outsideBoundBox(pos,bndBox,outside,Negh)
             found=.not.outside
          endif
       enddo
       !if(found) print*, "found a particle's locale!",i, blockID, int(sourceBuf(blk_og,i))
       if(.not.found) blockID = NONEXISTENT
       !!!!!!!!!!!!!!!!!
       
       if(blockID .eq. NONEXISTENT) then
          !if(numSource .gt. 500) then
          !   print*, "lone particle found NONEXISTENT", blockID, pos, gr_meshMe, i, blockID_save
          !end if
          numDest = numDest + 1
          destBuf(:,numDest) = sourceBuf(:,i)
       else
          if(m .gt. maxCount) then
             print*, "too many particles in smooth_localMatch", m
             call Driver_abortFlash("too many particles in smooth_localMatch")
          endif
          m=m+1
          dataBuf(:,m) = sourceBuf(:,i)
          dataBuf(blk_coarse,m) = blockID


          if(smooth_saveNodeType(blockID) .eq. LEAF) then
             print*, "smooth_localMatch: shouldn't be sending smooth particles to leaf blocks! wtf"
             print*, "source=", sourceBuf(:,i)
             print*, "blockID=", blockID
             print*, "save node type of blockid=", smooth_saveNodeType(blockID) 
             print*, "current nodetype of blockid=", nodetype(blockID)
             call abort(1)
          endif


       endif
       
    enddo
    
    localCount = m
    
    return
    
  end subroutine smooth_localMatch





  
  subroutine smooth_moveSieve(dataBuf,localCount,propCount,&
      maxCount,numDest)
    
    use smooth_Data
    use gr_ptSieveInterface, ONLY : gr_ptResetProcPair, gr_ptNextProcPair
    use gr_ptData, ONLY : gr_ptSieveCheckFreq, gr_ptSieveFreq
    
    
    implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

    integer,intent(INOUT) :: localCount
    integer,intent(IN) :: propCount
    integer,intent(IN) :: maxCount
    integer,intent(INOUT) :: numDest
    real, dimension(propCount, maxCount),intent(INOUT) :: dataBuf
    
    integer,dimension(MPI_STATUS_SIZE) :: status
    logical :: stillProcessing, mustCommunicate
    real :: nothingReal
    integer :: timesInLoop,  ierr
    integer :: numSource, src, dest
    integer :: sendCount, recvCount,sendTag,recvTag

    nothingReal = -1.0
    gr_ptSieveFreq=gr_ptSieveCheckFreq
    
    stillProcessing = numDest .gt. 0
    
    if(.not. stillProcessing) destBuf(1,1) = nothingReal
    
    call gr_ptResetProcPair(stillProcessing,mustCommunicate)
    
    timesInLoop = 0
    
    do while(mustCommunicate)
       
       timesInLoop = timesInLoop + 1
       
       sendTag = 1
       recvTag = 1
       sendCount = max(1,numDest*propCount)
       recvCount = maxCount*propCount
       
       call gr_ptNextProcPair(timesInLoop, src, dest, stillProcessing, mustCommunicate)
       if(mustCommunicate) then
          call MPI_SENDRECV(destBuf(1,1), sendCount, FLASH_REAL, &
               dest, sendTag, sourceBuf(1,1), &
               recvCount, FLASH_REAL, src, &
               recvTag, FLASH_COMM, status, ierr)
          
          call MPI_GET_COUNT(status,FLASH_REAL,numSource,ierr)
          
          !print*, "proc, localCount, times in loop=", gr_globalMe, localCount, timesInLoop
          
          if(numSource .eq. 1) then
             numDest = 0 
             destBuf(1,1) = nothingReal
             stillProcessing = .false.
          else
             numSource = numSource / propCount
             call smooth_localMatch(dataBuf,localCount,propCount,maxCount,&
                  sourceBuf,numSource,destBuf,numDest)
             stillProcessing = numDest .gt. 0
          endif
       
       endif  ! must communicate
       
    enddo  ! master do loop
    
    return

  end subroutine smooth_moveSieve
  
  
  


  
  subroutine smooth_verifySmoothParticleLocations(dataBuf,propCount,localCount,maxCount)
     
    use smooth_Data
    use tree, ONLY : lrefine
    use Driver_interface, ONLY : Driver_abortFlash
    use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_outsideBoundBox
    use Grid_data, ONLY : gr_meshMe

    implicit none

#include "constants.h"
#include "Flash.h"
    
    integer,intent(INOUT) :: localCount
    integer,intent(IN) :: propCount
    integer,intent(IN) :: maxCount
    real, dimension(propCount, maxCount),intent(INOUT) :: dataBuf   
    
    real, dimension(LOW:HIGH, MDIM) :: bndBox
    
    integer :: m, blockID, refinementLevel

    do m = 1, localCount

       blockID = int(dataBuf(blk_coarse, m))
       refinementLevel = int(dataBuf(lref_dest,m))
       
       if(blockID .eq. NONEXISTENT) then
          print*, "smooth particle = ", dataBuf(:,m)
          call Driver_abortFlash("smooth particle on NONEXISTENT block")
       endif
       
       if(lrefine(blockID) .ne. refinementLevel) then
          print*, "smooth particle = ", dataBuf(:,m)
          call Driver_abortFlash("smooth particle not on proper refinement level")
       endif
       
       if(smooth_saveNodeType(blockID) .eq. LEAF) then
          print*, "smoothing particle got mapped to a leaf block! this shouldn't happen..."
          print*, dataBuf(:,m)
          print*, "nodetype of coarse block = ", smooth_saveNodeType(blockID)
          print*, "this proc =", gr_meshMe
          call Driver_abortFlash("smooth particle mapped to leaf block")
       endif


       call Grid_getBlkBoundBox(blockID,bndBox)
       

        if(dataBuf(posx,m) < bndBox(LOW,IAXIS) .or. &
             dataBuf(posx,m) >= bndBox(HIGH,IAXIS) .or. &
             dataBuf(posy,m) < bndBox(LOW,JAXIS) .or. &
             dataBuf(posy,m) >= bndBox(HIGH,JAXIS) .or. &
             dataBuf(posz,m) < bndBox(LOW,KAXIS) .or. &
             dataBuf(posz,m) >= bndBox(HIGH,KAXIS)) then
           
           
           print*, "smooth particle = ", dataBuf(:,m)
           print*, 'bounding box=', bndBox(:,:)
           
           
           call Driver_abortFlash("smooth_verifyParticleLocations: particle in wrong location!")
           
        endif
       
       
    enddo

    return
    
  end subroutine smooth_verifySmoothParticleLocations
    
    
  





  subroutine TransferSmoothedAccelerations(particles,propCount,localCount,dataBuf,& 
       propCount_smooth,localCount_smooth,maxCount_smooth)
    
    ! take the particles_smooth array, send it back to its original processor, and give its
    ! accelerations to the original DM particles
    ! see gr_ptMove
    
    
    use smooth_Data
    use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs, gr_meshComm
    use Driver_interface, ONLY : Driver_abortFlash
    use ut_sortInterface, ONLY : ut_sortOnProcs

    implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
    
    integer,intent(IN) :: localCount
    integer,intent(INOUT) ::  localCount_smooth
    integer,intent(IN) :: propCount,propCount_smooth
    integer,intent(IN) :: maxCount_smooth
    real, dimension(propCount_smooth, maxCount_smooth),intent(INOUT) :: dataBuf  
    real,dimension(propCount,localCount),intent(INOUT) :: particles 
    
    integer, dimension(gr_meshNumProcs) :: perProc, toProcs,&
         fromProcs, maxCount
    
    integer, allocatable, dimension(:,:) :: status
    integer, allocatable, dimension(:) :: req
    
    integer :: blockID

    integer :: m, h, parent_id, parent_proc
    integer :: ierr, bufSize
    integer :: gettingFrom, sendingTo
    integer :: sendCount,recvCount, i,j,k
    
    integer, parameter :: tag_for_data = 33
    

    ! Loop over all particles. If that particle originated on local proc, 
    ! give its acceleration back to the original DM particle
    
    ! if not local, put in buffer for transfer,
    ! then receive transferred particles
    
   
    allocate(destBuf(propCount_smooth,maxCount_smooth))
    allocate(sourceBuf(propCount_smooth,maxCount_smooth))
    
    toProcs(:) = 0
    perProc(:) = 0
    sendingTo  = 0
    
    h = 0

    dataBuf(blk_coarse,localCount_smooth+1:maxCount_smooth) = NONEXISTENT

    do m = 1, localCount_smooth
       
       parent_id = int(dataBuf(pid,m))
       parent_proc = int(dataBuf(proc_og,m))
       
       if(parent_id .le. 0) then
          print*, "TransferSmoothedAcceleration, problem with parent_id", dataBuf(:,m)
          call Driver_abortFlash("TransferSmoothedAcceleration, problem with parent_id")
       endif
       if(parent_proc .lt. 0) then
          print*, "TransferSmoothedAcceleration, problem with parent_proc", dataBuf(:,m)
          call Driver_abortFlash("TransferSmoothedAcceleration, problem with parent_proc")
       endif

       if(parent_proc .eq. gr_meshMe) then
          ! local:
          particles(ACCX_PART_PROP,parent_id) = dataBuf(accx,m)
          particles(ACCY_PART_PROP,parent_id) = dataBuf(accy,m)
          particles(ACCZ_PART_PROP,parent_id) = dataBuf(accz,m)
       else
          ! non-local:
          h=h+1
          destBuf(:,h) = dataBuf(:,m)
       endif
       
    end do
    
    call ut_sortOnProcs(h,propCount_smooth,proc_og,gr_meshNumProcs,destBuf,sourceBuf,&
         perProc,toProcs,sendingTo)
    
    call MPI_ALLREDUCE(toProcs,fromProcs,gr_meshNumProcs,FLASH_INTEGER,MPI_SUM,gr_meshComm,ierr)
    ! fromProcs(procNo) now represents the number of other processors which are going to send data
    ! to procNo
    
    call MPI_ALLREDUCE(perProc,maxCount,gr_meshNumProcs,FLASH_INTEGER,MPI_MAX,gr_meshComm,ierr)
    ! maxCount(procNo) is the maximum number of messages that procNo should expect to receive
    ! from a given processor

    
    gettingFrom = fromProcs(gr_meshMe+1)
    recvCount = maxCount(gr_meshMe+1)
  
    allocate(status(MPI_STATUS_SIZE,gettingFrom))
    allocate(req(gettingFrom))
    req(:) = 0
    status(:,:) = 0
    
    
    if(gettingFrom .gt. 0) then
       sourceBuf(:,:) = NONEXISTENT
       bufSize = recvCount * propCount_smooth
       j = 1
       do i = 1, gettingFrom
          call MPI_IRECV(sourceBuf(1,j),bufSize,FLASH_REAL,&
               MPI_ANY_SOURCE,tag_for_data,gr_meshComm,&
               req(i),ierr)
          j = j + recvCount
       enddo
    endif
    

    if(sendingTo .gt. 0) then
       j = 0
       k = 1
       do i = 1, sendingTo
          do while(perProc(j+1) .eq. 0)
             j = j + 1
          enddo
          sendCount = perProc(j+1)
          bufSize   = sendCount * propCount_smooth
          call MPI_SEND(destBuf(1,k),bufSize,FLASH_REAL,j,tag_for_data,gr_meshComm,ierr)
          j = j + 1
          k = k + sendCount
       enddo
    endif


    if(gettingFrom .gt. 0) then
       call MPI_WAITALL(gettingFrom,req,status,ierr)
       
       ! loop over all the newly arrived smooth particles and give acceleration
       ! to parent particle

       do i = 1, recvCount*gettingFrom
          blockID = int(sourceBuf(blk_coarse,i))
          if(blockID .ne. NONEXISTENT) then
             parent_id = int(sourceBuf(pid,i))
             
             if(parent_id .gt. localCount) then
                print*, "parent_id > localCount - can't happen!", parent_id, localCount, localCount_smooth
                do k=1,propCount_smooth
                   print*, k, sourceBuf(k,i)
                end do
             endif
             
             if(parent_id .eq. 0) then
                print*, "parent_id of 0!", sourceBuf(pid,i)
                print*, "i, recvCount*gettingFrom=", i, recvCount*gettingFrom
                do k=1,propCount_smooth
                   print*, k, sourceBuf(k,i)
                end do
                call Driver_abortFlash("parent ID of 0")
             endif
             
             particles(ACCX_PART_PROP,parent_id) = sourceBuf(accx,i)
             particles(ACCY_PART_PROP,parent_id) = sourceBuf(accy,i)
             particles(ACCZ_PART_PROP,parent_id) = sourceBuf(accz,i)
  
         
          endif
          
       enddo
       
    endif
    

    deallocate(req)
    deallocate(status)
    deallocate(destBuf)
    deallocate(sourceBuf)
  
    ! finally, clear out particles_smooth array and number
    
    dataBuf(:,:) = 0.0
    localCount_smooth = 0
      
    return    
  end subroutine TransferSmoothedAccelerations
  






subroutine smooth_checkMassConservation(blocks)
  
  ! make sure sum(rho_del * volume) = sum(rho_dm * volume)
  
  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
       Grid_getListOfBlocks, Grid_getBlkPhysicalSize,  & 
       Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr, & 
       Grid_getBlkPhysicalSize
  use Driver_data, ONLY : dr_globalMe
  
  use Particles_data, ONLY : particles, pt_numLocal

  implicit none

  include "Flash_mpi.h"

  integer, intent(IN) :: blocks
  
  integer :: blockCount
  integer :: blockList(MAXBLOCKS)   
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real,pointer, dimension(:,:,:,: ) :: solnData
  real, dimension(MDIM)    :: size
  
  integer :: i,j,k, lb, blockID, ierr
  
  real :: mass_del, mass_dm, mass, vol
  real :: mass_dm_global, mass_del_global
  
  real :: part_mass_dm, part_mass_dm_global
  
  logical :: first_call = .true.
  
  if(first_call) then
     if(dr_globalMe .eq. MASTER_PE) &
     write(900,*), "Checking conservation of mass between smoothed and unsmoothed DM"
  end if
  
  call Grid_getListOfBlocks(blocks,blockList,blockCount)

  mass_dm = 0.0
  mass_del = 0.0
  
  do lb = 1, blockCount
     
     blockID = blockList(lb)
     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
     call Grid_getBlkPtr(blockID,solnData)
     call Grid_getBlkPhysicalSize(blockID,size)
     
     vol = size(1)*size(2)*size(3) / real(NXB*NYB*NZB)

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS) 
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS) 
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS) 
            
              mass = vol * solnData(PDEN_VAR,i,j,k)
              mass_dm = mass_dm + mass
              
              mass = vol * solnData(PDE_VAR,i,j,k)
              mass_del = mass_del + mass
              
           end do
        end do
     end do
     
     call Grid_releaseBlkPtr(blockID, solnData)
     
  end do     ! blocks
  
  call MPI_ALLREDUCE(mass_dm,mass_dm_global,1,FLASH_REAL,MPI_SUM,MPI_COMM_WORLD,ierr) 
  call MPI_ALLREDUCE(mass_del,mass_del_global,1,FLASH_REAL,MPI_SUM,MPI_COMM_WORLD,ierr) 
  
  ! mass of all particles
  part_mass_dm = 0.0
  if(pt_numLocal .gt. 0) then
     do i = 1, pt_numLocal
       if (particles(TYPE_PART_PROP,i).eq.real(DARK_PART_TYPE)) &
        part_mass_dm = part_mass_dm + particles(MASS_PART_PROP,i) 
     end do
  endif
  
  call MPI_ALLREDUCE(part_mass_dm,part_mass_dm_global,1,FLASH_REAL,MPI_SUM,MPI_COMM_WORLD,ierr) 

  if(dr_globalMe .eq. MASTER_PE) then
     write(900,*), "particle mass dm=", part_mass_dm_global
     write(900,*), "mass_dm:", mass_dm_global, part_mass_dm_global/mass_dm_global
     write(900,*), "mass_del:", mass_del_global, part_mass_dm_global/mass_del_global  
     write(900,*), "*************"
  end if

  return
end subroutine smooth_checkMassConservation
  
  






subroutine smooth_findTotalMass(var,totMass)
  ! find total mass on grid of variable var 
  ! (assumes var is a density)
  
  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
       Grid_getListOfBlocks, Grid_getBlkPhysicalSize,  & 
       Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr, & 
       Grid_getBlkPhysicalSize
  
  implicit none

  include "Flash_mpi.h"

  integer, intent(IN) :: var
  real, intent(OUT) :: totMass
  
  integer :: blockCount
  integer :: blockList(MAXBLOCKS)   
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real,pointer, dimension(:,:,:,: ) :: solnData
  real, dimension(MDIM)    :: size
  
  integer :: i,j,k, lb, blockID, ierr
  
  real :: mass, vol
  
  call Grid_getListOfBlocks(ALL_BLKS,blockList,blockCount)

  mass = 0.0
  
  do lb = 1, blockCount
     
     blockID = blockList(lb)
     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
     call Grid_getBlkPtr(blockID,solnData)
     call Grid_getBlkPhysicalSize(blockID,size)
     
     vol = size(1)*size(2)*size(3) / real(NXB*NYB*NZB)

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS) 
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS) 
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS) 
            
              mass = vol * solnData(var,i,j,k) + mass
              
           end do
        end do
     end do
     
     call Grid_releaseBlkPtr(blockID, solnData)
     
  end do     ! blocks
  
  call MPI_ALLREDUCE(mass,totMass,1,FLASH_REAL,MPI_SUM,MPI_COMM_WORLD,ierr) 


  return
  
end subroutine smooth_findTotalMass
  
  

  



end module Smooth_particles
