
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Check to see if it's time to zoom, or if zooming has completed and a SN needs
!! to be inserted

subroutine sim_checkStatus()

  use Simulation_data, ONLY : sim_num_stars, sim_zoom_status, sim_nova_time, sim_nova_started,&
                              sim_zoom_time, sim_ismaster
  use sim_interface, ONLY : sim_initNova, sim_setupNova
  use Driver_interface, ONLY: Driver_getSimTime
  
  implicit none
  
#include "Flash_mpi.h"
  
#include "zoom.h"
 
  integer, parameter :: ZERO = 0

  logical :: eventStarted
 
  integer :: n, ierr
  
  real :: time
  
  logical, dimension(sim_num_stars) :: local_ready, all_ready
   
  eventStarted = .FALSE.

  call Driver_getSimTime(time)
  
  ! see if the next SN is ready to start refining
  do n = 1, sim_num_stars
    if ((time.gt.sim_zoom_time(n)).and.(sim_zoom_status(n).eq.ZOOM_NOT_STARTED)) then
      sim_zoom_status(n) = ZOOM_REFINING
      
      ! pick a set of coordinates and calculate the radius/velocity from the local density
      
      if (sim_ismaster) print *, 'SIM: Starting supernova ', n, ' at time ', time, sim_zoom_time(n)
      
      call sim_setupNova(n)
    end if
  end do
  
  ! see if any of the supernovae are in the state ZOOM_LOCAL_READY
  local_ready = sim_zoom_status.eq.ZOOM_LOCAL_READY

  ! reduce to see if they are ready on all processors
  call MPI_ALLREDUCE(local_ready, all_ready, sim_num_stars, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
  
  ! mark any that are ready on all processors as being done and initialize the nova
  do n = 1,sim_num_stars
    ! if they're not all done, but they're ready, mark them
    if (sim_zoom_status(n).ne.ZOOM_ALL_DONE .and. all_ready(n)) then
      eventStarted = .TRUE.
      sim_zoom_status(n) = ZOOM_ALL_DONE
    end if
    
    ! if they are all done, but we haven't initialized the nova, go ahead and do that
    if (sim_zoom_status(n).eq.ZOOM_ALL_DONE .and. sim_nova_started(n).eq.ZERO) then
      eventStarted = .TRUE.
      call sim_initNova(n)
    end if
  end do

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

end subroutine sim_checkStatus
