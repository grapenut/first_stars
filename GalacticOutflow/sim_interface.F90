
Module sim_interface

implicit none

  interface
    subroutine sim_markInitial()
    end subroutine sim_markInitial
  end interface

  interface
    subroutine sim_markOverdensity()
    end subroutine sim_markOverdensity
  end interface

  interface
    subroutine sim_markZoom()
    end subroutine sim_markZoom
  end interface

  interface
    subroutine sim_setupNova(n)
      integer, intent(IN) :: n
    end subroutine sim_setupNova
  end interface

  interface
    subroutine sim_initNova(n)
      integer, intent(IN) :: n
    end subroutine sim_initNova
  end interface

  interface
    subroutine sim_checkStatus()
    end subroutine sim_checkStatus
  end interface


end Module sim_interface

