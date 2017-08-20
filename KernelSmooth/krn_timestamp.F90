
! print a log message with current time stamp

subroutine krn_timestamp(msg)

  use Kernel_data, ONLY : krn_mype

  implicit none
  character (len=*), intent(in) :: msg
  character (len=40) :: time_string

#include "constants.h"

  if (krn_mype.eq.MASTER_PE) then
    call current_date_time(time_string)
    print '(4A)', 'KERNEL: ', trim(time_string), ' | ', msg
  end if

end subroutine
  
