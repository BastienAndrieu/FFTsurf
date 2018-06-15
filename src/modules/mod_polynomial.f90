module mod_polynomial

  use mod_constants
  implicit none

  type type_polynomial
     integer                        :: nvar = 0
     integer                        :: base = 0
     integer                        :: degr(2) = -1
     integer                        :: dim = 0
     real(kind=fp), allocatable     :: coef(:,:,:)
  end type type_polynomial

contains

  subroutine free_polynomial( poly )
    implicit none
    type(type_polynomial), intent(inout) :: poly
    if ( allocated(poly%coef) ) deallocate( poly%coef )
  end subroutine free_polynomial


end module mod_polynomial
