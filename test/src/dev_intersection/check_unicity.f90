subroutine check_unicity( &
     vec, &
     dim, &
     array, &
     n, &
     tol, &
     id )
  use mod_math
  implicit none
  integer,       intent(in)  :: dim
  real(kind=fp), intent(in)  :: vec(dim)
  integer,       intent(in)  :: n
  real(kind=fp), intent(in)  :: array(dim,n)
  real(kind=fp), intent(in)  :: tol
  integer,       intent(out) :: id
  real(kind=fp)              :: tolsqr

  if ( n < 1 ) then
     id = 1
  else
     tolsqr = tol**2
     do id = 1,n
        if ( sum((vec(1:dim) - array(1:dim,id))**2) < tolsqr ) exit
     end do
  end if
  
end subroutine check_unicity
