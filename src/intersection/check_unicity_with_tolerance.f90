subroutine check_unicity_with_tolerance( &
     unew, &
     xnew, &
     tolunew, &
     dim, &
     u, &
     x, &
     tolu, &
     n, &
     id )
  use mod_math
  use mod_tolerances
  implicit none
  integer,       intent(in)  :: dim
  real(kind=fp), intent(in)  :: unew(dim)
  real(kind=fp), intent(in)  :: xnew(3)
  real(kind=fp), intent(in)  :: tolunew
  integer,       intent(in)  :: n
  real(kind=fp), intent(in)  :: u(dim,n)
  real(kind=fp), intent(in)  :: x(3,n)
  real(kind=fp), intent(in)  :: tolu(n)
  integer,       intent(out) :: id

  if ( n < 1 ) then
     id = 1
  else
     do id = 1,n
        if ( maxval(abs(unew - u(:,id))) < tolunew + tolu(id) .or. & 
           !sum((unew - u(:,id))**2) < (tolunew + tolu(id))**2 .or. &
             sum((xnew - x(:,id))**2) < EPSxyzsqr ) then
           exit
        end if
     end do
  end if
 
end subroutine check_unicity_with_tolerance
