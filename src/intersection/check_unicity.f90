subroutine check_unicity( &
     vec, &
     dim, &
     array, &
     n, &
     tol, &
     id )
  use mod_math
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  integer,       intent(in)  :: dim
  real(kind=fp), intent(in)  :: vec(dim)
  integer,       intent(in)  :: n
  real(kind=fp), intent(in)  :: array(dim,n)
  real(kind=fp), intent(in)  :: tol
  integer,       intent(out) :: id
  real(kind=fp)              :: tolsqr

  IF ( DEBUG ) PRINT *,'CHECK UNICITY OF',VEC(1:DIM)
  
  if ( n < 1 ) then
     id = 1
  else
     tolsqr = tol**2
     do id = 1,n
        IF ( DEBUG ) PRINT *,array(1:dim,id), NORM2(vec(1:dim) - array(1:dim,id))
        if ( sum((vec(1:dim) - array(1:dim,id))**2) < tolsqr ) then
           IF ( DEBUG ) PRINT '(A17,1x,I0,1x,A1)','DUPLICATE POINT (',ID,')'
           exit
        end if
     end do
  end if
  IF ( DEBUG ) PRINT '(A4,1X,I0)','ID =',ID
  
end subroutine check_unicity
