subroutine economization1( &
     c, &
     tol )
  implicit none
  type(type_chebyshev_series1), intent(inout) :: c
  real(kind=MATHpr),            intent(in)    :: tol
  real(kind=MATHpr)                           :: tolsqr
  integer                                     :: i, M

  tolsqr = tol**2

  M = c%degr+1
  do i = M,1,-1
     if ( sum( c%coef(i:M,:)**2 ) > tolsqr ) exit
  end do
  c%degr = i-1

end subroutine economization1
