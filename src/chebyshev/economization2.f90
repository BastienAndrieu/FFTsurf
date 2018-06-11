subroutine economization2( &
     c, &
     tol )
  implicit none
  type(type_chebyshev_series2), intent(inout) :: c
  real(kind=MATHpr),            intent(in)    :: tol
  real(kind=MATHpr)                           :: tolsqr
  real(kind=MATHpr), dimension(c%dim)         :: s, r
  integer                                     :: i, j

  tolsqr = tol**2

  s = sum( sum( abs(c%coef), 2), 1 )

  do i = c%degr(1)+1,1,-1
     r = sum( sum( abs(c%coef(1:i-1,1:c%degr(2)+1,:)), 2), 1 )
     if ( sum( (s - r)**2 ) > tolsqr ) exit
  end do
  c%degr(1) = i-1

  do j = c%degr(2)+1,1,-1
     r = sum( sum( abs(c%coef(1:i,1:j-1,:)), 2), 1 )
     if ( sum( (s - r)**2 ) > tolsqr ) exit
  end do
  c%degr(2) = j-1

end subroutine economization2
