subroutine chebfit1( &
     x, &
     y, &
     n, &
     dim, &
     degr, &
     c )
  implicit none
  integer,                      intent(in)  :: n, dim
  real(kind=MATHpr),            intent(in)  :: x(n)
  real(kind=MATHpr),            intent(in)  :: y(n,dim)
  integer,                      intent(in)  :: degr
  type(type_chebyshev_series1), intent(out) :: c
  real(kind=MATHpr)                         :: F(n,degr+1), twox(n), Ft(degr+1,n)
  integer                                   :: i, j

  F(:,:) = 0._MATHpr
  F(:,1) = 1._MATHpr
  F(:,2) = x

  if ( degr > 1 ) twox = 2._MATHpr * x

  do i = 3,degr+1
     F(:,i) = twox * F(:,i-1) - F(:,i-2)
  end do
  Ft = transpose(F)

  call reset_chebyshev_series1( c, degr, dim )

  do j = 1,dim
     call solve_NxN( c%coef(1:degr+1,j), matmul( Ft, F ), matmul( Ft, y(:,j) ) )
  end do

end subroutine chebfit1
