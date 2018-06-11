subroutine ifcht2_padded( &
     c, &
     f, &
     m, &
     n )
  implicit none
  type(type_chebyshev_series2), intent(in)  :: c
  integer,                      intent(in)  :: m, n
  real(kind=MATHpr),            intent(out) :: f(m,n,c%dim)
  real(kind=MATHpr)                         :: w(3*max(m,n)+15)
  integer                                   :: i, j, k, r, s

  r = min( m, c%degr(1)+1 )
  s = min( n, c%degr(2)+1 )

  f(1:r,1:s,:) = c%coef(1:r,1:s,:)
  if ( r < m ) f(r+1:m,:,:) = 0._MATHpr
  if ( s < n ) f(1:r,s+1:n,:) = 0._MATHpr

  f(2:m-1,:,:) = 0.5_MATHpr * f(2:m-1,:,:)
  call dcosti( m, w )
  do k = 1,c%dim
     do j = 1,n
        call dcost( m, f(:,j,k), w )
     end do
  end do

  f(:,2:n-1,:) = 0.5_MATHpr * f(:,2:n-1,:)
  if ( m /= n ) call dcosti( n, w )
  do k = 1,c%dim
     do i = 1,m
        call dcost( n, f(i,:,k), w )
     end do
  end do

end subroutine ifcht2_padded
