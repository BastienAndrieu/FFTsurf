subroutine fcht2( &
     f, &
     c, &
     m, &
     n, &
     dim, &
     epsilon )
  implicit none
  integer,                      intent(in)           :: m, n, dim
  real(kind=MATHpr),            intent(in)           :: f(m, n, dim)
  type(type_chebyshev_series2), intent(out)          :: c
  real(kind=MATHpr),            intent(in), optional :: epsilon
  real(kind=MATHpr)                                  :: w(3*max(m,n)+15)
  integer                                            :: degr(2), i, j, k

  degr(1) = m - 1
  degr(2) = n - 1

  call reset_chebyshev_series2( c, degr, dim )

  c%coef = f / real( degr(1), kind=MATHpr )
  call dcosti( m, w )
  do k = 1,dim
     do j = 1,n
        call dcost( m, c%coef(:,j,k), w )
     end do
  end do
  c%coef([1,m],:,:) = 0.5_MATHpr * c%coef([1,m],:,:)

  c%coef = c%coef / real( degr(2), kind=MATHpr )
  if ( degr(1) /= degr(2) ) call dcosti( n, w )
  do k = 1,dim
     do i = 1,m
        call dcost( n, c%coef(i,:,k), w )
     end do
  end do
  c%coef(:,[1,n],:) = 0.5_MATHpr * c%coef(:,[1,n],:)

  if ( present(epsilon) ) then
     where ( abs(c%coef) < epsilon ) c%coef = 0._MATHpr
  end if

end subroutine fcht2
