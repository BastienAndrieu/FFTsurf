subroutine fcht1( &
     f, &
     c, &
     m, &
     dim, &
     epsilon )
  ! Performs Fast Cosine Transform (using dfftpack routines)
  implicit none
  integer,                      intent(in)           :: m, dim
  real(kind=MATHpr),            intent(in)           :: f(m,dim)
  type(type_chebyshev_series1), intent(out)          :: c
  real(kind=MATHpr),            intent(in), optional :: epsilon
  real(kind=MATHpr), allocatable                     :: w(:)
  integer                                            :: j, p, degrprev
  save w, degrprev
  data degrprev /-1/

  call reset_chebyshev_series1( c, m-1, dim )

  c%coef = f / real( m-1, kind=MATHpr )

  p = sizeminw( m )
  if ( allocated(w) ) then
     if ( size(w) < p ) deallocate( w )
  end if
  if ( .not.allocated( w ) ) allocate( w(p) )

  if ( degrprev /= m-1 ) call dcosti( m, w(1:p) )
  degrprev = m-1

  do j = 1,dim
     call dcost( m, c%coef(:,j), w(1:p) )
  end do

  c%coef([1,m],:) = 0.5_MATHpr * c%coef([1,m],:)

  if ( present(epsilon) ) then
     where ( abs(c%coef) < epsilon ) c%coef = 0._MATHpr
  end if

end subroutine fcht1
