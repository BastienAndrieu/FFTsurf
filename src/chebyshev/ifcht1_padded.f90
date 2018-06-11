subroutine ifcht1_padded( &
     c, &
     f, &
     m )
  implicit none
  type(type_chebyshev_series1), intent(in)  :: c
  integer,                      intent(in)  :: m
  real(kind=MATHpr),            intent(out) :: f(m,c%dim)
  real(kind=MATHpr), allocatable            :: w(:)    
  integer                                   :: j, degrprev, r, p
  save w, degrprev
  data degrprev /0/

  r = min( m, c%degr+1 )
  f(1:r,:) = c%coef(1:r,:)
  if ( r < m ) f(r+1:m,:) = 0._MATHpr

  p = sizeminw( m )
  if ( allocated(w) ) then
     if ( size(w) < p ) deallocate( w )
  end if
  if ( .not.allocated( w ) ) allocate( w(p) )

  if ( degrprev /= m-1 ) call dcosti( m, w )
  degrprev = m-1

  f(2:m-1,:) = 0.5_MATHpr * f(2:m-1,:)
  do j = 1,c%dim
     call dcost( m, f(:,j), w )
  end do

end subroutine ifcht1_padded
