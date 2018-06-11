subroutine ifcht1( &
     c, &
     f )
  implicit none
  type(type_chebyshev_series1), intent(in)  :: c
  real(kind=MATHpr),            intent(out) :: f(c%degr+1,c%dim)
  real(kind=MATHpr), allocatable            :: w(:)    
  integer                                   :: j, degrprev, p
  save w, degrprev
  data degrprev /-1/

  f = c%coef(1:c%degr+1,:)
  f(2:c%degr,:) = 0.5_MATHpr * f(2:c%degr,:)

  p = sizeminw( c%degr+1 )
  if ( allocated(w) ) then
     if ( size(w) < p ) deallocate( w )
  end if
  if ( .not.allocated( w ) ) allocate( w(p) )

  if ( degrprev /= c%degr ) call dcosti( c%degr+1, w )
  degrprev = c%degr

  do j = 1,c%dim
     call dcost( c%degr+1, f(:,j), w )
  end do

end subroutine ifcht1
