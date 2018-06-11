subroutine berndiff( &
     b, &
     d, &
     degr, &
     dim )
  implicit none
  integer,           intent(in)  :: degr, dim
  real(kind=MATHpr), intent(in)  :: b(degr+1,dim)
  real(kind=MATHpr), intent(out) :: d(max(degr,1),dim)
  
  if ( degr < 1 ) then
     d(:,:) = 0._MATHpr
  else
     d = real( degr, kind=MATHpr ) * ( b(2:degr+1,:) - b(1:degr,:) )
  end if

end subroutine berndiff
