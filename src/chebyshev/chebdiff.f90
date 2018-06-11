subroutine chebdiff( &
     c, &
     d, &
     degr, &
     dim )
  implicit none
  integer,           intent(in)  :: degr, dim
  real(kind=MATHpr), intent(in)  :: c(degr+1,dim)
  real(kind=MATHpr), intent(out) :: d(max(degr,1),dim)
  integer                        :: i

  d(:,:) = 0._MATHpr
  if ( degr < 1 ) return

  d(degr,:) = real( 2*degr, kind=MATHpr ) * c(degr+1,:)

  if ( degr > 1 ) then
     d(degr-1,:) = real( 2*(degr-1), kind=MATHpr ) * c(degr,:)

     if ( degr > 2 ) then
        do i = degr-2,1,-1
           d(i,:) = real( 2*i, kind=MATHpr ) * c(i+1,:) + d(i+2,:)
        end do
     end if

  end if
  d(1,:) = 0.5_MATHpr * d(1,:)

end subroutine chebdiff
