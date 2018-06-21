subroutine rearrange_for_separability_test( &
     bcp, &
     nbcp, &
     pivot, &
     origin, &
     sep, &
     nsep )
  use mod_math
  use mod_tolerances
  implicit none
  integer,       intent(in)  :: nbcp
  real(kind=fp), intent(in)  :: bcp(nbcp,3)
  real(kind=fp), intent(in)  :: pivot(3)
  real(kind=fp), intent(in)  :: origin(3)
  real(kind=fp), intent(out) :: sep(nbcp,3)
  integer,       intent(out) :: nsep
  integer                    :: i

  nsep = 0
  do i = 1,nbcp
     if ( sum( (bcp(i,:) - pivot)**2 ) > EPSxyzsqr ) then
        nsep = nsep + 1
        sep(nsep,:) = bcp(i,:) - origin
     end if
  end do

end subroutine rearrange_for_separability_test
