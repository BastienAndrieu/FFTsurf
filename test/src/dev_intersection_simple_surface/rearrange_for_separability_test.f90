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
  real(kind=fp)              :: vec(3), vecsqr
  integer                    :: i

  nsep = 0
  do i = 1,nbcp
     if ( sum( (bcp(i,:) - pivot)**2 ) <= EPSxyzsqr ) cycle
     nsep = nsep + 1
     !sep(nsep,:) = bcp(i,:) - origin
     vec = bcp(i,:) - origin
     vecsqr = sum( vec**2 )
     if ( vecsqr > 3._fp*EPSfpsqr ) vec = vec / sqrt( vecsqr )
     sep(nsep,:) = vec
  end do

end subroutine rearrange_for_separability_test
