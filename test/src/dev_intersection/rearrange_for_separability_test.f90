subroutine rearrange_for_separability_test( &
     bcp, &
     nbcp, &
     vec, &
     sep, &
     nsep, &
     gaussmap )
  use mod_math
  use mod_tolerances
  implicit none
  integer,       intent(in)  :: nbcp
  real(kind=fp), intent(in)  :: bcp(nbcp,3)
  real(kind=fp), intent(in)  :: vec(3)
  logical,       intent(in)  :: gaussmap
  real(kind=fp), intent(out) :: sep(nbcp,3)
  integer,       intent(out) :: nsep
  real(kind=fp)              :: tmp(3), tmpsqr
  integer                    :: i

  if ( gaussmap ) tmpsqr = sum( vec**2 )
  nsep = 0
  do i = 1,nbcp

     if ( gaussmap ) then
        tmp = bcp(i,:) / norm2(bcp(i,:))
        if ( tmpsqr > EPSxyzsqr .and. abs(dot_product(tmp, vec)) >= 1._fp - EPSxyz ) cycle
     else
        tmp = bcp(i,:) - vec
        tmpsqr = sum( tmp**2 )
        if ( tmpsqr <= EPSxyzsqr ) cycle
        tmp = tmp / sqrt(tmpsqr)
     end if

     nsep = nsep + 1
     sep(nsep,:) = tmp

  end do

end subroutine rearrange_for_separability_test
