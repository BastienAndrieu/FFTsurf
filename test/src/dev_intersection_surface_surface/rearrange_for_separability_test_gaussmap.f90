subroutine rearrange_for_separability_test_gaussmap( &
     bcp, &
     nbcp, &
     norcol, &
     sep, &
     nsep )
  use mod_math
  use mod_tolerances
  implicit none
  integer,       intent(in)  :: nbcp
  real(kind=fp), intent(in)  :: bcp(nbcp,3)
  real(kind=fp), intent(in)  :: norcol(3)
  real(kind=fp), intent(out) :: sep(nbcp,3)
  integer,       intent(out) :: nsep
  real(kind=fp)              :: vec(3)
  integer                    :: i

  nsep = 0
  do i = 1,nbcp
     vec = bcp(i,:) / norm2( bcp(i,:) )
     if ( abs(dot_product(vec, norcol)) >= 1._fp - EPSxyz ) cycle
     nsep = nsep + 1
     sep(nsep,:) = vec
  end do

end subroutine rearrange_for_separability_test_gaussmap
