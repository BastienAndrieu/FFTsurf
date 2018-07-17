subroutine classify_endpoint( &
     uv, &
     tng, &
     uvbox, &
     stat )
  use mod_math
  use mod_tolerances
  !stat = -1 : exiting
  !        0 : interior
  !        1 : entering
  !        2 : ambiguous
  implicit none
  real(kind=fp), intent(in)  :: uv(2)
  real(kind=fp), intent(in)  :: tng(2) ! tangential direction
  real(kind=fp), intent(in)  :: uvbox(4)
  integer,       intent(out) :: stat
  real(kind=fp)              :: uvloc
  integer                    :: statvar(2)
  integer                    :: ivar

  do ivar = 1,2
     uvloc = ab2n1p1( uv(ivar), uvbox(2*ivar-1), uvbox(2*ivar) )
     if ( is_in_open_interval(uvloc, -1._fp, 1._fp, tolerance=EPSuv) ) then
        statvar(ivar) = 0
     else
        if ( abs(tng(ivar)) < EPSuv ) then
           statvar(ivar) = 2
        elseif ( tng(ivar)*uvloc > 0._fp ) then
           statvar(ivar) = -1
        else
           statvar(ivar) =  1
        end if
     end if
  end do
  
  if ( any(statvar == -1) ) then
     stat = -1 ! exiting
  else
     if ( any(statvar == 1) ) then
        stat = 1 ! entering
     else
        stat = sum(statvar) ! interior or ambiguous
     end if
  end if

end subroutine classify_endpoint
