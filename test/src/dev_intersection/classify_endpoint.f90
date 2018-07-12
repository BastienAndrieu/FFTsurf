subroutine classify_endpoint( &
     uv, &
     tng, &
     uvbox, &
     stat )
  use mod_math
  !stat = -1 : exiting
  !        0 : interior
  !        1 : entering
  implicit none
  real(kind=fp), intent(in)  :: uv(2)
  real(kind=fp), intent(in)  :: tng(2) ! tangential direction
  real(kind=fp), intent(in)  :: uvbox(4)
  integer,       intent(out) :: stat
  real(kind=fp)              :: uvloc
  integer                    :: ivar

  stat = 0
  do ivar = 1,2
     uvloc = ab2n1p1( uv(ivar), uvbox(2*ivar-1), uvbox(2*ivar) )
     if ( abs(uvloc) > 1._fp - EPSmath ) then
        if ( tng(ivar)*uvloc > 0._fp ) then
           stat = -1 ! exiting
           return
        else
           stat = 1 ! entering
        end if
     end if
  end do

end subroutine classify_endpoint
