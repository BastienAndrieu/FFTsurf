subroutine bernbivar2univar( &
    b2, &
    m, &
    n, &
    p, &
    b1, &
    ivar, &
    ival )
  use mod_constants
  implicit none
  integer,       intent(in)  :: m, n, p
  real(kind=fp), intent(in)  :: b2(m,n,p)
  integer,       intent(in)  :: ivar, ival
  real(kind=fp), intent(out) :: b1(size(b2,1+mod(ivar,2)),p)
  integer                    :: j

  j = 1 + (ival - 1)*size(b2,ivar)
  if ( ivar == 1 ) then
     b1 = b2(j,:,:)
  else
     b1 = b2(:,j,:)
  end if

end subroutine bernbivar2univar
