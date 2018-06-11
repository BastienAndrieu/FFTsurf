subroutine chebval2_n( &
     f, &
     c, &
     xy, &
     n )
  implicit none
  integer,                      intent(in)  :: n
  real(kind=MATHpr),            intent(in)  :: xy(2,n)
  type(type_chebyshev_series2), intent(in)  :: c
  real(kind=MATHpr),            intent(out) :: f(c%dim,n)
  integer                                   :: i

  do i = 1,n
     call chebval2( &
          f(:,i), &
          c, &
          xy(:,i) )
  end do

end subroutine chebval2_n
