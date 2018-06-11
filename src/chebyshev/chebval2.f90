subroutine chebval2( &
     f, &
     c, &
     xy )
  implicit none
  real(kind=MATHpr),            intent(in)  :: xy(2)
  type(type_chebyshev_series2), intent(in)  :: c
  real(kind=MATHpr),            intent(out) :: f(c%dim)
  real(kind=MATHpr)                         :: b(c%degr(2)+1,c%dim)
  integer                                   :: k

  do k = 1,c%dim
     call clenshaw( b(:,k), c%coef(:,:,k), [xy(1)], 1, c%degr(1), c%degr(2)+1 )
  end do
  call clenshaw( f, b, [xy(2)], 1, c%degr(2), c%dim )

end subroutine chebval2
