subroutine chebval1( &
     f, &
     c, &
     x )
  implicit none
  real(kind=MATHpr),            intent(in)  :: x
  type(type_chebyshev_series1), intent(in)  :: c
  real(kind=MATHpr),            intent(out) :: f(c%dim)

  call clenshaw( &
       f, &
       c%coef(1:c%degr+1,:), &
       [x], &
       1, &
       c%degr, &
       c%dim )

end subroutine chebval1
