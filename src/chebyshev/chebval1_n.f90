subroutine chebval1_n( &
     f, &
     c, &
     x, &
     n )
  implicit none
  integer,                      intent(in)  :: n
  real(kind=MATHpr),            intent(in)  :: x(n)
  type(type_chebyshev_series1), intent(in)  :: c
  real(kind=MATHpr),            intent(out) :: f(n,c%dim)

  call clenshaw( &
       f, &
       c%coef, &
       x, &
       n, &
       c%degr, &
       c%dim )

end subroutine chebval1_n
