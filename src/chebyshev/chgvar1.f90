subroutine chgvar1( &
     c, &
     s, &
     x0, &
     x1 )
  implicit none
  type(type_chebyshev_series1), intent(in)  :: c
  real(kind=MATHpr),            intent(in)  :: x0
  real(kind=MATHpr),            intent(in)  :: x1
  type(type_chebyshev_series1), intent(out) :: s

  call reset_chebyshev_series1( s, c%degr, c%dim )

  call chgvar( &
       c%coef(1:c%degr+1,:), &
       s%coef(1:c%degr+1,:), &
       x0, &
       x1, &
       c%degr, &
       c%dim )

end subroutine chgvar1
