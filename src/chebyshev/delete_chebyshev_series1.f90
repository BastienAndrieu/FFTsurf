subroutine delete_chebyshev_series1( &
     c )
  implicit none
  type(type_chebyshev_series1), intent(inout) :: c

  if (allocated(c%coef)) deallocate(c%coef)
  c%dim = 0
  c%degr = 0

end subroutine delete_chebyshev_series1
