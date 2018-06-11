subroutine berndiff1( &
     b, &
     d )
  implicit none
  type(type_bernstein_series1), intent(in)  :: b
  type(type_bernstein_series1), intent(out) :: d

  call reset_bernstein_series1( d, max( b%degr-1, 0 ), b%dim )
  call berndiff( b%coef(1:b%degr+1,:), d%coef(1:max(b%degr,1),:), b%degr, b%dim )

end subroutine berndiff1
