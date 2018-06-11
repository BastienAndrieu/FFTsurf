subroutine chebdiff1( &
     c, &
     d )
  implicit none
  type(type_chebyshev_series1), intent(in)  :: c
  type(type_chebyshev_series1), intent(out) :: d

  call reset_chebyshev_series1( d, max( c%degr-1, 0 ), c%dim )
  call chebdiff( c%coef(1:c%degr+1,:), d%coef(1:max(c%degr,1),:), c%degr, c%dim )

end subroutine chebdiff1
