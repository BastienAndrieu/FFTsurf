function ith_cgl_point( i, n ) result( x )
  ! Returns the i-th of the n+1 Chebyshev-Gauss-Lobatto points
  ! xi = cos( i\pi/n ) with 0 <= i <= n
  ! /!\ CGL points are ordered from +1 to -1 (x0 = +1, ..., xn = -1)
  implicit none
  integer, intent(in) :: i, n
  real(kind=MATHpr)   :: x

  x = cos( real( i, kind=MATHpr ) * MATHpi / real( n, kind=MATHpr ) )

end function ith_cgl_point
