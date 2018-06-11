subroutine chebval2_corners( &
     c, &
     f )
  implicit none
  type(type_chebyshev_series2), intent(in)              :: c
  real(kind=MATHpr),            intent(out)             :: f(c%dim,4)
  integer                                               :: i = 0, m, n
  real(kind=MATHpr), dimension(c%degr(1)+1,c%degr(2)+1) :: om, on

  om = spread( &
       source = real( (-1)**[ ( i, i=0,c%degr(1) ) ], kind=MATHpr ), &
       dim = 2, &
       ncopies = c%degr(2)+1 )
  on = spread( &
       source = real( (-1)**[ ( i, i=0,c%degr(2) ) ], kind=MATHpr ), &
       dim = 1, &
       ncopies = c%degr(1)+1 )

  m = c%degr(1)+1
  n = c%degr(2)+1

  f(:,:) = 0._MATHpr
  do i = 1,c%dim
     f(i,1) = sum( sum( c%coef(1:m,1:n,i) * om * on, 2 ), 1 )
     f(i,2) = sum( sum( c%coef(1:m,1:n,i) * on, 2 ), 1 )
     f(i,3) = sum( sum( c%coef(1:m,1:n,i) * om, 2 ), 1 )
     f(i,4) = sum( sum( c%coef(1:m,1:n,i), 2 ), 1 )
  end do

end subroutine chebval2_corners
