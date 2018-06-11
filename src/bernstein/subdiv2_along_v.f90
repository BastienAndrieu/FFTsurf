subroutine subdiv2_along_v( &
     b, &
     v, &
     bs, &
     bn )
  implicit none
  type(type_bernstein_series2), intent(in)              :: b
  real(kind=MATHpr),            intent(in)              :: v
  type(type_bernstein_series2), intent(out)             :: bs
  type(type_bernstein_series2), intent(out)             :: bn
  real(kind=MATHpr), dimension(b%degr(2)+1,b%degr(1)+1) :: bsT, bnT
  integer                                               :: m, n, k

  m = b%degr(1)+1
  n = b%degr(2)+1
  
  call reset_bernstein_series2( bs, b%degr, b%dim )
  call reset_bernstein_series2( bn, b%degr, b%dim )

  do k = 1,b%dim
     call de_casteljau( &
          transpose( b%coef(1:m,1:n,k) ), &
          v, &
          b%degr(2), &
          m, &
          bl=bsT, &
          br=bnT )
     bs%coef(1:m,1:n,k) = transpose(bsT)
     bn%coef(1:m,1:n,k) = transpose(bnT)
  end do

end subroutine subdiv2_along_v
