subroutine subdiv2_along_u( &
     b, &
     u, &
     bw, &
     be )
  implicit none
  type(type_bernstein_series2), intent(in)  :: b
  real(kind=MATHpr),            intent(in)  :: u
  type(type_bernstein_series2), intent(out) :: bw
  type(type_bernstein_series2), intent(out) :: be
  integer                                   :: m, n, k

  m = b%degr(1)+1
  n = b%degr(2)+1
  
  call reset_bernstein_series2( bw, b%degr, b%dim )
  call reset_bernstein_series2( be, b%degr, b%dim )

  do k = 1,b%dim
     call de_casteljau( &
          b%coef(1:m,1:n,k), &
          u, &
          b%degr(1), &
          n, &
          bl=bw%coef(1:m,1:n,k), &
          br=be%coef(1:m,1:n,k) )
  end do

end subroutine subdiv2_along_u
