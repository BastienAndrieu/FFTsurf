subroutine subdiv2( &
     b, &
     uv, &
     bsw, &
     bse, &
     bnw, &
     bne )
  use mod_math
  implicit none
  type(type_bernstein_series2), intent(in)              :: b
  real(kind=MATHpr),            intent(in)              :: uv(2)
  type(type_bernstein_series2), intent(out), optional   :: bsw
  type(type_bernstein_series2), intent(out), optional   :: bse
  type(type_bernstein_series2), intent(out), optional   :: bnw
  type(type_bernstein_series2), intent(out), optional   :: bne
  real(kind=MATHpr), dimension(b%degr(1)+1,b%degr(2)+1) :: bw, be
  real(kind=MATHpr), dimension(b%degr(2)+1,b%degr(1)+1) :: bsT, bnT
  integer                                               :: m, n, k

  m = b%degr(1)+1
  n = b%degr(2)+1

  if ( present(bsw) ) call reset_bernstein_series2( bsw, b%degr, b%dim )
  if ( present(bse) ) call reset_bernstein_series2( bse, b%degr, b%dim )
  if ( present(bnw) ) call reset_bernstein_series2( bnw, b%degr, b%dim )
  if ( present(bne) ) call reset_bernstein_series2( bne, b%degr, b%dim )

  do k = 1,b%dim
     call de_casteljau( &
          b%coef(1:m,1:n,k), &
          uv(1), &
          b%degr(1), &
          n, &
          bl=bw, &
          br=be )

     call de_casteljau( &
          transpose(bw), &
          uv(2), &
          b%degr(2), &
          m, &
          bl=bsT, &
          br=bnT )
     if ( present(bsw) ) bsw%coef(1:m,1:n,k) = transpose(bsT)
     if ( present(bnw) ) bnw%coef(1:m,1:n,k) = transpose(bnT)
     
     call de_casteljau( &
          transpose(be), &
          uv(2), &
          b%degr(2), &
          m, &
          bl=bsT, &
          br=bnT )
     if ( present(bse) ) bse%coef(1:m,1:n,k) = transpose(bsT)
     if ( present(bne) ) bne%coef(1:m,1:n,k) = transpose(bnT)     

  end do

end subroutine subdiv2
