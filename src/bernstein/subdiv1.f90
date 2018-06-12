subroutine subdiv1( &
     b, &
     t, &
     bl, &
     br )
  use mod_math
  implicit none
  type(type_bernstein_series1), intent(in)  :: b
  real(kind=MATHpr),            intent(in)  :: t
  type(type_bernstein_series1), intent(out) :: bl
  type(type_bernstein_series1), intent(out) :: br

  call reset_bernstein_series1( bl, b%degr, b%dim )
  call reset_bernstein_series1( br, b%degr, b%dim )

  call de_casteljau( &
       b%coef(1:b%degr+1,:), &
       t, &
       b%degr, &
       b%dim, &
       bl=bl%coef(1:b%degr+1,:), &
       br=br%coef(1:b%degr+1,:) )

end subroutine subdiv1
