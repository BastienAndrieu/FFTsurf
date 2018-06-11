subroutine chgvar2( &
     c, &
     s, &
     xy0, &
     xy1 )
  implicit none
  type(type_chebyshev_series2), intent(in)  :: c
  real(kind=MATHpr),            intent(in)  :: xy0(2)
  real(kind=MATHpr),            intent(in)  :: xy1(2)
  type(type_chebyshev_series2), intent(out) :: s
  real(kind=MATHpr)                         :: tmp(c%degr(2)+1,c%degr(1)+1)
  integer                                   :: k

  call reset_chebyshev_series2( s, c%degr, c%dim )
  do k = 1,c%dim
     call chgvar( &
          c%coef(1:c%degr(1)+1,1:c%degr(2)+1,k), &
          s%coef(1:c%degr(1)+1,1:c%degr(2)+1,k), &
          xy0(1), &
          xy1(1), &
          c%degr(1), &
          c%degr(2)+1 )

     call chgvar( &
          transpose( s%coef(1:c%degr(1)+1,1:c%degr(2)+1,k) ), &
          tmp, &
          xy0(2), &
          xy1(2), &
          c%degr(2), &
          c%degr(1)+1 )

     s%coef(1:c%degr(1)+1,1:c%degr(2)+1,k) = transpose(tmp)
  end do

end subroutine chgvar2
