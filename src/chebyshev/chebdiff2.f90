subroutine chebdiff2( &
     c, &
     du, &
     dv )
  implicit none
  type(type_chebyshev_series2), intent(in)            :: c
  type(type_chebyshev_series2), intent(out)           :: du
  type(type_chebyshev_series2), intent(out), optional :: dv
  real(kind=MATHpr)                                   :: tmp(max(c%degr(2),1),max(c%degr(1)+1,1))
  integer                                             :: k

  call reset_chebyshev_series2( du, max( [c%degr(1)-1, c%degr(2)], 0 ), c%dim )
  do k = 1,c%dim
     call chebdiff( &
          c%coef(1:c%degr(1)+1,1:c%degr(2)+1,k), &
          du%coef(1:max(c%degr(1),1),1:c%degr(2)+1,k), &
          c%degr(1), &
          c%degr(2)+1 )
  end do

  if ( present(dv) ) then
     call reset_chebyshev_series2( dv, max( [c%degr(1), c%degr(2)-1], 0 ), c%dim )
     do k = 1,c%dim
        call chebdiff( &
             transpose( c%coef(1:c%degr(1)+1,1:c%degr(2)+1,k) ), &
             tmp, &
             c%degr(2), &
             c%degr(1)+1 )
        dv%coef(1:size(tmp,2),1:size(tmp,1),k) = transpose(tmp)
     end do
  end if

end subroutine chebdiff2
