subroutine berndiff2( &
     b, &
     du, &
     dv )
  implicit none
  type(type_bernstein_series2), intent(in)            :: b
  type(type_bernstein_series2), intent(out)           :: du
  type(type_bernstein_series2), intent(out), optional :: dv
  real(kind=MATHpr)                                   :: tmp(max(b%degr(2),1),max(b%degr(1)+1,1))
  integer                                             :: k

  call reset_bernstein_series2( du, max( [b%degr(1)-1, b%degr(2)], 0 ), b%dim )
  do k = 1,b%dim
     call berndiff( &
          b%coef(1:b%degr(1)+1,1:b%degr(2)+1,k), &
          du%coef(1:max(b%degr(1),1),1:b%degr(2)+1,k), &
          b%degr(1), &
          b%degr(2)+1 )
  end do

  if ( present(dv) ) then
     call reset_bernstein_series2( dv, max( [b%degr(1), b%degr(2)-1], 0 ), b%dim )
     do k = 1,b%dim
        call berndiff( &
             transpose( b%coef(1:b%degr(1)+1,1:b%degr(2)+1,k) ), &
             tmp, &
             b%degr(2), &
             b%degr(1)+1 )
        dv%coef(1:size(tmp,2),1:size(tmp,1),k) = transpose(tmp)
     end do
  end if

end subroutine berndiff2
