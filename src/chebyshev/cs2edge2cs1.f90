subroutine cs2edge2cs1( &
     c2, &
     ivar, &
     ival, &
     c1 )
  implicit none
  type(type_chebyshev_series2), intent(in)  :: c2
  integer,                      intent(in)  :: ivar,ival
  type(type_chebyshev_series1), intent(out) :: c1
  integer                                   :: i, j, k=0

  call reset_chebyshev_series1( c1, c2%degr(1+mod(ivar,2)), c2%dim )

  if ( ival == 2 ) then
     do j = 1,c2%dim
        c1%coef(:,j) = sum( c2%coef(1:c2%degr(1)+1,1:c2%degr(2)+1,j), ivar )
     end do
  else
     if ( ivar == 1 ) then
        do j = 1,c2%dim
           do i = 1,c2%degr(2)+1
              c1%coef(i,j) = sum( &
                   c2%coef(1:c2%degr(1)+1,i,j) * real( [( (-1)**k, k=0,c2%degr(1) )], kind=MATHpr ) &
                   )
           end do
        end do
     else
        do j = 1,c2%dim
           do i = 1,c2%degr(1)+1
              c1%coef(i,j) = sum( &
                   c2%coef(i,1:c2%degr(2)+1,j) * real( [( (-1)**k, k=0,c2%degr(2) )], kind=MATHpr ) &
                   )
           end do
        end do
     end if
  end if

end subroutine cs2edge2cs1
