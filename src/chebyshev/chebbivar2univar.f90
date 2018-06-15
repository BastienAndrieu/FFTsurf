subroutine chebbivar2univar( &
     c2, &
     m, &
     n, &
     p, &
     c1, &
     ivar, &
     ival ) 
  implicit none
  integer,       intent(in)  :: m, n, p
  real(kind=fp), intent(in)  :: c2(m,n,p)
  integer,       intent(in)  :: ivar, ival
  real(kind=fp), intent(out) :: c1(size(c2,1+mod(ivar,2)),p)
  integer                    :: i, j, k

  if ( ival == 2 ) then

     do j = 1,p
        c1(:,j) = sum( c2(1:m,1:n,j), ivar )
     end do

  else

     if ( ivar == 1 ) then
        
        do j = 1,p
           do i = 1,n
              c1(i,j) = sum( c2(1:m,i,j) * real( [( (-1)**k, k=0,m-1 )], kind=fp ) )
           end do
        end do

     else

        do j = 1,p
           do i = 1,m
              c1(i,j) = sum( c2(i,1:n,j) * real( [( (-1)**k, k=0,n-1 )], kind=fp ) )
           end do
        end do

     end if

  end if

end subroutine chebbivar2univar
