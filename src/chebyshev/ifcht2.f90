subroutine ifcht2( &
     c, &
     f )
  implicit none
  type(type_chebyshev_series2), intent(in)  :: c
  real(kind=MATHpr),            intent(out) :: f(c%degr(1)+1,c%degr(2)+1,c%dim)
  real(kind=MATHpr)                         :: w(3*(max(c%degr(1),c%degr(2))+1)+15)
  integer                                   :: i, j, k

  f = c%coef(1:c%degr(1)+1,1:c%degr(2)+1,:)

  f(2:c%degr(1),:,:) = 0.5_MATHpr * f(2:c%degr(1),:,:)
  call dcosti( c%degr(1)+1, w )
  do k = 1,c%dim
     do j = 1,c%degr(2)+1
        call dcost( c%degr(1)+1, f(:,j,k), w )
     end do
  end do

  f(:,2:c%degr(2),:) = 0.5_MATHpr * f(:,2:c%degr(2),:)
  if ( c%degr(2) /= c%degr(1) ) call dcosti( c%degr(2)+1, w )
  do k = 1,c%dim
     do i = 1,c%degr(1)+1
        call dcost( c%degr(2)+1, f(i,:,k), w )
     end do
  end do

end subroutine ifcht2
