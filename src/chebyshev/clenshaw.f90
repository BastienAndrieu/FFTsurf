subroutine clenshaw( &
     f, &
     c, &
     x, &
     nx, &
     degr, &
     dim )
  implicit none
  integer,                      intent(in)  :: nx, degr, dim
  real(kind=MATHpr),            intent(in)  :: x(nx)
  real(kind=MATHpr),            intent(in)  :: c(degr+1,dim)
  real(kind=MATHpr),            intent(out) :: f(nx,dim)
  real(kind=MATHpr), dimension(nx)          :: bim1, bi, bip1
  integer                                   :: i, j

  do j = 1,dim
     bip1(:) = 0._MATHpr
     bi(:) = 0._MATHpr
     do i = degr+1,2,-1
        bim1 = c(i,j) + 2._MATHpr * x * bi - bip1
        if ( i > 2 ) then
           if (i < degr+1) bip1 = bi
           bi = bim1
        end if
     end do
     f(:,j) = c(1,j) + x * bim1 - bi
  end do

end subroutine clenshaw
