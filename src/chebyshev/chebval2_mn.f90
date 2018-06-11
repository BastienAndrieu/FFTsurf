subroutine chebval2_mn( &
     f, &
     c, &
     x, &
     y, &
     nx, &
     ny )
  ! Evaluation of bivariate Chebyshev series at nodes of a tensor-product grid
  implicit none
  integer,                      intent(in)  :: nx, ny
  real(kind=MATHpr),            intent(in)  :: x(nx), y(ny)
  type(type_chebyshev_series2), intent(in)  :: c
  real(kind=MATHpr),            intent(out) :: f(nx,ny,c%dim)
  real(kind=MATHpr)                         :: b(ny,c%degr(1)+1)
  integer                                   :: k

  do k = 1,c%dim
     call clenshaw( &
          b, &
          transpose( c%coef(1:c%degr(1)+1,1:c%degr(2)+1,k) ), &
          y, &
          ny, &
          c%degr(2), &
          c%degr(1)+1 )

     call clenshaw( &
          f(:,:,k), &
          transpose(b), &
          x, &
          nx, &
          c%degr(1), &
          ny )
  end do

end subroutine chebval2_mn
