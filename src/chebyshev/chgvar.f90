subroutine chgvar( c, s, x0, x1, degr, dim )
  implicit none
  integer,           intent(in)        :: degr, dim
  real(kind=MATHpr), intent(in)        :: c(1:degr+1,dim)
  real(kind=MATHpr), intent(in)        :: x0, x1
  real(kind=MATHpr), intent(out)       :: s(1:degr+1,dim)
  real(kind=MATHpr)                    :: a, b, twob
  real(kind=MATHpr), dimension(degr+1) :: col_nm2, col_nm1, col_n
  integer                              :: n, j

  a = 0.5_MATHpr * ( x1 - x0 )
  twob = x1 + x0
  b = 0.5_MATHpr * twob

  col_nm2(:) = 0._MATHpr
  col_nm1(:) = 0._MATHpr
  col_n(:) = 0._MATHpr

  s(:,:) = 0._MATHpr
  col_nm2(1) = 1._MATHpr
  if ( degr == 0 ) then
     s(1,:) = c(1,:)
     return
  end if
  col_nm1(1:2) = [ b , a ]
  s(1,:) = c(1,:) + b * c(2,:)
  s(2,:) = a * c(2,:)

  do n = 3,degr+1
     col_n(1:n-1) = twob * col_nm1(1:n-1)
     col_n(1:n-2) = col_n(1:n-2) - col_nm2(1:n-2)
     col_n(2) = col_n(2) + a * col_nm1(1)
     col_n(2:n) = col_n(2:n) + a * col_nm1(1:n-1)
     col_n(1:n-2) = col_n(1:n-2) + a * col_nm1(2:n-1)

     do j = 1,dim
        s(1:n,j) = s(1:n,j) + c(n,j) * col_n(1:n)
     end do

     if ( n > degr ) return 

     col_nm2 = col_nm1
     col_nm1 = col_n
  end do

end subroutine chgvar
