subroutine ch2be_matrix( A, degr )
  ! Transformation matrix from Chebyshev to Bernstein polynomial basis 
  ! "Transformation of Chebyshev–Bernstein polynomial basis", Rababah (2003) -- p.8 (615)
  implicit none
  integer,           intent(in)  :: degr
  real(kind=MATHpr), intent(out) :: A(1:degr+1,1:degr+1)
  real(kind=MATHpr)              :: fnumerator, fdegr, s
  integer                        :: i, j, k, m, imin, p

  if ( degr == 0 ) then
     A(1,1) = 1._MATHpr
     return
  end if

  m = degr / 2

  A(1:m+1,1) = 1._MATHpr
  A(1,2:degr+1:2) = -1._MATHpr
  A(1,3:degr+1:2) = 1._MATHpr

  A(2:m+1,2:degr+1) = 0._MATHpr
  do k = 1,degr
     fnumerator = factln(2*k) + factln(degr-k)
     do j = 1,m
        imin = max(0,j+k-degr)
        s = real( (-1)**(k-imin), kind=MATHpr )
        do i = imin,min(j,k)
           A(1+j,1+k) = A(1+j,1+k) + s * exp( fnumerator - &
                (factln(2*(k-i)) + factln(2*i) + factln(degr-k-j+i) + factln(j-i)) &
                )
           s = -s
        end do
     end do
  end do

  fdegr = factln( degr )
  do j = 1,m
     A(1+j,2:degr+1) = A(1+j,2:degr+1) / exp( fdegr - factln(degr-j) - factln(j) )
  end do

  p = m
  if ( mod(degr,2) /= 0 ) p = p + 1

  A(degr+1:p+1:-1,2:degr+1:2) = -A(1:m+1,2:degr+1:2)
  A(degr+1:p+1:-1,1:degr+1:2) = A(1:m+1,1:degr+1:2)

end subroutine ch2be_matrix
