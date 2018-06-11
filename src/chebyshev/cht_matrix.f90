subroutine cht_matrix( C, degr )
  ! Matrix for transformation from physical space to Chebyshev transform space 
  ! ( "Spectral Methods - Fundamentals in Single Domains", Canuto et al., p.86 )
  ! C_{ij} = \frac{2}{N \bar{c}_i \bar{c}_j} \cos \frac{\pi i j}{N} 
  ! \quad\text{with}\quad 
  ! \bar{c}_i = 
  ! \begin{cases}
  ! 2 & \text{ if } i = 0 \text{ or } N \\  
  ! 1 & \text{ else}.
  ! \end{cases}
  ! C_{ij} = C_{ji} = (-1)^j C_{(N-i)j} = (-1)^j C_{j(N-i)}
  implicit none
  integer,           intent(in)  :: degr
  real(kind=MATHpr), intent(out) :: C(1:degr+1,1:degr+1)
  integer                        :: m, p, i, j
  real(kind=MATHpr)              :: two_over_n, pi_over_n

  if ( degr == 0 ) then
     C(1,1) = 1._MATHpr
     return 
  end if

  m = degr / 2
  two_over_n = 1._MATHpr / real( degr, kind=MATHpr )
  pi_over_n = two_over_n * MATHpi
  two_over_n = 2._MATHpr * two_over_n

  C(:,:) = 0._MATHpr
  do j = 0,m
     do i = j,m
        C(1+i,1+j) = two_over_n * cos( pi_over_n * real(i*j, kind=MATHpr ) )
        if ( j == 0 ) C(1+i,1+j) = 0.5_MATHpr * C(1+i,1+j)
        if ( i == 0 ) C(1+i,1+j) = 0.5_MATHpr * C(1+i,1+j)
        if ( i > j ) C(1+j,1+i) = C(1+i,1+j) ! C is symmetric
     end do
  end do

  p = m
  if ( mod(degr,2) /= 0 ) p = p + 1

  C( m+2:degr+1:1, 1:m+1:2 ) = C( p:1:-1, 1:m+1:2 )
  C( m+2:degr+1:1, 2:m+1:2 ) = -C( p:1:-1, 2:m+1:2 )

  C( 1:degr+1:2, degr+1:m+2:-1 ) = C( 1:degr+1:2, 1:p:1 )
  C( 2:degr+1:2, degr+1:m+2:-1 ) = -C( 2:degr+1:2, 1:p:1 )

end subroutine cht_matrix
