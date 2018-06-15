module mod_chebyshev
  use mod_math

  implicit none


  type type_chebyshev_series1
     integer                        :: degr = 0
     integer                        :: dim = 0
     real(kind=MATHpr), allocatable :: coef(:,:)
  end type type_chebyshev_series1

  type type_chebyshev_series2
     integer                        :: degr(2) = 0
     integer                        :: dim = 0
     real(kind=MATHpr), allocatable :: coef(:,:,:)
  end type type_chebyshev_series2


  type ptr_chebyshev_series1
     type(type_chebyshev_series1), pointer :: ptr => null()
  end type ptr_chebyshev_series1

  type ptr_chebyshev_series2
     type(type_chebyshev_series2), pointer :: ptr => null()
  end type ptr_chebyshev_series2


  interface read_cs
     module procedure read_chebyshev_series1, read_chebyshev_series2
  end interface read_cs

  interface reset_cs
     module procedure reset_chebyshev_series1, reset_chebyshev_series2
  end interface reset_cs

  interface delete_cs
     module procedure delete_chebyshev_series1, delete_chebyshev_series2
  end interface delete_cs

  interface chebval
     module procedure chebval1, chebval1_n, chebval2, chebval2_n, chebval2_mn
  end interface chebval

contains



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



subroutine chebval1( &
     f, &
     c, &
     x )
  implicit none
  real(kind=MATHpr),            intent(in)  :: x
  type(type_chebyshev_series1), intent(in)  :: c
  real(kind=MATHpr),            intent(out) :: f(c%dim)

  call clenshaw( &
       f, &
       c%coef(1:c%degr+1,:), &
       [x], &
       1, &
       c%degr, &
       c%dim )

end subroutine chebval1



subroutine economization1( &
     c, &
     tol )
  implicit none
  type(type_chebyshev_series1), intent(inout) :: c
  real(kind=MATHpr),            intent(in)    :: tol
  real(kind=MATHpr)                           :: tolsqr
  integer                                     :: i, M

  tolsqr = tol**2

  M = c%degr+1
  do i = M,1,-1
     if ( sum( c%coef(i:M,:)**2 ) > tolsqr ) exit
  end do
  c%degr = i-1

end subroutine economization1



subroutine chebdiff( &
     c, &
     d, &
     degr, &
     dim )
  implicit none
  integer,           intent(in)  :: degr, dim
  real(kind=MATHpr), intent(in)  :: c(degr+1,dim)
  real(kind=MATHpr), intent(out) :: d(max(degr,1),dim)
  integer                        :: i

  d(:,:) = 0._MATHpr
  if ( degr < 1 ) return

  d(degr,:) = real( 2*degr, kind=MATHpr ) * c(degr+1,:)

  if ( degr > 1 ) then
     d(degr-1,:) = real( 2*(degr-1), kind=MATHpr ) * c(degr,:)

     if ( degr > 2 ) then
        do i = degr-2,1,-1
           d(i,:) = real( 2*i, kind=MATHpr ) * c(i+1,:) + d(i+2,:)
        end do
     end if

  end if
  d(1,:) = 0.5_MATHpr * d(1,:)

end subroutine chebdiff



subroutine ifcht2_padded( &
     c, &
     f, &
     m, &
     n )
  implicit none
  type(type_chebyshev_series2), intent(in)  :: c
  integer,                      intent(in)  :: m, n
  real(kind=MATHpr),            intent(out) :: f(m,n,c%dim)
  real(kind=MATHpr)                         :: w(3*max(m,n)+15)
  integer                                   :: i, j, k, r, s

  r = min( m, c%degr(1)+1 )
  s = min( n, c%degr(2)+1 )

  f(1:r,1:s,:) = c%coef(1:r,1:s,:)
  if ( r < m ) f(r+1:m,:,:) = 0._MATHpr
  if ( s < n ) f(1:r,s+1:n,:) = 0._MATHpr

  f(2:m-1,:,:) = 0.5_MATHpr * f(2:m-1,:,:)
  call dcosti( m, w )
  do k = 1,c%dim
     do j = 1,n
        call dcost( m, f(:,j,k), w )
     end do
  end do

  f(:,2:n-1,:) = 0.5_MATHpr * f(:,2:n-1,:)
  if ( m /= n ) call dcosti( n, w )
  do k = 1,c%dim
     do i = 1,m
        call dcost( n, f(i,:,k), w )
     end do
  end do

end subroutine ifcht2_padded



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



subroutine chebbivar2univar( &
     c2, &
     m, &
     n, &
     p, &
     c1, &
     ivar, &
     ival ) 
  use mod_constants
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



subroutine ifcht1( &
     c, &
     f )
  implicit none
  type(type_chebyshev_series1), intent(in)  :: c
  real(kind=MATHpr),            intent(out) :: f(c%degr+1,c%dim)
  real(kind=MATHpr), allocatable            :: w(:)    
  integer                                   :: j, degrprev, p
  save w, degrprev
  data degrprev /-1/

  f = c%coef(1:c%degr+1,:)
  f(2:c%degr,:) = 0.5_MATHpr * f(2:c%degr,:)

  p = sizeminw( c%degr+1 )
  if ( allocated(w) ) then
     if ( size(w) < p ) deallocate( w )
  end if
  if ( .not.allocated( w ) ) allocate( w(p) )

  if ( degrprev /= c%degr ) call dcosti( c%degr+1, w )
  degrprev = c%degr

  do j = 1,c%dim
     call dcost( c%degr+1, f(:,j), w )
  end do

end subroutine ifcht1



subroutine cgl_points( &
     x, &
     a, &
     b, &
     n )
    ! Returns the n Chebyshev-Gauss-Lobatto points
    ! mapped to the interval [a,b]
    implicit none
    integer,           intent(in)  :: n
    real(kind=MATHpr), intent(in)  :: a, b
    real(kind=MATHpr), intent(out) :: x(n)
    integer                        :: i = 0

    if (n == 1) then
       x = 0.5_MATHpr*(a+b)
    else
       x = a + (b-a) * 0.5_MATHpr * ( 1._MATHpr + &
            cos(MATHpi* [( real(i,kind=MATHpr)/real(n-1,kind=MATHpr) , &
            i = n-1,0,-1 )]) )
    end if

  end subroutine cgl_points



subroutine read_chebyshev_series2( &
     c, &
     filename )
  use mod_util
  implicit none
  character(*),                 intent(in)  :: filename
  type(type_chebyshev_series2), intent(out) :: c
  integer                                   :: degr(2), dim
  integer                                   :: fileunit, icoef, jcoef, idim

  call get_free_unit( fileunit )
  if ( fileunit == 0 ) STOP "read_chebyshev_series2 : could not find free unit"

  open( unit = fileunit, file = filename, action = "read" )
  read (fileunit,*) degr, dim
  degr = degr - 1

  if ( any(degr < 0) ) STOP "read_chebyshev_series2 : degr < 0"

  call reset_chebyshev_series2( &
       c, &
       degr, &
       dim)

  do idim = 1,dim
     do jcoef = 1,degr(2)+1
        do icoef = 1,degr(1)+1
           read (fileunit,*) c%coef(icoef,jcoef,idim)
        end do
     end do
  end do
  close(fileunit)

end subroutine read_chebyshev_series2



subroutine chebdiff1( &
     c, &
     d )
  implicit none
  type(type_chebyshev_series1), intent(in)  :: c
  type(type_chebyshev_series1), intent(out) :: d

  call reset_chebyshev_series1( d, max( c%degr-1, 0 ), c%dim )
  call chebdiff( c%coef(1:c%degr+1,:), d%coef(1:max(c%degr,1),:), c%degr, c%dim )

end subroutine chebdiff1



subroutine write_chebyshev_series2( &
     c, &
     filename )
  use mod_util
  implicit none
  character(*),                 intent(in) :: filename
  type(type_chebyshev_series2), intent(in) :: c
  integer                                  :: fileunit, icoef, jcoef, idim

  call get_free_unit( fileunit )
  if ( fileunit == 0 ) STOP "write_chebyshev_series2 : could not find free unit"

  open( unit = fileunit, file = filename, action = "write" )
  write (fileunit,*) c%degr+1, c%dim

  do idim = 1,c%dim
     do jcoef = 1,c%degr(2)+1
        do icoef = 1,c%degr(1)+1
           write (fileunit,"(ES22.15)") c%coef(icoef,jcoef,idim)
        end do
     end do
  end do
  close(fileunit)
end subroutine write_chebyshev_series2



subroutine delete_chebyshev_series1( &
     c )
  implicit none
  type(type_chebyshev_series1), intent(inout) :: c

  if (allocated(c%coef)) deallocate(c%coef)
  c%dim = 0
  c%degr = 0

end subroutine delete_chebyshev_series1



subroutine write_chebyshev_series1( &
     c, &
     filename )
  use mod_util
  implicit none
  character(*),                 intent(in) :: filename
  type(type_chebyshev_series1), intent(in) :: c
  integer                                  :: fileunit, icoef, idim

  call get_free_unit( fileunit )
  if ( fileunit == 0 ) STOP "write_chebyshev_series1 : could not find free unit"

  open( unit = fileunit, file = filename, action = "write" )
  write (fileunit,*) c%degr+1, c%dim

  do idim = 1,c%dim
     do icoef = 1,c%degr+1
        write (fileunit,"(ES22.15)") c%coef(icoef,idim)
     end do
  end do
  close(fileunit)
end subroutine write_chebyshev_series1



function ith_cgl_point( i, n ) result( x )
  ! Returns the i-th of the n+1 Chebyshev-Gauss-Lobatto points
  ! xi = cos( i\pi/n ) with 0 <= i <= n
  ! /!\ CGL points are ordered from +1 to -1 (x0 = +1, ..., xn = -1)
  implicit none
  integer, intent(in) :: i, n
  real(kind=MATHpr)   :: x

  x = cos( real( i, kind=MATHpr ) * MATHpi / real( n, kind=MATHpr ) )

end function ith_cgl_point



subroutine chebval2_corners( &
     c, &
     f )
  implicit none
  type(type_chebyshev_series2), intent(in)              :: c
  real(kind=MATHpr),            intent(out)             :: f(c%dim,4)
  integer                                               :: i = 0, m, n
  real(kind=MATHpr), dimension(c%degr(1)+1,c%degr(2)+1) :: om, on

  om = spread( &
       source = real( (-1)**[ ( i, i=0,c%degr(1) ) ], kind=MATHpr ), &
       dim = 2, &
       ncopies = c%degr(2)+1 )
  on = spread( &
       source = real( (-1)**[ ( i, i=0,c%degr(2) ) ], kind=MATHpr ), &
       dim = 1, &
       ncopies = c%degr(1)+1 )

  m = c%degr(1)+1
  n = c%degr(2)+1

  f(:,:) = 0._MATHpr
  do i = 1,c%dim
     f(i,1) = sum( sum( c%coef(1:m,1:n,i) * om * on, 2 ), 1 )
     f(i,2) = sum( sum( c%coef(1:m,1:n,i) * on, 2 ), 1 )
     f(i,3) = sum( sum( c%coef(1:m,1:n,i) * om, 2 ), 1 )
     f(i,4) = sum( sum( c%coef(1:m,1:n,i), 2 ), 1 )
  end do

end subroutine chebval2_corners



subroutine get_ch2be_from_collection( &
     collection, &
     N, &
     c2b )
  use mod_math
  implicit none
  integer,            intent(in)    :: N
  type(type_matrix),  intent(inout) :: collection(:)
  real(kind=MATHpr)                 :: c2b(N,N)

  if ( N > size(collection) ) then
     call ch2be_matrix( c2b, N-1 )
  else
     if ( .not.allocated(collection(N)%mat) ) then
        allocate( collection(N)%mat(N,N) )
        call ch2be_matrix( collection(N)%mat, N-1 )
     end if
     c2b = collection(N)%mat
  end if

end subroutine get_ch2be_from_collection



subroutine fcht2( &
     f, &
     c, &
     m, &
     n, &
     dim, &
     epsilon )
  implicit none
  integer,                      intent(in)           :: m, n, dim
  real(kind=MATHpr),            intent(in)           :: f(m, n, dim)
  type(type_chebyshev_series2), intent(out)          :: c
  real(kind=MATHpr),            intent(in), optional :: epsilon
  real(kind=MATHpr)                                  :: w(3*max(m,n)+15)
  integer                                            :: degr(2), i, j, k

  degr(1) = m - 1
  degr(2) = n - 1

  call reset_chebyshev_series2( c, degr, dim )

  c%coef = f / real( degr(1), kind=MATHpr )
  call dcosti( m, w )
  do k = 1,dim
     do j = 1,n
        call dcost( m, c%coef(:,j,k), w )
     end do
  end do
  c%coef([1,m],:,:) = 0.5_MATHpr * c%coef([1,m],:,:)

  c%coef = c%coef / real( degr(2), kind=MATHpr )
  if ( degr(1) /= degr(2) ) call dcosti( n, w )
  do k = 1,dim
     do i = 1,m
        call dcost( n, c%coef(i,:,k), w )
     end do
  end do
  c%coef(:,[1,n],:) = 0.5_MATHpr * c%coef(:,[1,n],:)

  if ( present(epsilon) ) then
     where ( abs(c%coef) < epsilon ) c%coef = 0._MATHpr
  end if

end subroutine fcht2



subroutine ifcht1_padded( &
     c, &
     f, &
     m )
  implicit none
  type(type_chebyshev_series1), intent(in)  :: c
  integer,                      intent(in)  :: m
  real(kind=MATHpr),            intent(out) :: f(m,c%dim)
  real(kind=MATHpr), allocatable            :: w(:)    
  integer                                   :: j, degrprev, r, p
  save w, degrprev
  data degrprev /0/

  r = min( m, c%degr+1 )
  f(1:r,:) = c%coef(1:r,:)
  if ( r < m ) f(r+1:m,:) = 0._MATHpr

  p = sizeminw( m )
  if ( allocated(w) ) then
     if ( size(w) < p ) deallocate( w )
  end if
  if ( .not.allocated( w ) ) allocate( w(p) )

  if ( degrprev /= m-1 ) call dcosti( m, w )
  degrprev = m-1

  f(2:m-1,:) = 0.5_MATHpr * f(2:m-1,:)
  do j = 1,c%dim
     call dcost( m, f(:,j), w )
  end do

end subroutine ifcht1_padded



subroutine delete_chebyshev_series2( &
     c )
  implicit none
  type(type_chebyshev_series2), intent(inout) :: c

  if (allocated(c%coef)) deallocate(c%coef)
  c%dim = 0
  c%degr(:) = 0

end subroutine delete_chebyshev_series2



subroutine chebfit1( &
     x, &
     y, &
     n, &
     dim, &
     degr, &
     c )
  implicit none
  integer,                      intent(in)  :: n, dim
  real(kind=MATHpr),            intent(in)  :: x(n)
  real(kind=MATHpr),            intent(in)  :: y(n,dim)
  integer,                      intent(in)  :: degr
  type(type_chebyshev_series1), intent(out) :: c
  real(kind=MATHpr)                         :: F(n,degr+1), twox(n), Ft(degr+1,n)
  integer                                   :: i, j

  F(:,:) = 0._MATHpr
  F(:,1) = 1._MATHpr
  F(:,2) = x

  if ( degr > 1 ) twox = 2._MATHpr * x

  do i = 3,degr+1
     F(:,i) = twox * F(:,i-1) - F(:,i-2)
  end do
  Ft = transpose(F)

  call reset_chebyshev_series1( c, degr, dim )

  do j = 1,dim
     call solve_NxN( c%coef(1:degr+1,j), matmul( Ft, F ), matmul( Ft, y(:,j) ) )
  end do

end subroutine chebfit1



subroutine chebval2_n( &
     f, &
     c, &
     xy, &
     n )
  implicit none
  integer,                      intent(in)  :: n
  real(kind=MATHpr),            intent(in)  :: xy(2,n)
  type(type_chebyshev_series2), intent(in)  :: c
  real(kind=MATHpr),            intent(out) :: f(c%dim,n)
  integer                                   :: i

  do i = 1,n
     call chebval2( &
          f(:,i), &
          c, &
          xy(:,i) )
  end do

end subroutine chebval2_n



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



subroutine reset_chebyshev_series1( &
     c, &
     degr, &
     dim )
  implicit none
  integer,                      intent(in)  :: degr
  integer,                      intent(in)  :: dim
  type(type_chebyshev_series1), intent(out) :: c

  if ( allocated(c%coef) ) then
     !if ( c%degr < degr .or. c%dim < dim ) deallocate(c%coef)
     if ( size(c%coef,1) <= degr .or. &
          size(c%coef,2) < dim ) deallocate(c%coef)
  end if
  if ( .not.allocated(c%coef) ) allocate(c%coef(1:degr+1, dim))

  c%degr = degr
  c%dim = dim
  c%coef(:,:) = 0._MATHpr

end subroutine reset_chebyshev_series1



subroutine fcht1( &
     f, &
     c, &
     m, &
     dim, &
     epsilon )
  ! Performs Fast Cosine Transform (using dfftpack routines)
  implicit none
  integer,                      intent(in)           :: m, dim
  real(kind=MATHpr),            intent(in)           :: f(m,dim)
  type(type_chebyshev_series1), intent(out)          :: c
  real(kind=MATHpr),            intent(in), optional :: epsilon
  real(kind=MATHpr), allocatable                     :: w(:)
  integer                                            :: j, p, degrprev
  save w, degrprev
  data degrprev /-1/

  call reset_chebyshev_series1( c, m-1, dim )

  c%coef = f / real( m-1, kind=MATHpr )

  p = sizeminw( m )
  if ( allocated(w) ) then
     if ( size(w) < p ) deallocate( w )
  end if
  if ( .not.allocated( w ) ) allocate( w(p) )

  if ( degrprev /= m-1 ) call dcosti( m, w(1:p) )
  degrprev = m-1

  do j = 1,dim
     call dcost( m, c%coef(:,j), w(1:p) )
  end do

  c%coef([1,m],:) = 0.5_MATHpr * c%coef([1,m],:)

  if ( present(epsilon) ) then
     where ( abs(c%coef) < epsilon ) c%coef = 0._MATHpr
  end if

end subroutine fcht1



function sizeminw( n ) 
  ! Minimal size of work array (see dfftpack documentation, dcosti subroutine)
  implicit none
  integer, intent(in) :: n
  integer             :: sizeminw
  sizeminw = 3*n + 15
end function sizeminw



subroutine chgvar2( &
     c, &
     s, &
     xy0, &
     xy1 )
  implicit none
  type(type_chebyshev_series2), intent(in)  :: c
  real(kind=MATHpr),            intent(in)  :: xy0(2)
  real(kind=MATHpr),            intent(in)  :: xy1(2)
  type(type_chebyshev_series2), intent(out) :: s
  real(kind=MATHpr)                         :: tmp(c%degr(2)+1,c%degr(1)+1)
  integer                                   :: k

  call reset_chebyshev_series2( s, c%degr, c%dim )
  do k = 1,c%dim
     call chgvar( &
          c%coef(1:c%degr(1)+1,1:c%degr(2)+1,k), &
          s%coef(1:c%degr(1)+1,1:c%degr(2)+1,k), &
          xy0(1), &
          xy1(1), &
          c%degr(1), &
          c%degr(2)+1 )

     call chgvar( &
          transpose( s%coef(1:c%degr(1)+1,1:c%degr(2)+1,k) ), &
          tmp, &
          xy0(2), &
          xy1(2), &
          c%degr(2), &
          c%degr(1)+1 )

     s%coef(1:c%degr(1)+1,1:c%degr(2)+1,k) = transpose(tmp)
  end do

end subroutine chgvar2



subroutine economization2( &
     c, &
     tol )
  implicit none
  type(type_chebyshev_series2), intent(inout) :: c
  real(kind=MATHpr),            intent(in)    :: tol
  real(kind=MATHpr)                           :: tolsqr
  real(kind=MATHpr), dimension(c%dim)         :: s, r
  integer                                     :: i, j

  tolsqr = tol**2

  s = sum( sum( abs(c%coef), 2), 1 )

  do i = c%degr(1)+1,1,-1
     r = sum( sum( abs(c%coef(1:i-1,1:c%degr(2)+1,:)), 2), 1 )
     if ( sum( (s - r)**2 ) > tolsqr ) exit
  end do
  c%degr(1) = i-1

  do j = c%degr(2)+1,1,-1
     r = sum( sum( abs(c%coef(1:i,1:j-1,:)), 2), 1 )
     if ( sum( (s - r)**2 ) > tolsqr ) exit
  end do
  c%degr(2) = j-1

end subroutine economization2



subroutine reset_chebyshev_series2( &
     c, &
     degr, &
     dim )
  implicit none
  integer,                      intent(in)  :: degr(2)
  integer,                      intent(in)  :: dim
  type(type_chebyshev_series2), intent(out) :: c

  if ( allocated(c%coef) ) then
     !if ( any(c%degr < degr) .or. c%dim < dim ) deallocate(c%coef)
     if ( size(c%coef,1) <= degr(1) .or. &
          size(c%coef,2) <= degr(2) .or. &
          size(c%coef,3) < dim ) deallocate(c%coef)
  end if
  if ( .not.allocated(c%coef) ) allocate(c%coef(1:degr(1)+1, 1:degr(2)+1, dim))

  c%degr = degr
  c%dim = dim
  c%coef(:,:,:) = 0._MATHpr

end subroutine reset_chebyshev_series2



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



subroutine chebval2( &
     f, &
     c, &
     xy )
  implicit none
  real(kind=MATHpr),            intent(in)  :: xy(2)
  type(type_chebyshev_series2), intent(in)  :: c
  real(kind=MATHpr),            intent(out) :: f(c%dim)
  real(kind=MATHpr)                         :: b(c%degr(2)+1,c%dim)
  integer                                   :: k

  do k = 1,c%dim
     call clenshaw( &
          b(:,k), &
          c%coef(1:c%degr(1)+1,1:c%degr(2)+1,k), &
          [xy(1)], &
          1, &
          c%degr(1), &
          c%degr(2)+1 )
  end do
  call clenshaw( &
       f, &
       b, &
       [xy(2)], &
       1, &
       c%degr(2), &
       c%dim )

end subroutine chebval2



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



subroutine chgvar1( &
     c, &
     s, &
     x0, &
     x1 )
  implicit none
  type(type_chebyshev_series1), intent(in)  :: c
  real(kind=MATHpr),            intent(in)  :: x0
  real(kind=MATHpr),            intent(in)  :: x1
  type(type_chebyshev_series1), intent(out) :: s

  call reset_chebyshev_series1( s, c%degr, c%dim )

  call chgvar( &
       c%coef(1:c%degr+1,:), &
       s%coef(1:c%degr+1,:), &
       x0, &
       x1, &
       c%degr, &
       c%dim )

end subroutine chgvar1



subroutine read_chebyshev_series1( &
     c, &
     filename )
  use mod_util
  implicit none
  character(*),                 intent(in)  :: filename
  type(type_chebyshev_series1), intent(out) :: c
  integer                                   :: degr, dim
  integer                                   :: fileunit, icoef, idim

  call get_free_unit( fileunit )
  if ( fileunit == 0 ) STOP "read_chebyshev_series1 : could not find free unit"

  open( unit = fileunit, file = filename, action = "read" )
  read (fileunit,*) degr, dim
  degr = degr - 1

  if ( degr < 0 ) STOP "read_chebyshev_series1 : degr < 0"

  call reset_chebyshev_series1( &
       c, &
       degr, &
       dim)

  do idim = 1,dim
     do icoef = 1,degr+1
        read (fileunit,*) c%coef(icoef,idim)
     end do
  end do
  close(fileunit)

end subroutine read_chebyshev_series1



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



subroutine chebval1_n( &
     f, &
     c, &
     x, &
     n )
  implicit none
  integer,                      intent(in)  :: n
  real(kind=MATHpr),            intent(in)  :: x(n)
  type(type_chebyshev_series1), intent(in)  :: c
  real(kind=MATHpr),            intent(out) :: f(n,c%dim)

  call clenshaw( &
       f, &
       c%coef(1:c%degr+1,:), &
       x, &
       n, &
       c%degr, &
       c%dim )

end subroutine chebval1_n
end module mod_chebyshev
