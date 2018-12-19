module mod_chebyshev

  use mod_math

  implicit none

contains

  subroutine clenshaw( &
       f, &
       c, &
       x, &
       nx, &
       degr, &
       dim )
    implicit none
    integer,       intent(in)    :: nx, degr, dim
    real(kind=fp), intent(in)    :: x(nx)
    real(kind=fp), intent(in)    :: c(degr+1,dim)
    real(kind=fp), intent(out)   :: f(nx,dim)
    real(kind=fp), dimension(nx) :: bim1, bi, bip1
    integer                      :: i, j

    do j = 1,dim
       bip1(:) = 0._fp
       bi(:)   = 0._fp
       bim1(:) = 0._fp
       do i = degr+1,2,-1
          bim1 = c(i,j) + 2._fp * x * bi - bip1
          if ( i > 2 ) then
             if (i < degr+1) bip1 = bi
             bi = bim1
          end if
       end do
       f(:,j) = c(1,j) + x * bim1 - bi
    end do

  end subroutine clenshaw







  subroutine chebdiff( &
       c, &
       d, &
       degr, &
       dim )
    implicit none
    integer,       intent(in)  :: degr, dim
    real(kind=fp), intent(in)  :: c(degr+1,dim)
    real(kind=fp), intent(out) :: d(max(degr,1),dim)
    integer                    :: i

    d(:,:) = 0._fp
    if ( degr < 1 ) return

    d(degr,:) = real( 2*degr, kind=fp ) * c(degr+1,:)

    if ( degr > 1 ) then
       d(degr-1,:) = real( 2*(degr-1), kind=fp ) * c(degr,:)

       if ( degr > 2 ) then
          do i = degr-2,1,-1
             d(i,:) = real( 2*i, kind=fp ) * c(i+1,:) + d(i+2,:)
          end do
       end if

    end if
    d(1,:) = 0.5_fp * d(1,:)

  end subroutine chebdiff






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











  ! Direct Fast Chebyshev Transform
  subroutine fcht1( &
       f, &
       c, &
       m, &
       dim, &
       epsilon )
    implicit none
    integer,       intent(in)           :: m, dim
    real(kind=fp), intent(in)           :: f(m,dim)
    real(kind=fp), intent(out)          :: c(m,dim)
    real(kind=fp), intent(in), optional :: epsilon
    real(kind=fp), allocatable          :: w(:)
    integer                             :: mprev, p, j
    save w, mprev
    data mprev /0/

    call manage_wsave( &
         w, &
         p, &
         m, &
         mprev )

    c = f / real( m-1, kind=fp )

    do j = 1,dim
       call dcost( m, c(1:m,j), w(1:p) )
    end do

    c([1,m],1:dim) = 0.5_fp * c([1,m],1:dim)

    if ( present(epsilon) ) then
       where ( abs(c) < epsilon ) c = 0._fp
    end if

  end subroutine fcht1




  ! Inverse Fast Chebyshev Transform
  subroutine ifcht1( &
       c, &
       f, &
       m, &
       dim )
    implicit none
    integer,       intent(in)  :: m, dim
    real(kind=fp), intent(in)  :: c(:,:)
    real(kind=fp), intent(out) :: f(m,dim)
    real(kind=fp), allocatable :: w(:)
    integer                    :: mprev, r, p, j
    save w, mprev
    data mprev /0/

    call manage_wsave( &
         w, &
         p, &
         m, &
         mprev )

    f(1:m,1:dim) = 0._fp
    r = min( m, size(c,1) )
    f(1:r,1:dim) = c(1:r,1:dim)
    f(2:m-1,1:dim) = 0.5_fp * f(2:m-1,1:dim)

    do j = 1,dim
       call dcost( m, f(1:m,j), w(1:p) )
    end do

  end subroutine ifcht1





  ! Direct double Fast Chebyshev Transform
  subroutine fcht2( &
      f, &
      c, &
      m, &
      n, &
      dim, &
      epsilon )
    implicit none
    integer,       intent(in)           :: m, n, dim
    real(kind=fp), intent(in)           :: f(m,n,dim)
    real(kind=fp), intent(out)          :: c(m,n,dim)
    real(kind=fp), intent(in), optional :: epsilon
    real(kind=fp)                       :: w(3*max(m,n)+15)
    integer                             :: i, j, k

    call dcosti( m, w )
    c = f / real( m-1, kind=fp )
    do k = 1,dim
       do j = 1,n
          call dcost( m, c(1:m,j,k), w )
       end do
    end do
    c([1,m],1:n,1:dim) = 0.5_fp * c([1,m],1:n,1:dim)

    
    if ( m /= n ) call dcosti( n, w )
    c = c / real( n-1, kind=fp )
    do k = 1,dim
       do i = 1,m
          call dcost( n, c(i,1:n,k), w )
       end do
    end do
    c(1:m,[1,n],1:dim) = 0.5_fp * c(1:m,[1,n],1:dim)

    if ( present(epsilon) ) then
       where ( abs(c) < epsilon ) c = 0._fp
    end if

  end subroutine fcht2
  




  ! Inverse double Fast Chebyshev Transform
  subroutine ifcht2( &
       c, &
       f, &
       m, &
       n, &
       dim )
    implicit none
    integer,       intent(in)  :: m, n, dim
    real(kind=fp), intent(in)  :: c(:,:,:)
    real(kind=fp), intent(out) :: f(m,n,dim)
    real(kind=fp)              :: w(3*max(m,n)+15)
    integer                    :: i, j, k, r, s

    call dcosti( m, w )
    f(1:m,1:n,1:dim) = 0._fp
    r = min( m, size(c,1) )
    s = min( n, size(c,2) )
    f(1:r,1:s,1:dim) = c(1:r,1:s,1:dim)
    f(2:m-1,1:n,1:dim) = 0.5_fp * f(2:m-1,1:n,1:dim)

    do k = 1,dim
       do j = 1,n
          call dcost( m, f(1:m,j,k), w )
       end do
    end do

    if ( m /= n ) call dcosti( n, w )
    f(1:m,2:n-1,1:dim) = 0.5_fp * f(1:m,2:n-1,1:dim)
    do k = 1,dim
       do i = 1,m
          call dcost( n, f(i,1:n,k), w )
       end do
    end do

  end subroutine ifcht2
  














  function sizeminwsave( n ) 
    ! Minimal size of work array (see dfftpack documentation, dcosti subroutine)
    implicit none
    integer, intent(in) :: n
    integer             :: sizeminwsave
    sizeminwsave = 3*n + 15

  end function sizeminwsave



  subroutine manage_wsave( &
       w, &
       p, &
       m, &
       mprev )
    implicit none
    real(kind=fp), allocatable, intent(inout) :: w(:)
    integer,                    intent(out)   :: p
    integer,                    intent(in)    :: m
    integer,                    intent(inout) :: mprev

    p = sizeminwsave( m )
    if ( allocated(w) ) then
       if ( size(w) < p ) deallocate( w )
    end if
    if ( .not.allocated(w) ) allocate( w(p) )

    if ( mprev /= m ) call dcosti( m, w(1:p) )
    mprev = m

  end subroutine manage_wsave






  
  subroutine chgvar( &
       c, &
       s, &
       x0, &
       x1, &
       degr, &
       dim )
    implicit none
    integer,       intent(in)        :: degr, dim
    real(kind=fp), intent(in)        :: c(1:degr+1,dim)
    real(kind=fp), intent(in)        :: x0, x1
    real(kind=fp), intent(out)       :: s(1:degr+1,dim)
    real(kind=fp)                    :: a, b, twob
    real(kind=fp), dimension(degr+1) :: col_nm2, col_nm1, col_n
    integer                          :: n, j

    a = 0.5_fp * ( x1 - x0 )
    twob = x1 + x0
    b = 0.5_fp * twob

    col_nm2(:) = 0._fp
    col_nm1(:) = 0._fp
    col_n(:) = 0._fp

    s(:,:) = 0._fp
    col_nm2(1) = 1._fp
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



  
  subroutine cgl_nodes( &
       x, &
       N, &
       a, &
       b )
    implicit none
    integer, intent(in)                  :: N
    real(kind=fp),  intent(out)          :: x(N+1)
    real(kind=fp),  intent(in), optional :: a, b
    integer                              :: i = 0

    if ( N == 0 ) then
       x(1) = 0._fp
    else
       x(1:N+1) = cos(CSTpi * [(real(i, kind=fp) / real(N, kind=fp), i=0,N)])
    end if

    if ( present(a) .and. present(b) ) then
       x = 0.5_fp * ((a-b)*x + a + b)
    end if

  end subroutine cgl_nodes
  

end module mod_chebyshev
