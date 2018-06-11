module mod_chebyshev_bak

  use mod_util
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!              RESET
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine reset_chebyshev_series1( &
       c, &
       degr, &
       dim)
    implicit none
    integer,                      intent(in)  :: degr
    integer,                      intent(in)  :: dim
    type(type_chebyshev_series1), intent(out) :: c

    if ( allocated(c%coef) ) then
       if ( c%degr < degr .or. c%dim < dim ) deallocate(c%coef)
    end if
    if ( .not.allocated(c%coef) ) allocate(c%coef(1:degr+1, dim))

    c%degr = degr
    c%dim = dim

  end subroutine reset_chebyshev_series1



  subroutine reset_chebyshev_series2( &
       c, &
       degr, &
       dim)
    implicit none
    integer,                      intent(in)  :: degr(2)
    integer,                      intent(in)  :: dim
    type(type_chebyshev_series2), intent(out) :: c

    if ( allocated(c%coef) ) then
       if ( any(c%degr < degr) .or. c%dim < dim ) deallocate(c%coef)
    end if
    if ( .not.allocated(c%coef) ) allocate(c%coef(1:degr(1)+1, 1:degr(2)+1, dim))
    c%degr = degr
    c%dim = dim

  end subroutine reset_chebyshev_series2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!              DELETE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine delete_chebyshev_series1( &
       c)
    implicit none
    type(type_chebyshev_series1), intent(inout) :: c

    if (allocated(c%coef)) deallocate(c%coef)
    c%dim = 0
    c%degr = 0

  end subroutine delete_chebyshev_series1



  subroutine delete_chebyshev_series2( &
       c)
    implicit none
    type(type_chebyshev_series2), intent(inout) :: c

    if (allocated(c%coef)) deallocate(c%coef)
    c%dim = 0
    c%degr(:) = 0

  end subroutine delete_chebyshev_series2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!              READ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_chebyshev_series1( &
       c, &
       filename)
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



  subroutine read_chebyshev_series2( &
       c, &
       filename)

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!              WRITE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_chebyshev_series1( &
       c, &
       filename)
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


  subroutine write_chebyshev_series2( &
       c, &
       filename)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!              FAST DIRECT TRANSFORM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sizeminw( n ) 
    ! Minimal size of work array (see dfftpack documentation, dcosti subroutine)
    implicit none
    integer, intent(in) :: n
    integer             :: sizeminw
    sizeminw = 3*n + 15
  end function sizeminw


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
    data degrprev /0/

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




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!              FAST INVERSE TRANSFORM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ifcht1( &
       c, &
       f )
    implicit none
    type(type_chebyshev_series1), intent(in)  :: c
    real(kind=MATHpr),            intent(out) :: f(c%degr+1,c%dim)
    real(kind=MATHpr), allocatable            :: w(:)    
    integer                                   :: j, degrprev, p
    save w, degrprev
    data degrprev /0/

    f = c%coef
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




  subroutine ifcht2( &
       c, &
       f )
    implicit none
    type(type_chebyshev_series2), intent(in)  :: c
    real(kind=MATHpr),            intent(out) :: f(c%degr(1)+1,c%degr(2)+1,c%dim)
    real(kind=MATHpr)                         :: w(3*(max(c%degr(1),c%degr(2))+1)+15)
    integer                                   :: i, j, k

    f = c%coef

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!              EVAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
         c%coef, &
         [x], &
         1, &
         c%degr, &
         c%dim )

  end subroutine chebval1


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
         c%coef, &
         x, &
         n, &
         c%degr, &
         c%dim )

  end subroutine chebval1_n



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
       call clenshaw( b(:,k), c%coef(:,:,k), [xy(1)], 1, c%degr(1), c%degr(2)+1 )
    end do
    call clenshaw( f, b, [xy(2)], 1, c%degr(2), c%dim )

  end subroutine chebval2



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
       call chebval2( f(:,i), c, xy(:,i) )
    end do

  end subroutine chebval2_n


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
            transpose( c%coef(1:c%degr(1),1:c%degr(2),k) ), &
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



  subroutine chebval2_corners( c, f )
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!              DIFF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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



  subroutine chebdiff1( &
       c, &
       d )
    implicit none
    type(type_chebyshev_series1), intent(in)  :: c
    type(type_chebyshev_series1), intent(out) :: d

    call reset_chebyshev_series1( d, max( c%degr-1, 0 ), c%dim )
    call chebdiff( c%coef, d%coef, c%degr, c%dim )

  end subroutine chebdiff1



  subroutine chebdiff2( &
       c, &
       du, &
       dv )
    implicit none
    type(type_chebyshev_series2), intent(in)            :: c
    type(type_chebyshev_series2), intent(out)           :: du
    type(type_chebyshev_series2), intent(out), optional :: dv
    real(kind=MATHpr)                                   :: tmp(max(c%degr(2),1),c%degr(1)+1)
    integer                                             :: k

    call reset_chebyshev_series2( du, max( [c%degr(1)-1, c%degr(2)], 0 ), c%dim )
    do k = 1,c%dim
       call chebdiff( &
            c%coef(:,:,k), &
            du%coef(:,:,k), &
            c%degr(1), &
            c%degr(2)+1 )
    end do

    if ( present(dv) ) then
       call reset_chebyshev_series2( dv, max( [c%degr(1), c%degr(2)-1], 0 ), c%dim )
       do k = 1,c%dim
          call chebdiff( &
               transpose( c%coef(:,:,k) ), &
               tmp, &
               c%degr(2), &
               c%degr(1)+1 )
          dv%coef(:,:,k) = transpose(tmp)
       end do
    end if

  end subroutine chebdiff2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!              LINEAR CHANGE OF VARIABLES ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    call chgvar( c%coef, s%coef, x0, x1, c%degr, c%dim )

  end subroutine chgvar1



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
            c%coef(:,:,k), &
            s%coef(:,:,k), &
            xy0(1), &
            xy1(1), &
            c%degr(1), &
            c%degr(2)+1 )

       call chgvar( &
            transpose( s%coef(:,:,k) ), &
            tmp, &
            xy0(2), &
            xy1(2), &
            c%degr(2), &
            c%degr(1)+1 )

       s%coef(:,:,k) = transpose(tmp)
    end do

  end subroutine chgvar2


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
          c1%coef(:,j) = sum( c2%coef(:,:,j), ivar )
       end do
    else
       if ( ivar == 1 ) then
          do j = 1,c2%dim
             do i = 1,c2%degr(2)+1
                c1%coef(i,j) = sum( &
                     c2%coef(:,i,j) * real( [( (-1)**k, k=0,c2%degr(1) )], kind=MATHpr ) &
                     )
             end do
          end do
       else
          do j = 1,c2%dim
             do i = 1,c2%degr(1)+1
                c1%coef(i,j) = sum( &
                     c2%coef(i,:,j) * real( [( (-1)**k, k=0,c2%degr(2) )], kind=MATHpr ) &
                     )
             end do
          end do
       end if
    end if

  end subroutine cs2edge2cs1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cgl_points(x, a, b, n)
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


  function ith_cgl_point( i, n ) result( x )
    ! Returns the i-th of the n+1 Chebyshev-Gauss-Lobatto points
    ! xi = cos( i\pi/n ) with 0 <= i <= n
    ! /!\ CGL points are ordered from +1 to -1 (x0 = +1, ..., xn = -1)
    implicit none
    integer, intent(in) :: i, n
    real(kind=MATHpr)   :: x
    
    x = cos( real( i, kind=MATHpr ) * MATHpi / real( n, kind=MATHpr ) )

  end function ith_cgl_point


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!              TRANSFORMATION MATRICES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine zero_padding1( c, degr )
    ! Padds with zeros the coeffcient vector of a univariate Chebyshev series
    implicit none
    integer,                      intent(in)    :: degr
    type(type_chebyshev_series1), intent(inout) :: c
    real(kind=MATHpr)                           :: tmp(degr+1,c%dim)
    integer                                     :: dim

    if ( degr <= c%degr ) return
    tmp(:,:) = 0._MATHpr
    tmp(1:c%degr+1,:) = c%coef
    dim = c%dim
    call reset_chebyshev_series1( c, degr, dim )
    c%coef = tmp

  end subroutine zero_padding1


  subroutine zero_padding2( c, degr )
    ! Padds with zeros the coeffcient matrix of a bivariate Chebyshev series
    implicit none
    integer,                      intent(in)    :: degr(2)
    type(type_chebyshev_series2), intent(inout) :: c
    real(kind=MATHpr)                           :: tmp(degr(1)+1,degr(2)+1,c%dim)
    integer                                     :: dim

    if ( any(degr <= c%degr) ) return
    tmp(:,:,:) = 0._MATHpr
    tmp(1:c%degr(1)+1,1:c%degr(2)+1,:) = c%coef
    dim = c%dim
    call reset_chebyshev_series2( c, degr, dim )
    c%coef = tmp

  end subroutine zero_padding2



  subroutine economization1( c, tol )
    implicit none
    type(type_chebyshev_series1), intent(inout)        :: c
    real(kind=MATHpr),            intent(in), optional :: tol
    real(kind=MATHpr)                                  :: tol2
    integer                                            :: i, M
    
    if ( present(tol) ) then
       tol2 = tol**2
    else
       tol2 = 0._MATHpr
    end if

    M = c%degr+1
    do i = M,1,-1
       if ( sum( c%coef(i:M,:)**2 ) > tol2 ) exit
    end do
    c%degr = i-1

  end subroutine economization1



  subroutine economization2( c, tol )
    implicit none
    type(type_chebyshev_series2), intent(inout)        :: c
    real(kind=MATHpr),            intent(in), optional :: tol
    real(kind=MATHpr)                                  :: tol2
    real(kind=MATHpr), dimension(c%dim)                :: s, r
    integer                                            :: i, j
    
    if ( present(tol) ) then
       tol2 = tol**2
    else
       tol2 = 0._MATHpr
    end if

    s = sum( sum( abs(c%coef), 2), 1 )

    do i = c%degr(1)+1,1,-1
       r = sum( sum( abs(c%coef(1:i-1,:,:)), 2), 1 )
       if ( sum( (s - r)**2 ) > tol2 ) exit
    end do
    c%degr(1) = i-1

    do j = c%degr(2)+1,1,-1
       r = sum( sum( abs(c%coef(1:i,1:j-1,:)), 2), 1 )
       if ( sum( (s - r)**2 ) > tol2 ) exit
    end do
    c%degr(2) = j-1
    
  end subroutine economization2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!              ARRAY TO CS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine array2chebyshev_series1( &
     a, &
     m, &
     n, &
     c )
  implicit none
  integer,                      intent(in)  :: m, n
  real(kind=MATHpr),            intent(in)  :: a(m,n)
  type(type_chebyshev_series1), intent(out) :: c

  call reset_chebyshev_series1( c, m-1, n )
  c%coef(1:m,1:n) = a(1:m,1:n)

end subroutine array2chebyshev_series1



subroutine array2chebyshev_series2( &
     a, &
     m, &
     n, &
     p, &
     c )
  implicit none
  integer,                      intent(in)  :: m, n, p
  real(kind=MATHpr),            intent(in)  :: a(m,n,p)
  type(type_chebyshev_series2), intent(out) :: c

  call reset_chebyshev_series2( c, [m-1,n-1], p )
  c%coef(1:m,1:n,1:p) = a(1:m,1:n,1:p)

end subroutine array2chebyshev_series2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!              LEAST SQUARE FITTING ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
       call solve_NxN( c%coef(:,j), matmul( Ft, F ), matmul( Ft, y(:,j) ) )
    end do

  end subroutine chebfit1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !subroutine fcht1( &
  !     f, &
  !     c, &
  !     n, &
  !     dim, &
  !     trunc )
  !  ! Performs Fast Cosine Transform 
  !  implicit none
  !  integer,                      intent(in)  :: n, dim
  !  real(kind=MATHpr),            intent(in)  :: f(n,dim)
  !  type(type_chebyshev_series1), intent(out) :: c
  !  logical, intent(in), optional             :: trunc
  !  real(kind=MATHpr)                         :: w(3*n+15)
  !  integer                                   :: j
  !  call reset_chebyshev_series1( c, n-1, dim )
  !  c%coef = f / real( n-1, kind=MATHpr )
  !  call dcosti( n, w )
  !  do j = 1,dim
  !     call dcost( n, c%coef(:,j), w )
  !  end do
  !  c%coef([1,n],:) = 0.5_MATHpr * c%coef([1,n],:)
  !  if ( present(trunc) ) then
  !     if ( trunc ) then
  !        where ( abs(c%coef) < MATHeps ) c%coef = 0._MATHpr
  !     end if
  !  end if
  !end subroutine fcht1

  !subroutine fcht2( &
  !     f, &
  !     c, &
  !     trunc )
  !  implicit none
  !  real(kind=MATHpr),            intent(in)  :: f(:,:,:)
  !  type(type_chebyshev_series2), intent(out) :: c
  !  logical, intent(in), optional             :: trunc
  !  real(kind=MATHpr)                         :: w(3*max(size(f,1),size(f,2))+15)
  !  integer                                   :: degr(2), dim, i, j, k
  !  degr(1) = size(f,1) - 1
  !  degr(2) = size(f,2) - 1
  !  dim = size(f,3)
  !  call reset_chebyshev_series2( c, degr, dim )
  !  c%coef = f / real( degr(1), kind=MATHpr )
  !  call dcosti( degr(1)+1, w )
  !  do k = 1,dim
  !     do j = 1,degr(2)+1
  !        call dcost( degr(1)+1, c%coef(:,j,k), w )
  !     end do
  !  end do
  !  c%coef([1,degr(1)+1],:,:) = 0.5_MATHpr * c%coef([1,degr(1)+1],:,:)
  !  c%coef = c%coef / real( degr(2), kind=MATHpr )
  !  if ( degr(1) /= degr(2) ) call dcosti( degr(2)+1, w )
  !  do k = 1,dim
  !     do i = 1,degr(1)+1
  !        call dcost( degr(2)+1, c%coef(i,:,k), w )
  !     end do
  !  end do
  !  c%coef(:,[1,degr(2)+1],:) = 0.5_MATHpr * c%coef(:,[1,degr(2)+1],:)
  !  if ( present(trunc) ) then
  !     if ( trunc ) then
  !        where ( abs(c%coef) < MATHeps ) c%coef = 0._MATHpr
  !     end if
  !  end if
  !end subroutine fcht2


  !subroutine ifcht1( &
  !     c, &
  !     f )
  !  implicit none
  !  type(type_chebyshev_series1), intent(in)  :: c
  !  real(kind=MATHpr),            intent(out) :: f(c%degr+1,c%dim)
  !  real(kind=MATHpr)                         :: w(3*(c%degr+1)+15)
  !  integer                                   :: j
  !  f = c%coef
  !  f(2:c%degr,:) = 0.5_MATHpr * f(2:c%degr,:)
  !  call dcosti( c%degr+1, w )
  !  do j = 1,c%dim
  !     call dcost( c%degr+1, f(:,j), w )
  !  end do
  !end subroutine ifcht1


  !subroutine fcht2b( &
  !     f, &
  !     c, &
  !     m, &
  !     n, &
  !     dim, &
  !     trunc )
  !  implicit none
  !  integer,                      intent(in)           :: m, n, dim
  !  real(kind=MATHpr),            intent(in)           :: f(m,n,dim)
  !  type(type_chebyshev_series2), intent(out)          :: c
  !  logical,                      intent(in), optional :: trunc
  !  real(kind=MATHpr), allocatable                     :: wu(:), wv(:)
  !  integer                                       :: degr(2), i, j, k, degrprev(2)
  !  save wu, wv, degrprev
  !  data degrprev /0, 0/
  !  degr(1) = m - 1
  !  degr(2) = n - 1
  !  call reset_chebyshev_series2( c, degr, dim )
  !  c%coef = f / real( degr(1), kind=MATHpr )
  !  if ( allocated(wu) ) then
  !     if ( size(wu) < 3*m+15 ) deallocate( wu )
  !  end if
  !  if ( .not.allocated(wu) ) allocate( wu(3*m+15) )
  !  if ( degrprev(1) /= degr(1) ) call dcosti( m, wu )
  !  degrprev(1) = degr(1)
  !  do k = 1,dim
  !     do j = 1,n
  !        call dcost( m, c%coef(:,j,k), wu )
  !     end do
  !  end do
  !  c%coef([1,m],:,:) = 0.5_MATHpr * c%coef([1,m],:,:)
  !  c%coef = c%coef / real( degr(2), kind=MATHpr )
  !  if ( allocated(wv) ) then
  !     if ( size(wv) < 3*n+15 ) deallocate( wv )
  !  end if
  !  if ( .not.allocated(wv) ) allocate( wv(3*n+15) )
  !  if ( m == n ) then
  !     wv(1:3*n+15) = wu(1:3*n+15) 
  !  elseif ( degrprev(2) /= degr(2) ) then
  !     call dcosti( n, wv )
  !  end if
  !  degrprev(2) = degr(2)
  !  do k = 1,dim
  !     do i = 1,m
  !        call dcost( n, c%coef(i,:,k), wv )
  !     end do
  !  end do
  !  c%coef(:,[1,n],:) = 0.5_MATHpr * c%coef(:,[1,n],:)
  !  if ( present(trunc) ) then
  !     if ( trunc ) then
  !        where ( abs(c%coef) < MATHeps ) c%coef = 0._MATHpr
  !     end if
  !  end if
  !end subroutine fcht2b

  !subroutine clenshaw_alt( &
  !     f, &
  !     c, &
  !     x, &
  !     nx, &
  !     degr, &
  !     dim )
  !  implicit none
  !  integer,                      intent(in)  :: nx, degr, dim
  !  real(kind=MATHpr),            intent(in)  :: x(nx)
  !  real(kind=MATHpr),            intent(in)  :: c(degr+1,dim)
  !  real(kind=MATHpr),            intent(out) :: f(nx,dim)
  !  real(kind=MATHpr), dimension(nx)          :: c0, c1, tmp, twox
  !  integer                                   :: i, j
  !  if ( degr < 0 ) return
  !  do j = 1,dim
  !     if ( degr == 0 ) then
  !        c0(:) = c(1,j)
  !        c1(:) = 0._MATHpr
  !     elseif ( degr == 1 ) then
  !        c0(:) = c(1,j)
  !        c1(:) = c(2,j)
  !     else
  !        twox = 2._MATHpr * x
  !        c0(:) = c(degr,j)
  !        c1(:) = c(degr+1,j)
  !        do i = degr-1,1,-1
  !           tmp = c0
  !           c0 = c(i,j) - c1
  !           c1 = tmp + c1  * twox
  !        end do
  !     end if
  !     f(:,j) = c0 + c1 * x
  !  end do
  !end subroutine clenshaw_alt


  !subroutine chebdiff_col( &
  !     c, &
  !     d, &
  !     degr ) 
  !  implicit none
  !  integer,           intent(in)  :: degr
  !  real(kind=MATHpr), intent(in)  :: c(degr+1)
  !  real(kind=MATHpr), intent(out) :: d(degr)
  !  integer                        :: i
  !  d(:) = 0._MATHpr
  !  if ( degr < 1 ) return
  !  d(degr) = real( 2*degr, kind=MATHpr ) * c(degr+1)
  !  if ( degr > 1 ) then
  !     d(degr-1) = real( 2*(degr-1), kind=MATHpr ) * c(degr)
  !     if ( degr > 2 ) then
  !        do i = degr-2,1,-1
  !           d(i) = real( 2*i, kind=MATHpr ) * c(i+1) + d(i+2)
  !        end do
  !     end if
  !  end if
  !  d(1) = 0.5_MATHpr * d(1)
  !end subroutine chebdiff_col
end module mod_chebyshev_bak
