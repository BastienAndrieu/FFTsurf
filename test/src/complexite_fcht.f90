program complexite_fcht
  
  use mod_util
  use mod_math
  use mod_chebyshev

  implicit none

  integer, parameter :: pow_min = 2
  integer, parameter :: pow_max = 7
  real,    parameter :: pow_step = 0.5!1.0
  integer, parameter :: lenN = 1 + nint(real(pow_max-pow_min)/pow_step)
  
  !integer, parameter :: nrepeat1 = 10000
  integer, parameter :: nrepeat2 = 1000
  !integer, parameter :: dim = 3

  !type(type_chebyshev_series1) :: c1
  type(type_chebyshev_series2) :: c2
  integer                      :: N, iN

  !real(kind=MATHpr), allocatable :: x(:), f1(:,:)
  real(kind=MATHpr), allocatable :: f2a(:,:,:), f2b(:,:,:)
  real(kind=MATHpr), allocatable :: w(:)

  integer*8 :: tic, toc, count_rate
  integer   :: rep, file_unit
  real      :: data(3,lenN)


  
  !! UNIVARIATE
  !do iN = 1,lenN
  !   N = nint( 2.0**( real(pow_min) + real(iN-1)*pow_step ) )
  !   data(1,iN) = N
  !   allocate( x(N), f1(N,dim) )
  !   call cgl_points( &
  !        x, &
  !        1._MATHpr, &
  !        -1._MATHpr, &
  !        N )  
  !   f1 = fun1( x, dim )
  !   call system_clock( tic, count_rate )
  !   call system_clock( toc )
  !   data(2,iN) = real( toc - tic ) / real( count_rate * nrepeat1 )   
  !end do


  !! BIVARIATE
  do iN = 1,lenN

     N = nint( 2.0**( real(pow_min) + real(iN-1)*pow_step ) )
     print *,'N =',N
     data(1,iN) = N
     
     call random_surface( N, c2 )

     allocate( f2a(N,N,c2%dim), f2b(N,N,c2%dim) )


     ! IFCHT2 subroutine
     call system_clock( tic, count_rate )
     do rep = 1,nrepeat2
        call ifcht2( c2, f2a )
     end do
     call system_clock( toc )
     data(2,iN) = real( toc - tic ) / real( count_rate * nrepeat2 )   


     ! keep 'global' wsave
     allocate( w(3*N+15) )
     call dcosti( N, w )

     call system_clock( tic, count_rate )
     do rep = 1,nrepeat2
        call ifcht2_globalw( c2, f2b, w, w )
     end do
     call system_clock( toc )
     data(3,iN) = real( toc - tic ) / real( count_rate * nrepeat2 )   

     deallocate( w )



     ! check similarity
     print *,'difference =', maxval( abs(f2a - f2b) )
     


     deallocate( f2a, f2b )

  end do


  call get_free_unit( file_unit )
  open( unit=file_unit, file='complexite_fcht/2var.dat', action='write' )
  do iN = 1,lenN
     write (file_unit,*) data(:,iN)
  end do
  close(file_unit)




contains

  
  function fun1( x, dim ) result( f )
    implicit none
    integer           :: dim
    real(kind=MATHpr) :: x(:), f(size(x),dim)
    integer           :: j
    
    f(:,1) = exp( x ) * sin( 5._MATHpr * x )
    do j = 2,dim
       f(:,j) = f(:,1) * cos( real(j-1,kind=MATHpr) * x )
    end do

  end function fun1
  
  
  function fun2( x, y, dim ) result( f )
    implicit none

    integer :: dim
    real(kind=MATHpr) :: x(:), y(:), f(size(x),size(y),dim)
    real(kind=MATHpr), dimension(size(x),size(y)) :: xx, yy
    integer :: k

    xx = spread( x, 2, size(y) )
    yy = spread( y, 1, size(x) )

    f(:,:,1) = exp( xx ) * sin( 5._MATHpr * yy )
    do k = 2,dim
       f(:,:,k) = f(:,:,1) * cos( real(k-1,kind=MATHpr) * xx ) * sin( real(k-1,kind=MATHpr) * yy )
    end do

  end function fun2





  subroutine random_surface( N, c )
    implicit none
    real(kind=MATHpr), parameter :: noise = 0.2_MATHpr
    integer, intent(in) :: N
    type(type_chebyshev_series2), intent(out) :: c
    real(kind=MATHpr) :: bcp(N,N,3)
    real(kind=MATHpr) :: c2b(N,N)
    integer           :: i

    call random_number( bcp )
    bcp = noise * ( 2._MATHpr * bcp - 1._MATHpr ) 
    
    bcp(:,:,1) = bcp(:,:,1) + spread( linspace( -1._MATHpr, 1._MATHpr, N ), dim=2, ncopies=N )
    bcp(:,:,2) = bcp(:,:,2) + spread( linspace( -1._MATHpr, 1._MATHpr, N ), dim=1, ncopies=N )    
    
    call ch2be_matrix( c2b, N-1 )

    call reset_chebyshev_series2( c, [N-1,N-1], 3 )
    do i = 1,3
       c%coef(:,:,i) = matmul( c2b, matmul( bcp(:,:,i), transpose( c2b ) ) )
    end do

  end subroutine random_surface










  subroutine ifcht2_globalw( &
       c, &
       f, &
       wu, &
       wv )
    implicit none
    type(type_chebyshev_series2), intent(in)  :: c
    real(kind=MATHpr),            intent(in)  :: wu(:), wv(:)
    real(kind=MATHpr),            intent(out) :: f(c%degr(1)+1,c%degr(2)+1,c%dim)
    integer                                   :: i, j, k

    f = c%coef(1:c%degr(1)+1,1:c%degr(2)+1,:)

    f(2:c%degr(1),:,:) = 0.5_MATHpr * f(2:c%degr(1),:,:)
    do k = 1,c%dim
       do j = 1,c%degr(2)+1
          call dcost( c%degr(1)+1, f(:,j,k), wu )
       end do
    end do

    f(:,2:c%degr(2),:) = 0.5_MATHpr * f(:,2:c%degr(2),:)
    do k = 1,c%dim
       do i = 1,c%degr(1)+1
          call dcost( c%degr(2)+1, f(i,:,k), wv )
       end do
    end do

  end subroutine ifcht2_globalw
  

end program complexite_fcht
