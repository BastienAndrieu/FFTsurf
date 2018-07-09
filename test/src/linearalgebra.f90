program linearalgebra

  use mod_math
  use mod_util
  use mod_linalg

  implicit none

  type string
     character(30) :: str = ''
  end type string


  LOGICAL :: VERBOSE = .FALSE.
  LOGICAL :: TEST_LINSOLVE = .false.
  LOGICAL :: TEST_SCALING = .true.
  
  integer, parameter :: repeat = 100
  integer, parameter :: pow_min = 1
  integer, parameter :: pow_max = 7
  real,    parameter :: pow_step = 0.5!1.0
  integer, parameter :: len = 1 + nint(real(pow_max-pow_min)/pow_step)

  integer                    :: fid, m, n, i, j, rank, in
  real(kind=fp), allocatable :: a(:,:), u(:,:), v(:,:), w(:), q(:,:), r(:,:), p(:,:), b(:), x(:)
  real(kind=fp)              :: cond
  integer*8                  :: tic, toc, count_rate
  integer                    :: listn(len)
  real                       :: timer(len,5)
  type(string)               :: label(5)
  character(200)             :: header


  if ( TEST_LINSOLVE ) then
  call get_free_unit( fid )
  open( unit=fid, file='linearalgebra/A.dat', action='read' )
  read (fid,*) m, n
  allocate( a(m,n), u(m,n), v(n,n), w(n) )
  do j = 1,n
     do i = 1,m
        read (fid,*) a(i,j)
     end do
  end do
  close(fid)

  IF (VERBOSE) THEN
     print *,'A ='
     call print_mat( a )
     print *,''
  END IF

  u = a
  !call system_clock( tic, count_rate )
  !call svdcmp( u, m, n, w, v )
  !call system_clock( toc )
  !print *,'elapsed :', real(toc - tic) / real(count_rate)

  IF (VERBOSE) THEN
     print *,'U ='
     call print_mat( a )
     print *,''

     print *,'W ='
     call print_mat( diag(w) )
     print *,''

     print *,'V ='
     call print_mat( v )
     print *,''
  END IF

  !print *,'max( |A - U*W*Vt| ) ='
  !call print_mat( a - matmul( matmul( u, diag(w) ), transpose(v) ) )
  !print *,maxval( abs( a - matmul( matmul( u, diag(w) ), transpose(v) ) ) )
  !print *,''

  !print *,'cond =',maxval(w) / minval(w)



  allocate( q(m,m), r(m,n), p(n,n) )
  !call system_clock( tic, count_rate )
  !call qr_colpiv( a, q, r, p, rank, m, n )
  !call system_clock( toc )
  !print *,'elapsed :', real(toc - tic) / real(count_rate)



  ! linear solve
  allocate( b(m), x(n) )
  call get_free_unit( fid )
  open( unit=fid, file='linearalgebra/b.dat', action='read' )
  do i = 1,m
     read (fid,*) b(i)
  end do
  close(fid)

  u = a
  print *,''
  print *,'linsolve_svd'
call system_clock( tic, count_rate )
  call linsolve_svd( &
       x, &
       u, &
       b, &
       m, &
       n, &
       1, &
       cond, &
       rank )
  call system_clock( toc )
  print *,'elapsed :', real(toc - tic) / real(count_rate)

  print *,'rank = ',rank,', cond =',cond
  

  print *,'x = ',x
  print *,'err =', maxval( abs( matmul(A,x) - b ) )

  print *,''
  print *,'linsolve_qr'
  call system_clock( tic, count_rate )
  call linsolve_qr( &
       x, &
       A, &
       b, &
       m, &
       n, &
       1, &
       rank )
  call system_clock( toc )
  print *,'elapsed :', real(toc - tic) / real(count_rate)
  print *,'rank = ',rank
  print *,'x = ',x
  print *,'err =', maxval( abs( matmul(A,x) - b ) )

end if



if ( TEST_SCALING ) then

   if ( allocated( a ) ) deallocate( a )
   if ( allocated( u ) ) deallocate( u )
   if ( allocated( v ) ) deallocate( v )
   if ( allocated( w ) ) deallocate( w )
   if ( allocated( q ) ) deallocate( q )
   if ( allocated( r ) ) deallocate( r )
   if ( allocated( p ) ) deallocate( p )
   if ( allocated( b ) ) deallocate( b )
   if ( allocated( x ) ) deallocate( x )


   do in = 1,len
      n = nint( 2.0**( real(pow_min) + real(iN-1)*pow_step ) )
      m = n
      listn(in) = n
      PRINT *,'N =',N

      allocate( a(m,n), u(m,n), v(n,n), w(n), q(m,m), r(m,n), p(n,n), b(n), x(n) )

      call random_number( a )
      a = 2._fp * a - 1._fp
      call random_number( b )
      b = 2._fp * b - 1._fp

      call system_clock( tic, count_rate )
      do i = 1,repeat
         u = a
      end do
      call system_clock( toc )
      timer(in,1) = real( toc - tic ) / real( count_rate )
      label(1)%str = 'copy'

      call system_clock( tic, count_rate )
      do i = 1,repeat
         u = a
         call svdcmp( u, m, n, w, v )
      end do
      call system_clock( toc )
      timer(in,2) = real( toc - tic ) / real( count_rate )
      label(2)%str = 'SVD'
      
      call system_clock( tic, count_rate )
      do i = 1,repeat
         u = a
         call linsolve_svd( &
              x, &
              u, &
              b, &
              m, &
              n, &
              1 )
      end do
      call system_clock( toc )
      timer(in,3) = real( toc - tic ) / real( count_rate )
      label(3)%str = 'linsolve SVD'


      
      call system_clock( tic, count_rate )
      do i = 1,repeat
         call qr_colpiv( &
              a, &
              q, &
              r, &
              p, &
              rank, &
              m, &
              n )
      end do
      call system_clock( toc )
      timer(in,4) = real( toc - tic ) / real( count_rate )
      label(4)%str = 'QR'

      
      call system_clock( tic, count_rate )
      do i = 1,repeat
         call linsolve_qr( &
              x, &
              A, &
              b, &
              m, &
              n, &
              1 )
      end do
      call system_clock( toc )
      timer(in,5) = real( toc - tic ) / real( count_rate )
      label(5)%str = 'linsolve QR'

      deallocate( a, u, v, w, q, r, p, b, x )

   end do

   timer = timer / real(repeat)

   call get_free_unit( fid )
   open( unit=fid, file='linearalgebra/complexite.csv', action='write' )
   header = 'n'
   
   !write ( fid, fmt='(A1,A1,1x)', advance='no' ) 'n', ','
   do i = 1,size(label)
      header = trim(header)// ', ' // trim(label(i)%str)
      !write ( fid, fmt='(A0)', advance='no' ) trim(label(i)%str)
      !if ( i < size(label) ) then
      !   write ( fid, fmt='(A1,1x)' ) ','
      !else
      !   write ( fid, * )
      !end if
   end do
   write ( fid, * ) trim(header)
   do in = 1,len
      write ( fid, fmt='(I0,A1,1x)', advance='no' ) listn(in), ','
      do i = 1,size(timer,2)
         write ( fid, fmt='(E13.6)', advance='no' ) timer(in,i)
         if ( i < size(timer,2) ) then
            write ( fid, fmt='(A1,1x)', advance='no' ) ','
         else
            write ( fid, * )
         end if
      end do
   end do
   close(fid)

end if










end program linearalgebra
