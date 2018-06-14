program perfo_subdivision

  use mod_util
  use mod_math
  use mod_chebyshev
  use mod_bernstein

  implicit none

  integer, parameter             :: pow_min = 2
  integer, parameter             :: pow_max = 6
  real,    parameter             :: pow_step = 0.5!1.0
  integer, parameter             :: lenN = 1 + nint(real(pow_max-pow_min)/pow_step)

  integer, parameter             :: nrep1 = 10000
  integer, parameter             :: nrep2 = 1000

  type(type_chebyshev_series1)   :: c1, c1s
  type(type_chebyshev_series2)   :: c2, c2s
  type(type_bernstein_series1)   :: b1, bl, br
  type(type_bernstein_series2)   :: b2, bsw, bse, bnw, bne

  real(kind=MATHpr)              :: t, tc, uv(2), uvc(2)

  integer                        :: N
  real(kind=MATHpr), allocatable :: c2bu(:,:), c2bv(:,:)

  real*8                         :: timer(2)

  integer                        :: file_unit
  integer*8                      :: tic, toc, count_rate
  integer                        :: rep, iN, k

  !! CURVE
  PRINT *,'CURVE'
  call random_number( t )
  tc = 2._MATHpr * t - 1._MATHpr

  call get_free_unit( file_unit )
  open( unit=file_unit, file='perfo_subdivision/timer_curve.dat', action='write' )
  do iN = 1,lenN
     N = nint( 2.0**( real(pow_min) + real(iN-1)*pow_step ) )
     PRINT *,'N =',N

     call random_curve( b1, c1, N-1 )

     call system_clock( tic, count_rate )
     do rep = 1,nrep1
        call subdiv1( &
             b1, &
             t, &
             bl, &
             br )
     end do
     call system_clock( toc )
     timer(1) = dble( toc - tic ) / dble( count_rate * nrep1 )
     
     

     allocate( c2bu(N,N) )
     call ch2be_matrix( c2bu, N-1 )

     call system_clock( tic, count_rate )
     do rep = 1,nrep1
        call reset_bernstein_series1( bl, c1%degr, c1%dim )
        call chgvar1( &
             c1, &
             c1s, &
             -1._MATHpr, &
             tc )
        bl%coef(1:c1%degr+1,1:c1%dim) = matmul( c2bu, c1s%coef(1:c1%degr+1,1:c1%dim) )
        
        call reset_bernstein_series1( br, c1%degr, c1%dim )
        call chgvar1( &
             c1, &
             c1s, &
             tc, &
             1._MATHpr )
        br%coef(1:c1%degr+1,1:c1%dim) = matmul( c2bu, c1s%coef(1:c1%degr+1,1:c1%dim) )
     end do
     call system_clock( toc )
     timer(2) = dble( toc - tic ) / dble( count_rate * nrep1 )
     deallocate( c2bu )
     
     write ( file_unit, * ) N, timer

  end do
  close( file_unit )




  !! SURFACE
  PRINT *,'SURFACE'
  call random_number( uv )
  uvc = 2._MATHpr * uv - 1._MATHpr

   call get_free_unit( file_unit )
   open( unit=file_unit, file='perfo_subdivision/timer_surface.dat', action='write' )
   do iN = 1,lenN
      N = nint( 2.0**( real(pow_min) + real(iN-1)*pow_step ) )
      PRINT *,'N =',N

      call random_surface( b2, c2, N-1 )

      call system_clock( tic, count_rate )
      do rep = 1,nrep2
         call subdiv2( &
              b2, &
              uv, &
              bsw, &
              bse, &
              bnw, &
              bne )
      end do
      call system_clock( toc )
      timer(1) = dble( toc - tic ) / dble( count_rate * nrep2 )


      allocate( c2bu(N,N), c2bv(N,N) )
      call ch2be_matrix( c2bu, N-1 )
      c2bv = transpose( c2bu )
      
      call system_clock( tic, count_rate )
      do rep = 1,nrep2
         call reset_bernstein_series2( bsw, c2%degr, c2%dim )
         call chgvar2( &
              c2, &
              c2s, &
              [-1._MATHpr, -1._MATHpr], &
              uvc )
         do k = 1,c2%dim
            bsw%coef(:,:,k) = matmul( matmul( c2bu, c2s%coef(:,:,k) ), c2bv )
         end do
         

         call reset_bernstein_series2( bse, c2%degr, c2%dim )
         call chgvar2( &
              c2, &
              c2s, &
              [uvc(1), -1._MATHpr], &
              [1._MATHpr, uvc(2)] )
         do k = 1,c2%dim
            bse%coef(:,:,k) = matmul( matmul( c2bu, c2s%coef(:,:,k) ), c2bv )
         end do

         
         call reset_bernstein_series2( bnw, c2%degr, c2%dim )
         call chgvar2( &
              c2, &
              c2s, &
              [-1._MATHpr, uvc(2)], &
              [uvc(1), 1._MATHpr] )
         do k = 1,c2%dim
            bnw%coef(:,:,k) = matmul( matmul( c2bu, c2s%coef(:,:,k) ), c2bv )
         end do

         
         call reset_bernstein_series2( bne, c2%degr, c2%dim )
         call chgvar2( &
              c2, &
              c2s, &
              uvc, &
              [1._MATHpr, 1._MATHpr] )
         do k = 1,c2%dim
            bne%coef(:,:,k) = matmul( matmul( c2bu, c2s%coef(:,:,k) ), c2bv )
         end do

      end do
      call system_clock( toc )
      timer(2) = dble( toc - tic ) / dble( count_rate * nrep2 )
      deallocate( c2bu, c2bv )

      write ( file_unit, * ) N, timer

   end do
   close( file_unit )




contains
  
  
  subroutine random_curve( b, c, degr )
    implicit none
    integer,                      intent(in)  :: degr
    type(type_bernstein_series1), intent(out) :: b
    type(type_chebyshev_series1), intent(out) :: c
    real(kind=MATHpr)                         :: c2b(degr+1,degr+1)
    integer                                   :: i
    
    call reset_bernstein_series1( b, degr, 3 )
    call reset_chebyshev_series1( c, degr, 3 )

    call random_number( c%coef )
    c%coef = c%coef / spread( real( [(i**2, i=1,degr+1)], kind=MATHpr ), dim=2, ncopies=c%dim )


    call ch2be_matrix( c2b, degr )
    b%coef = matmul( c2b, c%coef )

  end subroutine random_curve
  
  

  
  subroutine random_surface( b, c, degr )
    implicit none
    integer,                      intent(in)  :: degr
    type(type_bernstein_series2), intent(out) :: b
    type(type_chebyshev_series2), intent(out) :: c
    real(kind=MATHpr)                         :: c2b(degr+1,degr+1)
    integer                                   :: i

    call reset_bernstein_series2( b, [degr, degr], 3 )
    call reset_chebyshev_series2( c, [degr, degr], 3 )

    call random_number( c%coef )
    do k = 1,c%dim
       c%coef(:,:,k) = c%coef(:,:,k) / ( &
            spread( real( [(i**2, i=1,degr+1)], kind=MATHpr ), dim=2, ncopies=degr+1 ) * &
            spread( real( [(i**2, i=1,degr+1)], kind=MATHpr ), dim=1, ncopies=degr+1 ) )
    end do


    call ch2be_matrix( c2b, degr )
    do k = 1,c%dim
       b%coef(:,:,k) = matmul( matmul( c2b, c%coef(:,:,k) ), transpose(c2b) )
    end do

  end subroutine random_surface


end program perfo_subdivision
