program complexite_chebyshev

  use mod_util
  use mod_math
  use mod_chebyshev
  use mod_polynomial

  implicit none

  !integer, parameter         :: N0 = 10, Nm = 120, m = 12
  integer, parameter         :: N0 = 2, Nm = 512, m = 9
  integer, parameter         :: nrep = 5000

  real(kind=fp), allocatable :: xcgl(:), f0cgl(:), f1cgl(:), DN(:,:)
  type(type_polynomial)      :: INf, DNf

  real*8                     :: timer(6)
  real*8                     :: stride, ratio
  integer*8                  :: count_rate, start, finish
  real*8                     :: tstart, tfinish
  integer                    :: N, i, r, fid

  stride = (Nm - N0)/dble(m-1)
  ratio = dble(Nm/N0)**(1.d0/dble(m-1))

  call get_free_unit(fid)
  open( unit=fid, &
       file='chebyshev/complexite.dat', &
       action='write' )
  do i = 1,m
     !N = nint(dble(N0) + dble(i-1)*stride)
     N = nint(dble(N0) * ratio**(i-1))
     print *,'N =',N
     allocate(xcgl(N+1), f0cgl(N+1), f1cgl(N+1), DN(N+1,N+1))

     call cgl_nodes(xcgl, N)

     call fun1(xcgl, N+1, f0cgl)
     call reset_polynomial( &
             INf, &
             1, &
             1, &
             [N], &
             1 )
     call fcht1( &
             f0cgl(1:N+1), &
             INf%coef(1:N+1,1,1), &
             N+1, &
             1 )

     ! construction matrice differentiation
     call system_clock(start, count_rate)
     call cpu_time(tstart)
     do r = 1,nrep
        call cheb_diffmatrix_cgl(N, DN)
     end do
     call cpu_time(tfinish)
     call system_clock(finish)
     !timer(1) = dble(finish - start)/dble(nrep*count_rate)
     timer(1) = (tfinish - tstart)/dble(nrep)
     
     ! produit matrice-vecteur
     call system_clock(start, count_rate)
     call cpu_time(tstart)
     do r = 1,nrep
        f1cgl = matmul(DN, f0cgl)
     end do
     call cpu_time(tfinish)
     call system_clock(finish)
     !timer(2) = dble(finish - start)/dble(nrep*count_rate)
     timer(2) = (tfinish - tstart)/dble(nrep)
     
     ! produit matrice-vecteur
     call system_clock(start, count_rate)
     call cpu_time(tstart)
     do r = 1,nrep
        call matmuldot(DN, f0cgl, f1cgl, N+1, N+1)
     end do
     call cpu_time(tfinish)
     call system_clock(finish)
     !timer(5) = dble(finish - start)/dble(nrep*count_rate)
     timer(5) = (tfinish - tstart)/dble(nrep)
     
     ! produit matrice-vecteur
     call system_clock(start, count_rate)
     call cpu_time(tstart)
     do r = 1,nrep
        call matmulloop(DN, f0cgl, f1cgl, N+1, N+1)
     end do
     call cpu_time(tfinish)
     call system_clock(finish)
     !timer(6) = dble(finish - start)/dble(nrep*count_rate)
     timer(6) = (tfinish - tstart)/dble(nrep)

     
     ! 'chebdiff'
     call system_clock(start, count_rate)
     call cpu_time(tstart)
     do r = 1,nrep
        call diff1(INf, DNf)
     end do
     call cpu_time(tfinish)
     call system_clock(finish)
     !timer(3) = dble(finish - start)/dble(nrep*count_rate)
     timer(3) = (tfinish - tstart)/dble(nrep)

     ! ifcht
     call system_clock(start, count_rate)
     call cpu_time(tstart)
     do r = 1,nrep
        call ifcht1( &
             DNf%coef(1:N,1:1,1), &
             f1cgl(1:N+1), &
             N+1, &
             1 )
     end do
     call cpu_time(tfinish)
     call system_clock(finish)
     !timer(4) = dble(finish - start)/dble(nrep*count_rate)
     timer(4) = (tfinish - tstart)/dble(nrep)

     write (fid,*) N, timer, timer(1)+timer(2), timer(3)+timer(4)

     deallocate(xcgl, f0cgl, f1cgl, DN)
  end do
  close(fid)


contains

  subroutine fun1(x, n, f)
    implicit none
    integer, intent(in)  :: n
    real*8,  intent(in)  :: x(n)
    real*8,  intent(out) :: f(n)

    f = 1._fp / ( 1._fp + 25._fp*x**2 )
    
  end subroutine fun1


  subroutine matmuldot(A, x, b, m, n)
    implicit none
    integer, intent(in)  :: m, n
    real*8,  intent(in)  :: A(m,n), x(n)
    real*8,  intent(out) :: b(m)
    integer              :: i

    do i = 1,m
       b(i) = dot_product(A(i,1:n), x(1:n))
    end do
    
  end subroutine matmuldot


  subroutine matmulloop(A, x, b, m, n)
    implicit none
    integer, intent(in)  :: m, n
    real*8,  intent(in)  :: A(m,n), x(n)
    real*8,  intent(out) :: b(m)
    integer              :: i, j

    do i = 1,m
       b(i) = 0.d0
       do j = 1,n
          b(i) = b(i) + A(i,j)*x(j)
       end do
    end do
    
  end subroutine matmulloop
  
end program complexite_chebyshev
