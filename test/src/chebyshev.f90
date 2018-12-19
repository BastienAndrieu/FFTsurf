program chebyshev

  use mod_util
  use mod_math
  use mod_chebyshev
  use mod_polynomial

  implicit none

  real*8, parameter          :: eps = 1.d1 * epsilon(1.d0)
  integer, parameter         :: N0 = 10, Nm = 80, m = 8
  integer, parameter         :: n_verif = 10000

  real(kind=fp)              :: x(n_verif), f(n_verif,3), y(n_verif,3), z(n_verif,2)
  real(kind=fp), allocatable :: xcgl(:), fcgl(:,:)
  type(type_polynomial)      :: INf(3), DNf(2)
  real*8                     :: stride
  integer                    :: N, i, j, fid, fid2, fid3
  
  stride = (Nm - N0)/dble(m-1)
  x = linspace(-1._fp, 1._fp, n_verif)
  call fun1(x, n_verif, f)

  call get_free_unit(fid)
  open(unit=fid, file='../../These/memoire/figures/data/xf012.dat', action='write')
  do i = 1,n_verif
     write (fid,*) x(i), f(i,1), f(i,2), f(i,3)
  end do
  close(fid)
  
  open(unit=fid, file='../../These/memoire/figures/data/convergence_INf012.dat', action='write')

  call get_free_unit(fid2)
  open(unit=fid2, file='../../These/memoire/figures/data/convergence_DNf012.dat', action='write')
  do i = 1,m
     N = nint(dble(N0) + dble(i-1)*stride)

     allocate(xcgl(N+1), fcgl(N+1,3))
     call cgl_nodes(xcgl, N)

     call fun1(xcgl, N+1, fcgl)

     do j = 1,3 ! <------------------------------+
        call reset_polynomial( &                 !
             INf(j), &                           !
             1, &                                !
             1, &                                !
             [N], &                              !
             1 )                                 !
        !                                        !
        call fcht1( &                            !
             fcgl(1:N+1,j), &                    !
             INf(j)%coef(1:N+1,1,1), &           !
             N+1, &                              !
             1)                                  !
        !                                        !
        call polyval1( &                         !
             y(1:n_verif,j), &                   !
             INf(j), &                           !
             x(1:n_verif), &                     !
             n_verif )                           !
        !                                        !
        if ( j > 1 ) then ! <------------------+ !
           if ( j == 2 ) then ! <------------+ ! !
              call diff1(INf(1), DNf(1))     ! ! !
           else ! ---------------------------+ ! !
              call diff1(DNf(j-2), DNf(j-1)) ! ! !
           end if ! <------------------------+ ! !
           call polyval1( &                    ! !
                z(1:n_verif,j-1), &            ! !
                DNf(j-1), &                    ! !
                x(1:n_verif), &                ! !
                n_verif )                      ! !
        end if ! <-----------------------------+ !
     end do ! <----------------------------------+
     write (fid,'(i0,1x,e22.15,1x,e22.15,1x,e22.15)') N, maxval(abs(y - f), dim=1) / maxval(abs(f), dim=1)
     write (fid2,'(i0,1x,e22.15,1x,e22.15)') N, maxval(abs(z - f(:,2:3)), dim=1) / maxval(abs(f(:,2:3)), dim=1)

     if ( i == m ) then
        call get_free_unit(fid3)
        open(unit=fid3, file='../../These/memoire/figures/data/coeffs_f012.dat', action='write')
        do j = 1,N-1
           write (fid3,'(i0,1x,e22.15,1x,e22.15,1x,e22.15)') j-1, &
                abs(INf(1)%coef(j,1,1)), abs(DNf(1)%coef(j,1,1)), abs(DNf(2)%coef(j,1,1))
        end do
        write (fid3,'(i0,1x,e22.15,1x,e22.15,1x,e22.15)') N-1, &
             abs(INf(1)%coef(N,1,1)), abs(DNf(1)%coef(N,1,1)), 0._fp
        write (fid3,'(i0,1x,e22.15,1x,e22.15,1x,e22.15)') N, &
             abs(INf(1)%coef(N+1,1,1)), 0._fp, 0._fp
        close(fid3)
     end if

     deallocate(xcgl, fcgl)
  end do
  close(fid)
  close(fid2)
  
contains

  subroutine fun1(x, n, f)
    implicit none
    real*8, parameter              :: a = 3.d0, b = 3.d0
    integer, intent(in)            :: n
    real*8,  intent(in)            :: x(n)
    real*8,  intent(out)           :: f(n,3)
    real*8, dimension(n)           :: g, dg, dg2, h, dh, dh2

    h = a*x**b! + 1
    dh = b*a*x**(b-1)
    dh2 = (b-1)*b*a*x**(b-2)

    g = sin(h)!cos(h)
    dg = cos(h)!-sin(h)
    dg2 = -sin(h)!-cos(h)

    f(1:n,1) = exp(g)
    f(1:n,2) = f(1:n,1)*dh*dg
    f(1:n,3) = f(1:n,1)*(dh2*dg + dg2*dh**2 + (dh*dg)**2)
  end subroutine fun1

end program chebyshev
