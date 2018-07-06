module mod_linalg
  ! Linear Algebra routines

  use mod_math

  implicit none

contains


  subroutine linsolve_svd( &
       x, &
       A, &
       b, &
       m, &
       n, &
       p, &
       cond, &
       rank )
    implicit none
    integer,                 intent(in)    :: m, n, p
    real(kind=fp),           intent(inout) :: A(m,n)
    real(kind=fp),           intent(in)    :: b(m,p)
    real(kind=fp),           intent(out)   :: x(n,p)
    real(kind=fp), optional, intent(out)   :: cond
    integer,       optional, intent(out)   :: rank
    real(kind=fp)                          :: V(n,n), W(n), wmax, wmin
    integer                                :: j

    ! compute Singular Value Decomposition of A (A is replaced by U)
    call svdcmp( A, m, n, W, V )

    ! get maximal singular value
    wmax = maxval( w )
    if ( present(cond) ) cond = wmax / minval(W) ! condition number
    wmin = real( max(m,n), kind=fp ) * wmax * epsilon( 1._fp )
    if ( present(rank) ) rank = count( w >= wmin )
    where( w < wmin ) w = 0._fp

    ! back substitution
    do j = 1,p
       call svbksb( &
            A, &
            W, &
            V, &
            b(:,j), &
            x(:,j), &
            m, &
            n )
    end do

  end subroutine linsolve_svd






  subroutine svbksb( &
       u, &
       w, &
       v, &
       b, &
       x, &
       m, &
       n )
    ! from "Numerical recipes in Fortran 90" (1996)
    ! Solves Ax = b for a vector x, where A is specified by the arrays u, v, w as returned by svdcmp.
    ! Here u is M × N , v is N × N ,and w is of length N. b is the M-dimensional input right-hand side. 
    ! x is the N-dimensional output solution vector. No input quantities are destroyed, 
    ! so the routine may be called sequentially with different b’s.
    implicit none
    integer,       intent(in)  :: m, n
    real(kind=fp), intent(in)  :: u(m,n), w(n), v(n,n), b(m)
    real(kind=fp), intent(out) :: x(n)
    real(kind=fp)              :: tmp(n)

    where ( .not.is_zero(w) )
       tmp = matmul( b, u ) / w
    elsewhere
       tmp = 0._fp
    end where
    x = matmul( v, tmp )

  end subroutine svbksb



  subroutine svdcmp( &
       a, &
       m, &
       n, &
       w, &
       v )
    ! Singular Value Decomposition
    ! from "Numerical recipes in Fortran 90" (1996)
    ! Given an M × N matrix a , this routine computes its singular value decomposition, 
    ! A =U · W · Vt . The matrix U replaces a on output. The diagonal matrix of singular values
    ! W is output as the N -dimensional vector w . The N × N matrix V (not the transpose Vt )
    ! is output as v .
    implicit none
    integer,       intent(in)    :: m, n
    real(kind=fp), intent(inout) :: a(m,n)
    real(kind=fp), intent(out)   :: v(n,n), w(n)
    real(kind=fp)                :: tempm(m), rv1(n), tempn(n)
    real(kind=fp)                :: anorm, c, f, g, h, s, scale, x, y, z
    integer                      :: i, its, j, k, l, nm

    g = 0._fp
    scale = 0._fp
    do i = 1,n ! Householder reduction to bidiagonal form.
       l = i + 1
       rv1(i) = scale * g
       g = 0._fp
       scale = 0._fp
       if (i <= m) then
          scale = sum( abs(a(i:m,i)) )
          if ( .not.is_zero(scale) ) then !if (scale /= 0.0) then
             a(i:m,i) = a(i:m,i) / scale
             s = dot_product( a(i:m,i), a(i:m,i) )
             f = a(i,i)
             g = -sign( sqrt(s), f )
             h = f * g - s
             a(i,i) = f - g
             tempn(l:n) = matmul( a(i:m,i), a(i:m,l:n) ) / h
             a(i:m,l:n) = a(i:m,l:n) + outer_product( a(i:m,i), tempn(l:n) )
             a(i:m,i) = scale * a(i:m,i)
          end if
       end if
       w(i) = scale * g
       g = 0._fp
       scale = 0._fp
       if ( i <= m .and. i /= n ) then
          scale = sum( abs(a(i,l:n)) )
          if ( .not.is_zero(scale) ) then !if (scale /= 0.0) then
             a(i,l:n) = a(i,l:n) / scale
             s = dot_product( a(i,l:n), a(i,l:n) )
             f = a(i,l)
             g = -sign( sqrt(s), f )
             h = f * g - s
             a(i,l) = f - g
             rv1(l:n) = a(i,l:n) / h
             tempm(l:m) = matmul( a(l:m,l:n), a(i,l:n) )
             a(l:m,l:n) = a(l:m,l:n) + outer_product( tempm(l:m), rv1(l:n) )
             a(i,l:n) = scale * a(i,l:n)
          end if
       end if
    end do

    anorm = maxval( abs(w) + abs(rv1) )

    do i=n,1,-1 ! Accumulation of right-hand transformations.
       if (i < n) then
          if ( .not.is_zero(g) ) then !if (g /= 0.0) then
             v(l:n,i) = ( a(i,l:n) / a(i,l) ) / g
             tempn(l:n) = matmul( a(i,l:n), v(l:n,l:n) )
             v(l:n,l:n) = v(l:n,l:n) + outer_product( v(l:n,i), tempn(l:n) )
          end if
          v(i,l:n)=0._fp
          v(l:n,i)=0._fp
       end if
       v(i,i) = 1._fp
       g = rv1(i)
       l = i
    end do

    do i = min(m,n),1,-1 ! Accumulation of left-hand transformations.
       l = i + 1
       g = w(i)
       a(i,l:n) = 0._fp
       if ( .not.is_zero(g) ) then !if (g /= 0.0) then
          g = 1._fp / g
          tempn(l:n) = ( matmul( a(l:m,i), a(l:m,l:n) ) / a(i,i) ) * g
          a(i:m,l:n) = a(i:m,l:n) + outer_product( a(i:m,i), tempn(l:n) )
          a(i:m,i) = a(i:m,i) * g
       else
          a(i:m,i) = 0._fp
       end if
       a(i,i) = a(i,i) + 1._fp
    end do

    do k = n,1,-1 ! Diagonalization of the bidiagonal form: Loop over singular values, and over allowed iterations.
       do its = 1,30
          do l = k,1,-1 ! Test for splitting.
             nm = l - 1
             if ( abs(rv1(l)) + anorm <= anorm + epsilon(1._fp) ) exit !if ( abs(rv1(l)) + anorm == anorm) exit
             ! Note that rv1(1) is always zero, so can never fall through bottom of loop.
             if ( abs(w(nm)) + anorm <= anorm + epsilon(1._fp) ) then !if ( abs(w(nm)) + anorm == anorm) then
                c = 0._fp ! Cancellation of rv1(l), if l > 1.
                s = 1._fp
                do i = l,k
                   f = s * rv1(i)
                   rv1(i) = c * rv1(i)
                   if ( abs(f) + anorm < anorm + epsilon(1._fp) ) exit !if ((abs(f)+anorm) == anorm) exit
                   g = w(i)
                   h = hypot( f, g )
                   w(i) = h
                   h = 1._fp / h
                   c = ( g * h )
                   s = -( f * h )
                   tempm(1:m) = a(1:m,nm)
                   a(1:m,nm) = a(1:m,nm)*c + a(1:m,i)*s
                   a(1:m,i) = -tempm(1:m)*s + a(1:m,i)*c
                end do
                exit
             end if
          end do
          z = w(k)
          if ( l == k ) then ! Convergence.
             if ( z < 0._fp ) then ! Singular value is made nonnegative.
                w(k) = -z
                v(1:n,k) = -v(1:n,k)
             end if
             exit
          end if
          if (its == 30) STOP "svdcmp_dp: no convergence in svdcmp"
          x = w(l) ! Shift from bottom 2-by-2 minor.
          nm = k - 1
          y = w(nm)
          g = rv1(nm)
          h = rv1(k)
          f = ( (y-z)*(y+z) + (g-h)*(g+h) ) / ( 2._fp * h * y )
          g = hypot( f, 1._fp )
          f = ( (x-z)*(x+z) + h*((y/(f+sign(g,f)))-h) ) / x
          c = 1._fp ! Next QR transformation:
          s = 1._fp
          do j = l,nm
             i = j + 1
             g = rv1(i)
             y = w(i)
             h = s*g
             g = c*g
             z = hypot( f, h )
             rv1(j) = z
             c = f/z
             s = h/z
             f =  (x*c) + (g*s)
             g = -(x*s) + (g*c)
             h = y*s
             y = y*c
             tempn(1:n) = v(1:n,j)
             v(1:n,j) = v(1:n,j)*c + v(1:n,i)*s
             v(1:n,i) = -tempn(1:n)*s + v(1:n,i)*c
             z = hypot( f, h )
             w(j) = z ! Rotation can be arbitrary if z = 0.
             if ( .not.is_zero(z) ) then !if (z /= 0.0) then
                z = 1._fp / z
                c = f*z
                s = h*z
             end if
             f =  (c*g) + (s*y)
             x = -(s*g) + (c*y)
             tempm(1:m) = a(1:m,j)
             a(1:m,j) =  a(1:m,j)*c + a(1:m,i)*s
             a(1:m,i) = -tempm(1:m)*s + a(1:m,i)*c
          end do
          rv1(l) = 0._fp
          rv1(k) = f
          w(k) = x
       end do
    end do

  end subroutine svdcmp












  subroutine linsolve_qr2( &
       x, &
       A, &
       b, &
       m, &
       n, &
       p, &
       rank )
    implicit none
    integer,           intent(in)  :: m, n, p
    real(kind=fp),     intent(in)  :: A(m,n)
    real(kind=fp),     intent(in)  :: b(m,p)
    real(kind=fp),     intent(out) :: x(n,p)
    integer, optional, intent(out) :: rank
    real(kind=fp)                  :: Q(m,m), R(m,n), Perm(n,n), c(m,p), y(n,p)
    integer                        :: rnk, j

    if (m < n) STOP 'linsolve_qr : underdetermined systems not supported (m < n)'

    ! QR factorization with column pivoting
    call QR_colpiv( &
         A, &
         Q, &
         R, &
         Perm, &
         rnk )
    if ( present(rank) ) rank = rnk

    c = matmul( transpose(Q), b )

    if ( rnk == n ) then ! rnk == min(m,n)
       ! A has full column rank (
       do j = 1,p
          call solve_up_tri( &
               y(1:n,j), &
               R(1:n,1:n), &
               c(1:n,j), &
               n, &
               n )
       end do
    else
       ! A is rank-deficient
       call linsolve_qr_rankdef( &
            y, &
            R(1:rnk,1:n), &
            c(1:rnk,1:p), &
            n, &
            p, &
            rnk )
    end if

    x = matmul( Perm, y )

  end subroutine linsolve_qr2



  subroutine linsolve_qr_rankdef( &
       y, &
       R, &
       c, &
       n , &
       p, &
       rank )
    implicit none
    integer,       intent(in)  :: n, p, rank
    real(kind=fp), intent(in)  :: R(rank,n)
    real(kind=fp), intent(in)  :: c(rank,p)
    real(kind=fp), intent(out) :: y(n,p)
    real(kind=fp)              :: y1(rank)
    real(kind=fp)              :: y2(n-rank,p)
    real(kind=fp)              :: S(n,n-rank)
    real(kind=fp)              :: t(n,p)
    integer                    :: rankS
    integer                    :: j

    S(:,:) = 0._fp
    do j = 1,n-rank
       call solve_up_tri( &
            S(1:rank,j), &
            R(1:rank,1:rank), &
            R(1:rank,rank+j), &
            rank, &
            rank )
       S(rank+j,j) = 1._fp
    end do

    do j = 1,p
       call solve_up_tri( &
            t(1:rank,j), &
            R(1:rank,1:rank), &
            c(1:rank,j), &
            rank, &
            rank )
    end do
    t(rank+1:n,1:p) = 0._fp

    call linsolve_qr2( &
         y2, &
         S, &
         t, &
         n, &
         n-rank, &
         p, &
         rankS )

    if ( rankS < n ) STOP 'linsolve_rank_deficient : rank(S) < n'

    do j = 1,p
       call solve_up_tri( &
            y1, &
            R(1:rank,1:rank), &
            c(1:rank,j) - matmul( R(1:rank,rank+1:n), y2(:,j) ), &
            rank, &
            rank )
       y(1:rank,j) = y1
    end do
    y(rank+1:n,1:p) = y2

  end subroutine linsolve_qr_rankdef



end module mod_linalg
