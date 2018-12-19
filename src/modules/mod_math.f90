module mod_math

  implicit none

  
  integer,       parameter :: fp = SELECTED_REAL_KIND(15,307)!kind(1.d0)

  real(kind=fp), parameter :: EPSfp = epsilon( 1._fp )
  real(kind=fp), parameter :: EPSfpsqr = EPSfp**2
  real(kind=fp), parameter :: EPSmath = 10._fp * EPSILON( 1._fp )
  real(kind=fp), parameter :: CSTpi = 4._fp * atan(1._fp)

  ! ------- à supprimer -------
  integer,           parameter :: MATHpr = SELECTED_REAL_KIND(15,307)!kind(1.d0) 
  real(kind=MATHpr), parameter :: MATHeps = 10._MATHpr * epsilon(1._MATHpr)!2.d-15
  real(kind=MATHpr), parameter :: MATHpi = 3.14159265358979323846d0
! ---------------------------
  
  ! --------------------------------------------------------------------------------------------------------------------------- !
  ! Workaround for making arrays of matrices (rank-2 arrays) with different sizes (useful for sequences of matrices)            !
  type type_matrix                                                                                                              !
     real(kind=fp), allocatable :: mat(:,:)                                                                                     !
  end type type_matrix                                                                                                          !
  ! --------------------------------------------------------------------------------------------------------------------------- !




  interface cross
     module procedure cross_1!, cross_n, cross_mn
  end interface

  !interface norm2
  !   module procedure norm2_1, norm2_n!, norm2_mn
  !end interface

  interface print_mat
     module procedure print_mat_I, print_mat_R, print_mat_D
  end interface print_mat

  
contains

  function cross_1( &
       u, &
       v) result(w)

    implicit none

    real(kind=fp) :: u(3)
    real(kind=fp) :: v(3)
    real(kind=fp) :: w(3)
    integer             :: i, j, k

    do i = 1,3
       j = 1 + mod(i,3)
       k = 1 + mod(j,3)
       w(i) = u(j)*v(k) - u(k)*v(j)
    end do

  end function cross_1




  function norm2_1( &
       u) result(m)

    implicit none

    real(kind=fp) :: u(:)
    real(kind=fp) :: m

    m = sqrt(sum(u**2))

  end function norm2_1

  

  function norm2_n( &
       u) result(m)

    implicit none

    real(kind=fp) :: u(:,:)
    real(kind=fp) :: m(size(u,2))

    m = sqrt(sum(u**2,1))

  end function norm2_n



  function triple_product( a, b, c ) result(d)
    implicit none
    real(kind=fp), dimension(3) :: a, b, c
    real(kind=fp)               :: d

    d = dot_product( a, cross( b, c ) )

  end function triple_product
  

  



  ! =======================================
  !           Algèbre linéaire
  ! =======================================


  ! -----------------------------------------
  ! solve_up_tri
  ! Résolution de système triangulaire supérieur
  subroutine solve_up_tri(x, U, b, M, N)
    implicit none

    integer, intent(in) :: M, N
    real(kind=fp), dimension(M,N), intent(in) :: U
    real(kind=fp), dimension(M), intent(in) :: b
    real(kind=fp), dimension(N), intent(out) :: x

    integer :: i

    x(N) = b(N)/U(N,N)

    do i = N-1,1,-1
       x(i) = ( b(i) - sum( x(i+1:N)*U(i,i+1:N) ) )/U(i,i)
    end do


  end subroutine solve_up_tri
  ! -----------------------------------------


  ! -----------------------------------------
  ! solve_NxN
  ! Résolution de système linéaire avec pivot de Gauss
  subroutine solve_NxN(x,A,b,singular)

    implicit none

    real(kind=fp), dimension(:,:), intent(in) :: A
    real(kind=fp), dimension(:), intent(in) :: b
    real(kind=fp), dimension(size(A,2)), intent(out) :: x
    logical, intent(out), optional :: singular

    integer :: M, N
    real(kind=fp), dimension(size(A,1), size(A,2)+1) :: C
    REAL(kind=fp) :: DET, res
    INTEGER :: I


    M = size(A,1)
    N = size(A,2)

    C(:,1:N) = A
    C(:,N+1) = b

    call gauss_elim(C, singular)

    if (present(singular)) then
       DET = 1.0D0
       DO I = 1,MIN(M,N)
          DET = DET * C(I,I)
       END DO
       IF (ABS(DET) < MATHeps) THEN
          SINGULAR = .TRUE.
       END IF

       !if (singular) print *,'solve_NxN : matrice singulière, det =', DET
    end if

    call solve_up_tri(x, C(:,1:N), C(:,N+1), M, N)


    ! Vérif
    res = maxval(abs(matmul(A,x)-b))
    !if (res > 1.0d-4) PRINT *,'/!\ solve_NxN : res =',res

  end subroutine solve_NxN
  ! -----------------------------------------


  ! -----------------------------------------
  ! gauss_elim
  ! Pivot de Gauss
  subroutine gauss_elim(A, singular)

    implicit none

    real(kind=fp), dimension(:,:), intent(inout) :: A
    logical, intent(out), optional :: singular

    integer :: m, n, i, k, i_max

    m = size(A,1) ; n = size(A,2)

    if (present(singular)) singular = .false.

    do k = 1,min(m,n)
       i_max = k-1 + maxloc(abs(A(k:m,k)),1) ! pivot

       if (present(singular)) then ! vérification singularité
          if (abs(A(i_max,k)) < 1.0d-15) singular = .true.
       end if

       call swap_rows(A,i_max,k)

       do i = k+1,m
          A(i,k+1:) = A(i,k+1:) - A(k,k+1:)*A(i,k)/A(k,k)

          A(i,k) = 0.0d0
       end do
    end do

  end subroutine gauss_elim
  ! -----------------------------------------

  ! -----------------------------------------
  ! swap_rows
  ! Permutation de lignes
  subroutine swap_rows(A, i, j)

    implicit none

    real(kind=fp), intent(inout) :: A(:,:)
    integer,       intent(in)    :: i, j
    real(kind=fp)                :: r(size(A,2))

    if (i == j) return

    r = A(i,:)
    A(i,:) = A(j,:)
    A(j,:) = r

  end subroutine swap_rows
  ! -----------------------------------------


  ! -----------------------------------------
  ! solve_2x2
  ! Résolution d'un système linéaire 2x2 
  ! par la méthode de Cramer
  subroutine solve_2x2( &
       x, &
       a11, &
       a12, &
       a21, &
       a22, &
       b, &
       singular)
    implicit none
    real(kind=fp), intent(in)  :: a11, a12, a21, a22, b(2)
    real(kind=fp), intent(out) :: x(2)
    logical,       intent(out) :: singular
    real(kind=fp)              :: det

    det = a11*a22 - a12*a21

    singular = ( abs(det) < MATHeps )
    
    x(1) = (b(1)*a22 - a12*b(2)) / det
    x(2) = (a11*b(2) - b(1)*a21) / det

  end subroutine solve_2x2
  ! -----------------------------------------




  
  ! -----------------------------------------
  ! outer_product
  ! Produit tensoriel de deux vecteurs colonnes
  function outer_product(u,v)

    implicit none

    real(kind=fp), dimension(:) :: u, v
    real(kind=fp)               :: outer_product(size(u),size(v))

    outer_product = spread(u, dim=2, ncopies=size(v)) * spread(v, dim=1, ncopies=size(u))

  end function outer_product
  ! -----------------------------------------







  


  
  ! -----------------------------------------
  subroutine nd_box_constraint( &
       x, &
       lowerb, &
       upperb, &
       dx, &
       lambda )
    ! All possible cases :
    !-------+-------+-----------+
    !   x   ! x+dx  !  lambda   !
    !-------+-------+-----------!
    !  in   !  in   !     1     !
    !  in   !  out  ! 0 < * < 1 !
    !  out  !  in   !     1     !
    !  out  !  out  !    ???    !
    !-------+-------+-----------+
    implicit none
    real(kind=fp), intent(in)  :: x(:)
    real(kind=fp), intent(in)  :: lowerb(:)
    real(kind=fp), intent(in)  :: upperb(:)
    real(kind=fp), intent(in)  :: dx(:)
    real(kind=fp), intent(out) :: lambda
    integer                    :: idim
  
    lambda = 1._fp
    do idim = 1,size(x)
       if ( abs(dx(idim)) < EPSmath ) cycle
       if ( dx(idim) < 0._fp ) then
          lambda = min( lambda, (lowerb(idim) - x(idim)) / dx(idim) )
       else
          lambda = min( lambda, (upperb(idim) - x(idim)) / dx(idim) )
       end if
    end do
    

  end subroutine nd_box_constraint
  ! -----------------------------------------

  ! -----------------------------------------
  subroutine nd_box_reflexions( &
       x, &
       lowerb, &
       upperb, &
       dx, &
       dim, &
       stat )
    implicit none
    LOGICAL, PARAMETER :: DEBUG = .FALSE.
    INTEGER, PARAMETER :: ITMAX = 10
    integer,       intent(in)    :: dim
    real(kind=fp), intent(in)    :: x(dim)
    real(kind=fp), intent(in)    :: lowerb(dim)
    real(kind=fp), intent(in)    :: upperb(dim)
    real(kind=fp), intent(inout) :: dx(dim)
    integer,       intent(out)   :: stat
    real(kind=fp)                :: xtmp(dim), dxtmp(dim), xtmp1(dim)
    real(kind=fp)                :: lambda, lambda_i
    integer                      :: ib, idim, it

    IF ( DEBUG ) THEN
       PRINT *,'IN'
       PRINT *,'LOWERB =',LOWERB
       PRINT *,'UPPERB =',UPPERB
       PRINT *,'     X =',X
       PRINT *,'    DX =',DX
    END IF
    xtmp = x
    dxtmp = dx

    ! check validity of bounds
    if ( any(lowerb - upperb > EPSmath) ) then
       PRINT *,'nd_box_reflexions : invalid bounds'
       stat = 3
       return
    end if

    ! check if x is inside feasible region...
    do idim = 1,dim
       if ( .not.is_in_closed_interval( &
            x(idim), &
            lowerb(idim), &
            upperb(idim), &
            tolerance=EPSmath) ) then
          PRINT *, x(idim), ' not in interval', lowerb(idim), upperb(idim)
          PRINT *,'nd_box_reflexions : starting point out of bounds'
          stat = 2
          return
       end if
    end do

    stat = 1
    it = 0
    do
       it = it + 1
       if ( it > itmax ) then
          xtmp = xtmp1
          dx = xtmp - x
          !STOP 'nd_box_reflexions : it > itmax'
       end if
       lambda = 1._fp
       ib = 0
       do idim = 1,dim
          if ( abs(dxtmp(idim)) < EPSmath ) cycle
          if ( dxtmp(idim) < 0._fp ) then
             lambda_i = min( lambda, (lowerb(idim) - xtmp(idim)) / dxtmp(idim) )
          else
             lambda_i = min( lambda, (upperb(idim) - xtmp(idim)) / dxtmp(idim) )
          end if
          if ( lambda > lambda_i ) then
             ib = idim
             lambda = lambda_i
          end if
       end do

       if ( ib == 0 ) then
          xtmp = xtmp + dxtmp
          dx = xtmp - x
          stat = 0
          return
       else
          xtmp = xtmp + lambda*dxtmp
          if ( it == 1 ) xtmp1 = xtmp
          dxtmp = (1._fp - lambda) * dxtmp
          dxtmp(ib) = -dxtmp(ib)
       end if
       
    end do

  end subroutine nd_box_reflexions
  ! -----------------------------------------










  ! -----------------------------------------
  ! print_mat
  ! print une matrice rectangulaire
  subroutine print_mat_D(A)
    implicit none
    real*8, dimension(:,:), intent(in) :: A
    integer                            :: i

    do i = 1,size(A,1)
       print *,A(i,:)
    end do
    
  end subroutine print_mat_D

  subroutine print_mat_R(A)
    implicit none
    real, dimension(:,:), intent(in) :: A
    integer                          :: i

    do i = 1,size(A,1)
       print *,A(i,:)
    end do

  end subroutine print_mat_R

  subroutine print_mat_I(A)
    implicit none
    integer, dimension(:,:), intent(in) :: A
    integer                             :: i

    do i = 1,size(A,1)
       print *,A(i,:)
    end do

  end subroutine print_mat_I
  ! -----------------------------------------




  
  function Dnchoosek( &
       n, &
       k )
    ! Numerical Recipes in Fortran 77
    implicit none
    integer :: n, k
    real(kind=fp)  :: Dnchoosek

    Dnchoosek = exp( factln(n) - factln(k) - factln(n-k) )

  end function Dnchoosek


  
  function factln( &
       n )

    ! Numerical Recipes in Fortran 77
    implicit none
    
    integer :: n
    real(kind=fp)  :: factln
    real(kind=fp)  :: a(100)

    save a
    data a /100*-1.d0/

    if (n < 0) STOP 'factln : n < 0'
    
    if (n < size(a)) then
       if (a(n+1) < 0) a(n+1) = gammaln( dble(n+1) )
       factln = a(n+1)
    else
       factln = gammaln( dble(n+1) )
    end if
    
  end function factln



  function gammaln( &
       xx )
    
    ! Numerical Recipes in Fortran 77
    implicit none
    
    real(kind=fp)  :: xx
    real(kind=fp)  :: gammaln
    real(kind=fp)  :: ser, stp, tmp, x, y, cof(6)
    integer :: j

    save cof, stp
    data cof, stp /76.18009172947146d0,-86.50532032941677d0, &
         24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
         -.5395239384953d-5,2.5066282746310005d0/

    x = xx
    y = x

    tmp = x + 5.5d0
    tmp = (x + 0.5d0)*log(tmp) - tmp
    ser = 1.000000000190015d0
    do j = 1,6
       y = y + 1.d0
       ser = ser + cof(j)/y
    end do
    gammaln = tmp + log(stp*ser/x)

  end function gammaln


  recursive subroutine myfact( n, f )
    
    implicit none
    
    integer,   intent(in)  :: n
    integer*8, intent(out) :: f
    integer*8              :: fnm1

    if (n < 0) STOP 'myfact : n < 0'

    if (n == 0) then
       f = 1
       return 
    else
       call myfact( n-1, fnm1 )
       f = int8(n) * fnm1
    end if

  end subroutine myfact


  


  function myfactorial( n ) result( f )
    
    implicit none

    integer   :: n
    integer*8 :: f

    if (n > 170) PRINT *,'myfactorial : N > 170 --> OVERFLOW !!!'
    if (n < 21) then
       call myfact( n, f )
    else
       f = int8( exp( factln(n) ) )
    end if

  end function myfactorial




  
  function mean_angle( ang ) 
    implicit none
    real(kind=fp), intent(in) :: ang(:)
    real(kind=fp)             :: mean_angle
    complex(kind=fp)          :: s

    s = sum( dcmplx( cos(ang), sin(ang) ) ) / dble( size(ang) )

    mean_angle = atan2( imagpart(s), realpart(s) )
    
  end function mean_angle


  function diff_angle( a1, a2 )
    implicit none
    real(kind=fp), intent(in) :: a1, a2
    real(kind=fp)             :: diff_angle
    real(kind=fp)             :: c1, s1, c2, s2

    c1 = cos(a1)
    s1 = sin(a1)
    c2 = cos(a2)
    s2 = sin(a2)
    diff_angle = atan2( s1*c2 - c1*s2, c1*c2 + s1*s2 )

  end function diff_angle




  function ab2n1p1( x, a, b )
    ! Applies to a scalar x the linear change of variable that 
    ! maps the interval [a,b] to [-1,1]
    implicit none
    real(kind=fp), intent(in) :: x, a, b
    real(kind=fp)             :: ab2n1p1
    ab2n1p1 = -1._fp + 2._fp * (x - a) / (b - a)
  end function ab2n1p1





  function n1p12ab( x, a, b )
    ! Applies to a scalar x the linear change of variable that 
    ! maps the interval [-1,1] to [a,b]
    implicit none
    real(kind=fp), intent(in) :: x, a, b
    real(kind=fp)             :: n1p12ab
    n1p12ab = a + 0.5_fp * (b-a) * (x + 1._fp)
  end function n1p12ab




  
  function linspace( a, b, n )
    ! Returns a vector of n values linearly ranging from a to b
    implicit none
    real(kind=fp), intent(in) :: a, b
    integer,       intent(in) :: n
    real(kind=fp)             :: linspace(n)
    integer                   :: i = 0

    linspace = a + (b-a) * &
         [( real(i-1, kind=fp) / real(n-1, kind=fp), i=0,n-1 )]
    
  end function linspace





  function is_in_closed_interval( x, a, b, tolerance )
    ! Returns .true. if a <= x <= b, .false. if x < a or x > b
    implicit none
    real(kind=fp),           intent(in) :: x, a, b
    real(kind=fp), optional, intent(in) :: tolerance
    real(kind=fp)                       :: tol
    logical                             :: is_in_closed_interval
    
    if ( present(tolerance) ) then
       tol = tolerance
    else
       tol = EPSfp
    end if

    is_in_closed_interval = ( x >= min(a,b) - tol .and. x <= max(a,b) + tol )

  end function is_in_closed_interval




  function is_in_open_interval( x, a, b, tolerance )
    ! Returns .true. if a < x < b, .false. if x <= a or x >= b
    implicit none
    real(kind=fp), intent(in)           :: x, a, b
    real(kind=fp), optional, intent(in) :: tolerance
    real(kind=fp)                       :: tol
    logical                             :: is_in_open_interval

    if ( present(tolerance) ) then
       tol = tolerance
    else
       tol = EPSfp
    end if

    is_in_open_interval = ( x > min(a,b) + tol .and. x < max(a,b) - tol )

  end function is_in_open_interval

  


  function identity_matrix( n )
    implicit none
    integer, intent(in) :: n
    real(kind=fp)       :: identity_matrix(n,n)
    integer             :: i
    identity_matrix(:,:) = 0._fp
    do i = 1,n
       identity_matrix(i,i) = 1._fp
    end do
  end function identity_matrix



  function diag( v )
    implicit none
    real(kind=fp), intent(in) :: v(:)
    real(kind=fp)             :: diag(size(v), size(v))
    integer                   :: i
    diag(:,:) = 0._fp
    do i = 1,size(v)
       diag(i,i) = v(i)
    end do
  end function diag



  
  elemental function is_zero( x )
    implicit none
    real(kind=fp), intent(in) :: x
    logical                   :: is_zero

    is_zero = (abs(x) < EPSfp)
  end function is_zero




  function overlap_intervals( &
       ab, &
       cd, &
       tolerance )
    implicit none
    real(kind=fp), dimension(2), intent(in) :: ab, cd
    real(kind=fp),     optional, intent(in) :: tolerance
    logical                                 :: overlap_intervals
    real(kind=fp)                           :: tol

    if ( present(tolerance) ) then
       tol = tolerance
    else
       tol = EPSfp
    end if

    overlap_intervals = ( min(maxval(ab), maxval(cd)) >= max(minval(ab), minval(cd)) - tolerance )
    
  end function overlap_intervals


  
end module mod_math
