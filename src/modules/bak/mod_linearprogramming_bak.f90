module mod_linearprogramming

  use mod_math

  implicit none

  real(kind=fp), parameter :: EPSlp = epsilon( 1._fp )
  real(kind=fp), parameter :: EPSlpsqr = EPSlp**2
  real(kind=fp), parameter :: BIGlp = real( 1.e6, kind=fp )
  LOGICAL, PARAMETER :: DEBUG = .false.

contains





  recursive subroutine lpp_solve( &
       x, &
       stat, &
       A, &
       c, &
       dim, &
       n )
    ! Solves the d-dimensional linear programming problem
    ! min c.x subject to A(:,1:dim)*x + A(:,dim+1) >= 0
    ! returns 'stat' = 0 if the problem is bounded and feasible
    !         'stat' =-1 if the problem is unbounded
    !         'stat' = 1 if the problem is not feasible
    !         'stat' = 2 if an unexpected error is encountered
    implicit none
    integer,       intent(in)    :: dim, n
    real(kind=fp), intent(in)    :: A(n,dim+1)
    real(kind=fp), intent(in)    :: c(dim)
    real(kind=fp), intent(inout) :: x(dim)
    integer,       intent(out)   :: stat
    logical                      :: mask(dim)
    integer                      :: j(dim-1)
    real(kind=fp)                :: Ai(n,dim), ci(dim-1), xj(dim-1), inv_Ail
    integer                      :: i, k, l
    
    stat = 0

    !PRINT *,'-----------lpp_solve-----------'
    !PRINT *,'X =',X

    if ( dim == 1 ) then
       call lpp_solve_1d( &
            x(1), &
            stat, &
            A, &
            c(1), &
            n )
       !call lpsolve_1d( &
       !     x(1), &
       !     stat, &
       !     A, &
       !     c(1) )
       return
    end if

    do i = dim+1,n
       !IF ( I == 638 .AND. dot_product( A(i,1:dim), x ) + A(i,dim+1) > (dim+1)*EPSlp ) THEN
       !   PRINT *,'*I =',I,', X =',X
       !   PRINT *,'*I =',I,', Ai =',A(i,:)
       !   PRINT *,'*I =',I,', AX + B =',dot_product( A(i,1:dim), x ) + A(i,dim+1)
       !END IF
       !PRINT *,'I =',I,', Ai =',A(i,:)
       !PRINT *,'I =',I,', AX + B =',dot_product( A(i,1:dim), x ) + A(i,dim+1)
       if ( dot_product( A(i,1:dim), x ) + A(i,dim+1) > (dim+1)*EPSlp ) cycle
       ! the provisional optimal point x is on the wrong side of the hyperplane 
       ! associated to the i-th constraint
       !PRINT *,'I =',I,', X =',X
       !PRINT *,'I =',I,', Ai =',A(i,:)
       !PRINT *,'I =',I,', AX + B =',dot_product( A(i,1:dim), x ) + A(i,dim+1)
       ! find the largest coefficient of that hyperplane's equation
       l = maxloc( abs(A(i,1:dim)), DIM=1 )
       !PRINT *,'L =',L

       if ( abs(A(i,l)) < EPSlp ) then
          !PRINT *,'lpp_solve : /!\ AI =',A(I,:)
          cycle
       end if
       
       ! eliminate that variable
       mask(:) = .true.
       mask(l) = .false.
       j = pack( [(k,k=1,dim)], mask )
       !PRINT *,'J =',J

       ! project constraints up to i-1 and the objective function to lower dimension
       inv_Ail = 1._fp / A(i,l)
       do k = 1,i-1
          Ai(k,:) = A(k,[j,dim+1]) - A(k,l) * A(i,[j,dim+1]) * inv_Ail
       end do
       ci = c(j) - c(l) * A(i,j) * inv_Ail

       !PRINT *,'CI =',CI
       !PRINT *,'AI ='
       !CALL PRINT_MAT( AI(1:I-1,:) )

       ! solve lower-dimensional LP subproblem
       xj = x(j)
       !PRINT *,'XJ =',XJ
       
       call lpp_solve( &
            xj, &
            stat, &
            Ai(1:i-1,1:dim), &
            ci, &
            dim-1, &
            i-1 )
       !PRINT *,'XJ =',XJ
       
       if ( stat > 0 ) return ! the problem is unfeasible
       ! back substitution
       x(j) = xj
       x(l) = - ( A(i,dim+1) + dot_product( A(i,j), xj ) ) * inv_Ail
       
    end do


  end subroutine lpp_solve







  
  subroutine lpp_solve_1d( &
       x, &
       stat, &
       A, &
       c, &
       n)
    ! Solves the one-dimensional linear programming problem
    ! min c*x subject to A(:,1)*x + A(:,2) >= 0
    ! returns 'stat' = 0 if the problem is bounded, feasible
    !         'stat' =-1 if the problem is unbounded, feasible
    !         'stat' = 1 if the problem is not feasible
    implicit none
    integer,       intent(in)  :: n
    real(kind=fp), intent(in)  :: A(n,2)
    real(kind=fp), intent(in)  :: c
    real(kind=fp), intent(out) :: x
    integer,       intent(out) :: stat
    real(kind=fp)              :: L, R
    integer                    :: i, np

    !PRINT *,' --------------LPP_SOLVE_1D --------------'
    !PRINT *,'C =',C
    !PRINT *,'A ='
    !CALL PRINT_MAT( A )
    

    stat = 0

    L = -BIGlp
    R =  BIGlp

    np = 0

    do i = 1,n
       if ( A(i,1) < -EPSlp ) then
          np = np + 1
          R = min( R, -A(i,2)/A(i,1) )
       elseif ( A(i,1) > EPSlp ) then
          L = max( L, -A(i,2)/A(i,1) )
       end if
    end do
    
    !PRINT *,'LPP_SOLVE_1D : L=',L,', R=',R
    if ( np == n ) then
       x = L
       if ( c < -EPSlp ) stat = -1
    elseif ( np == 0 ) then
       x = R
       if ( c > EPSlp ) stat = -1
    else
       if ( L > R + EPSlp ) then
          !PRINT *,'LPP_SOLVE_1D : UNFEASBILE'
          stat = 1
       else
          if ( c < -EPSlp ) then
             x = R
          elseif ( c > EPSlp ) then
             x = L
          else
             if ( abs(L) > abs(R) + EPSlp ) then
                x = R
             else
                x = L
             end if
          end if
       end if
    end if
    
  end subroutine lpp_solve_1d





















  recursive subroutine lpsolve( &
       x, &
       stat, &
       A, &
       c )
    ! Solves the d-dimensional linear programming problem
    ! min c.x subject to A(:,1:dim)*x + A(:,dim+1) >= 0
    ! returns 'stat' = 0 if the problem is bounded and feasible
    !         'stat' =-1 if the problem is unbounded
    !         'stat' = 1 if the problem is not feasible
    !         'stat' = 2 if an unexpected error is encountered
    implicit none
    real(kind=fp), intent(in)  :: A(:,:)
    real(kind=fp), intent(in)  :: c(:)
    real(kind=fp), intent(out) :: x(size(c))
    integer,       intent(out) :: stat
    integer                    :: dim, i, j(size(c)-1), k, l
    logical                    :: singular, mask(size(c))
    real(kind=fp)              :: normxsqr, Ai(size(A,1),size(c)), ci(size(c)-1), xj(size(c)-1)
    real(kind=fp)              :: invA_il

    dim = size( c )

    if ( dim == 1 ) then
       call lpsolve_1d( x(1), stat, A, c(1) )
       return
    end if

    stat = 0

    ! deal with secondary constraints
    call solve_NxN( &
         x, &
         A(1:dim,1:dim), &
         -A(1:dim,dim+1), &
         singular )
    if ( singular ) then
       !STOP 'lpsolve : singular secondary constraints.'
       stat = 2
       return
    end if
    normxsqr = sum( x**2 )
    if ( normxsqr > real(dim,kind=fp)*EPSlpsqr ) then
       stat = 2
       return
    else
       x = x / sqrt( normxsqr )
    end if
    !x(:) = 0._fp
    !do i = 1,dim
    !   if ( abs(A(i,i)) > EPSmath ) x(i) = -A(i,dim+1) / A(i,i)
    !end do



    ! deal with primary constraints
    !PRINT *,'A ='
    !CALL PRINT_MAT( A )
    do i = dim+1,size(A,1)
       IF ( DEBUG ) THEN
          PRINT *,'I= ',I
          PRINT *,'X =',X
          PRINT *,A(I,:)
          PRINT *,'AX =',dot_product( A(i,1:dim), x ) + A(i,dim+1)
       END IF
       if ( dot_product( A(i,1:dim), x ) + A(i,dim+1) >= 0._fp ) cycle
       ! the provisional optimum x does not satisfy the i-th constraint
       l = maxloc( abs(A(i,1:dim)), dim=1 )
       IF ( DEBUG ) PRINT *,'KMAX =',L
       if ( abs(A(i,l)) < EPSmath ) then
          !PRINT *,'lpsolve : A = 0'
          stat = 2
          return
       end if
       mask(:) = .true.
       mask(l) = .false.
       j = pack( [(k,k=1,dim)], mask )
       IF ( DEBUG ) PRINT *,'VARS =',J

       invA_il = 1._fp / A(i,l)
       do k = 1,i-1
          Ai(k,:) = A(k,[j,dim+1]) - A(k,l) * A(i,[j,dim+1]) * invA_il
       end do
       ci = c(j) - c(l) * A(i,j) * invA_il

       ! solve lower-dimensional linear programming subproblem
       call lpsolve( &
            xj, &
            stat, &
            Ai(1:i-1,1:dim), &
            ci )

       if ( stat > 0 ) return ! unfeasible problem

       x(j) = xj
       x(l) = - ( A(i,dim+1) + dot_product( A(i,j), xj ) ) * invA_il

    end do

  end subroutine lpsolve








  subroutine lpsolve_1d( &
       x, &
       stat, &
       A, &
       c )
    ! Solves the one-dimensional linear programming problem
    ! min c*x subject to A(:,1)*x + A(:,2) >= 0
    ! returns 'stat' = 0 if the problem is bounded, feasible
    !         'stat' =-1 if the problem is unbounded, feasible
    !         'stat' = 1 if the problem is not feasible
    implicit none
    real(kind=fp), intent(in)  :: A(:,:)
    real(kind=fp), intent(in)  :: c
    real(kind=fp), intent(out) :: x
    integer,       intent(out) :: stat
    real(kind=fp)              :: L, R
    integer                    :: i, n

    IF ( DEBUG ) THEN
       PRINT *,'LP1D'
       print *,'A ='
       call print_mat( A )
       print *,'c ='
       print *,c
    END IF

    R = BIGlp
    L = -BIGlp

    n = 0
    do i = 1,size(A,1)
       if ( A(i,1) < 0._fp ) then
          n = n + 1
          R = min( R, -A(i,2)/A(i,1) )
       elseif ( A(i,1) > 0._fp ) then
          L = max( L, -A(i,2)/A(i,1) )
       end if
    end do

    IF ( DEBUG ) PRINT *,'L =',L,', R =',R

    stat = 0
    if ( n == size(A,1) ) then
       x = R
       if ( c < 0._fp ) stat = -1 ! unbounded, feasible problem
    elseif ( n == 0 ) then
       x = L
       if ( c > 0._fp ) stat = -1 ! unbounded, feasible problem
    else
       if ( L > R ) then
          PRINT *,'LP1D UNFEASIBLE : L, R =',L,R
          stat = 1 ! unfeasible problem
       else
          ! bounded, feasible problem
          if ( c < 0._fp ) then
             x = R
          else
             x = L
          end if
       end if
    end if

  end subroutine lpsolve_1d



end module mod_linearprogramming
