module mod_linprog

  use mod_math

  implicit none

  real(kind=fp), parameter :: EPSlp = epsilon( 1._fp )
  real(kind=fp), parameter :: EPSlpsqr = EPSlp**2
  real(kind=fp), parameter :: BIGlp = real( 1.e6, kind=fp )
  LOGICAL, PARAMETER :: DEBUG = .false.

contains

  recursive subroutine lp_solve( &
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
       call lp_solve_1d( &
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
       
       call lp_solve( &
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


  end subroutine lp_solve







  
  subroutine lp_solve_1d( &
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
    
  end subroutine lp_solve_1d


end module mod_linprog
