module mod_linearprogramming

  use mod_math

  implicit none

  real(kind=MATHpr), parameter :: OMGlp = real( 1.e6, kind=MATHpr )

contains


  recursive subroutine lpsolve( &
       x, &
       stat, &
       A, &
       c )
    ! Solves the d-dimensional linear programming problem
    ! min c.x subject to A(:,1)*x + A(:,2) > 0
    ! returns 'stat' = 0 if the problem is bounded and feasible
    !         'stat' =-1 if the problem is unbounded
    !         'stat' = 1 if the problem is not feasible
    !         'stat' = 2 if an unexpected error is encountered
    implicit none
    real(kind=MATHpr), intent(in)  :: A(:,:)
    real(kind=MATHpr), intent(in)  :: c(:)
    real(kind=MATHpr), intent(out) :: x(size(c))
    integer,           intent(out) :: stat
    integer                        :: dim, i, j(size(c)-1), k, l
    logical                        :: singular, mask(size(c))
    real(kind=MATHpr)              :: Ai(size(A,1),size(c)), ci(size(c)-1), xj(size(c)-1)
    real(kind=MATHpr)              :: invA_il

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

    ! deal with primary constraints
    do i = dim+1,size(A,1)
       if ( dot_product( A(i,1:dim), x ) + A(i,dim+1) <= 0._MATHpr ) cycle
       
       ! the provisional optimum x does not satisfy the i-th constraint
       l = maxloc( abs(A(i,1:dim)), 1 )
       if ( abs(A(i,l)) < MATHeps ) then
          !PRINT *,'lpsolve : A = 0'
          stat = 2
          return
       end if
       mask(:) = .true.
       mask(l) = .false.
       j = pack( [(k,k=1,dim)], mask )

       invA_il = 1._MATHpr / A(i,l)
       do k = 1,i-1
          Ai(k,:) = A(k,[j,dim+1]) - A(k,l) * A(i,[j,dim+1]) * invA_il
       end do
       ci = c(j) - c(l) * A(i,j) * invA_il

       ! solve lower-dimensional linear programming subproblem
       call lpsolve( &
            xj, &
            stat, &
            Ai(1:i-1,:), &
            ci )

       if ( stat > 0 ) return ! unfeasible problem

       x(j) = xj
       x(l) = - ( A(i,dim+1) + dot_product( A(i,j), xj ) )  * invA_il

    end do

  end subroutine lpsolve








  subroutine lpsolve_1d( &
       x, &
       stat, &
       A, &
       c )
    ! Solves the one-dimensional linear programming problem
    ! min c*x subject to A(:,1)*x + A(:,2) > 0
    ! returns 'stat' = 0 if the problem is bounded, feasible
    !         'stat' =-1 if the problem is unbounded, feasible
    !         'stat' = 1 if the problem is not feasible
    implicit none
    real(kind=MATHpr), intent(in)  :: A(:,:)
    real(kind=MATHpr), intent(in)  :: c
    real(kind=MATHpr), intent(out) :: x
    integer,           intent(out) :: stat
    real(kind=MATHpr)              :: L, R
    integer                        :: i, n

    R = OMGlp
    L = -OMGlp

    n = 0
    do i = 1,size(A,1)
       if ( A(i,1) > 0._MATHpr ) then
          n = n + 1
          R = min( R, -A(i,2)/A(i,1) )
       elseif ( A(i,1) < 0._MATHpr ) then
          L = max( L, -A(i,2)/A(i,1) )
       end if
    end do

    stat = 0
    if ( n == size(A,1) ) then
       x = R
       if ( c < 0._MATHpr ) stat = -1 ! unbounded, feasible problem
    elseif ( n == 0 ) then
       x = L
       if ( c > 0._MATHpr ) stat = -1 ! unbounded, feasible problem
    else
       if ( L > R ) then
          !PRINT *,'LP1D UNFEASIBLE : L, R =',L,R
          stat = 1 ! unfeasible problem
       else
          ! bounded, feasible problem
          if ( c < 0._MATHpr ) then
             x = R
          else
             x = L
          end if
       end if
    end if

  end subroutine lpsolve_1d

end module mod_linearprogramming
