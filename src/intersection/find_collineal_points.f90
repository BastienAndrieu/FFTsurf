subroutine find_collineal_points( &
    surf, &
    uv, &
    lowerb, &
    upperb, &
    stat )
  use mod_math
  use mod_diffgeom
  use mod_tolerances
  ! Searches for a pair of collineal points on two rectangular parametric surfaces
  ! using a box-constrained Newton-Raphson algorithm. 
  ! The lower (resp. upper) bounds of the 4-dimensional feasible domain are stored 
  ! in 'lowerb' (resp. 'upperb').
  ! Returns 'stat' > 0 if no such pair has been found
  !                  ('stat' = 2 if the Jacobian matrix becomes (near-)singular);
  !         'stat' =-1 if a tangential contact point has been found;
  !         'stat' = 0 else (pair of non-coincident collineal points);
  ! as well as the uv-coordniates of those potential points.
  implicit none
  logical,                      parameter     :: VERBOSE = .false.
  integer,                      parameter     :: nitmax = 10
  type(ptr_parametric_surface), intent(in)    :: surf(2)
  real(kind=MATHpr),            intent(inout) :: uv(2,2)
  real(kind=MATHpr),            intent(in)    :: lowerb(4)
  real(kind=MATHpr),            intent(in)    :: upperb(4)
  integer,                      intent(out)   :: stat
  real(kind=MATHpr)                           :: s(3,2), s1(3,2,2), s2(3,3,2)
  real(kind=MATHpr)                           :: n(3), r(3)
  real(kind=MATHpr)                           :: f(4), jac(4,4), duv(4)
  logical                                     :: singular
  real(kind=MATHpr)                           :: lambda
  integer                                     :: it, isurf, ivar

  stat = 1


  !PRINT *,'LOWERB =', LOWERB
  !PRINT *,'UPPERB =', UPPERB
  !PRINT *,'   UV0 =', UV

  do it = 1,nitmax
     IF (VERBOSE) THEN
        PRINT *,''
        PRINT '(A4,1x,I0)','IT.#',IT
     END IF
     
     do isurf = 1,2
        ! position vector s
        call eval( s(:,isurf), surf(isurf)%ptr, uv(:,isurf) )

        ! tangent vectors su, sv
        do ivar = 1,2
           call evald( s1(:,ivar,isurf), surf(isurf)%ptr, uv(:,isurf), ivar )
        end do
     end do

     n = cross( s1(:,1,1), s1(:,2,1) ) ! (pseudo-)normal to surface 1
     r = s(:,1) - s(:,2)

     ! right-hand side
     f(1:2) = matmul( transpose(s1(:,:,2)), n )
     f(3:4) = matmul( transpose(s1(:,:,1)), r )
     IF (VERBOSE) PRINT *,'    |F| =',REAL(NORM2(F))

     ! convergence criterion
     if ( sum(f**2) < EPScollinealsqr ) then
        IF (VERBOSE) PRINT *,'    CONVERGE |F|'
        if ( sum(r**2) < EPSxyzsqr ) then
           stat = -1
        else
           stat = 0
        end if
        return
     end if

     ! terminate if Newton step becomes too small
     if ( it > 1 ) then
        if ( sum(duv(1:2)**2) < EPSuvsqr .and. sum(duv(3:4)**2) < EPSuvsqr ) then
           IF (VERBOSE) PRINT *,'    |DUV| << 1'
           return
        end if
     end if

     ! second-order partial derivatives
     do isurf = 1,2
        do ivar = 1,3
           call evald2( s2(:,ivar,isurf), surf(isurf)%ptr, uv(:,isurf), ivar )
        end do
     end do

     ! compute Jacobian matrix
     do ivar = 1,2
        jac(1:2,ivar) = matmul( &
             transpose(s1(:,:,2)), &
             cross( s2(:,ivar,1), s1(:,2,1) ) + cross( s1(:,1,1), s2(:,ivar+1,1) ) &
             )
        jac(3:4,ivar) = matmul( &
             transpose(s2(:,ivar:ivar+1,1)), r ) + &
             matmul( transpose(s1(:,:,1)), s1(:,ivar,1) &
             )
        jac(1:2,2+ivar) = matmul( transpose(s2(:,ivar:ivar+1,2)), n )
        jac(3:4,2+ivar) = -matmul( transpose(s1(:,:,1)), s1(:,ivar,2) )
     end do

     ! solve for Newton step
     call linsolve_QR( &
          duv, &
          jac, &
          -f, &
          4, &
          4, &
          singular )
     !PRINT *,'JAC='
     !CALL PRINT_MAT( REAL(JAC) )
     !PRINT *,'  F=', REAL(F)
     !PRINT *,'DUV=', REAL(DUV)

     if ( singular ) then
        IF (VERBOSE) PRINT *,' find_collineal_points : singular jacobian matrix'
        stat = 2
        return
     end if

     ! scale down Newton step to keep the solution inside feasible region
     call nd_box_constraint( &
          reshape( uv, [4] ), &
          lowerb, &
          upperb, &
          duv, &
          lambda )
     !PRINT *,'LAMBDA=',LAMBDA
     !PRINT *,''

     if ( lambda < -MATHeps ) then
        IF (VERBOSE) PRINT *,' find_collineal_points : negative Newton step rescaling, &
             & lambda =', lambda
        stat = 3
        return
     end if

     ! update solution
     duv = lambda * duv
     uv(:,1) = uv(:,1) + duv(1:2)
     uv(:,2) = uv(:,2) + duv(3:4)
     IF (VERBOSE) PRINT *,'    UV =', REAL(UV)

  end do

end subroutine find_collineal_points
