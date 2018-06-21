program newton_curve_surface_singular
  
  use mod_math
  use mod_polynomial
  use mod_diffgeom2

  implicit none
  
  type(type_curve)   :: curv
  type(type_surface) :: surf
  real(kind=fp)      :: tuv(3), xyz(3), dxyz_dt(3), dxyz_duv(3,2), n(3)
  integer            :: stat, ivar


  
  call read_polynomial( &
       curv%x, &
       'newton_curve_surface_singular/curv.cheb', &
       1, &
       1 )

  !CALL PRINT_MAT( CURV%X%COEF(:,:,1) )
  !STOP




  call read_polynomial( &
       surf%x, &
       'newton_curve_surface_singular/surf.cheb', &
       2, &
       1 )

  call compute_deriv1( curv )
  call compute_deriv2( curv )

  call compute_deriv1( surf )
  call compute_deriv2( surf )

  tuv(:) = 0._fp
  call newton_curve_surface( &
       curv, &
       surf, &
       [-1._fp, 1._fp], &
       [-1._fp, 1._fp, -1._fp, 1._fp], &
       tuv, &
       stat, &
       xyz )

  if ( stat == 0 ) then
     PRINT *,'TUV =',tuv

     call evald1( &
          dxyz_dt, &
          curv, &
          tuv(1) )

     do ivar = 1,2
        call evald1( &
             dxyz_duv(:,ivar), &
             surf, &
             tuv(2:3), &
             ivar )
     end do
     n = cross( dxyz_duv(:,1), dxyz_duv(:,2) )

     PRINT *,'   XYZ =',XYZ
     PRINT *,'|GT.N| =',abs( dot_product( dxyz_dt, n ) )

  end if
  

contains

  subroutine newton_curve_surface( &
       curv, &
       surf, &
       tbox, &
       uvbox, &
       tuv, &
       stat, &
       xyz )
    use mod_math
    use mod_diffgeom2
    use mod_tolerances    
    ! stat = 0 : converged
    !        1 : not converged
    !        2 : degeneracy
    implicit none
    real(kind=fp), parameter          :: THRESHOLD = real(1.e-2, kind=fp)
    real(kind=fp), parameter          :: EPS = EPSmath !real(1.e-12, kind=fp)
    real(kind=fp), parameter          :: TOL = 1.d-15!EPSxyz
    real(kind=fp), parameter          :: EPSsqr = EPS**2
    real(kind=fp), parameter          :: TOLsqr = TOL**2
    integer, parameter                :: nitmax = 30!ceiling(-log10(EPS))
    integer, parameter                :: nitcheck = 5

    type(type_curve),   intent(in)    :: curv
    type(type_surface), intent(in)    :: surf
    real(kind=fp),      intent(in)    :: tbox(2)
    real(kind=fp),      intent(in)    :: uvbox(4)
    real(kind=fp),      intent(inout) :: tuv(3)
    integer,            intent(out)   :: stat
    real(kind=fp),      intent(out)   :: xyz(3)
    real(kind=fp), dimension(3)       :: lowerb, upperb, rng
    real(kind=fp), dimension(3)       :: xyz_c, xyz_s, r
    real(kind=fp), dimension(3)       :: tuvtmp, dtuv
    real(kind=fp)                     :: rescheck, res, restmp
    real(kind=fp)                     :: jac(3,3), lambda
    logical                           :: singular
    integer                           :: rank
    integer                           :: it

    ! Feasible t,u,v-domain
    lowerb = [ tbox(1), uvbox([1,3]) ]
    upperb = [ tbox(2), uvbox([2,4]) ]
    rng = upperb - lowerb
    lowerb = lowerb - EPsuv*rng
    upperb = upperb + EPsuv*rng

    stat = 1
    restmp = huge(1._fp)
    newton_iteration : do it = 1,nitmax

       ! position vector
       call eval( xyz_c, curv, tuv(1) )   ! curve
       call eval( xyz_s, surf, tuv(2:3) ) ! surface

       r = xyz_s - xyz_c ! residual vector
       res = sum( r**2 ) ! squared norm of residual vector
       !PRINT *,'NEWTON, IT#',IT,', RES =',REAL(NORM2(R))
       PRINT *,NORM2(R)

       ! check signs of convergence
       if ( it == 1 ) rescheck = THRESHOLD * res
       if ( it > nitcheck .and. res > rescheck ) then
          !   PRINT *,'NO SIGN OF CONVERGENCE, STOP NEWTON ITERATION'
          return
       end if

       ! convergence criterion
       if ( res < TOLsqr ) then
          stat = 0
          if ( res < restmp ) then
             restmp = res
             tuvtmp = tuv
             xyz = 0.5_fp * ( xyz_s + xyz_c )
          end if
          if ( restmp < EPSsqr ) then
             tuv = tuvtmp
             return
          end if
       end if

       ! Jacobian matrix
       call evald1( jac(:,1), curv, tuv(1) )
       jac(:,1) = -jac(:,1)
       call evald1( jac(:,2), surf, tuv(2:3), 1 )
       call evald1( jac(:,3), surf, tuv(2:3), 2 )

       ! solve for Newton step
       call linsolve_QR( &
            dtuv, &
            jac, &
            -r, &
            3, &
            3, &
            singular, &
            rank )

       PRINT *,'DET(J) =',triple_product( jac(:,1),jac(:,2),jac(:,3) )

       if ( rank < 3 ) then ! degeneracy
          PRINT *,''
          PRINT *,'IT #',IT
          PRINT *,'TUV =',TUV
          PRINT *,'RES =',SQRT(RES)
          PRINT *,'JACOBIAN ='
          CALL PRINT_MAT( JAC )
          PRINT *,'RANK =',RANK
          PRINT *,'DTUV=',DTUV
          if ( stat == 0 ) then
             stat = 2
             return
          end if
       end if

       ! scale down Newton step to keep the solution inside feasible region
       call nd_box_constraint( &
            tuv, &
            lowerb, &
            upperb, &
            dtuv, &
            lambda )

       if ( lambda < -EPSmath ) return ! negative damped factor

       dtuv = lambda * dtuv
       if ( abs(dtuv(1)) < EPSuv .or. sum(dtuv(2:3)**2) < EPSuvsqr ) then
          ! damped Newton step is too small
          !PRINT *,'|DT| =',ABS(DTUV(1)),' |DUV| =',NORM2( DTUV(2:3) )
          return
       end if

       ! update solution
       tuv = tuv + lambda * dtuv

    end do newton_iteration

  end subroutine newton_curve_surface

end program
