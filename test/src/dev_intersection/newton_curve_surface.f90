subroutine newton_curve_surface( &
     curv, &
     surf, &
     lowerb, &
     upperb, &
     stat, &
     tuv, &
     xyz )
  use mod_math
  use mod_linalg
  use mod_diffgeom2
  use mod_tolerances    
  ! stat = 0 : converged
  !        1 : not converged
  !        2 : degeneracy
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  real(kind=fp), parameter          :: THRESHOLD = real(1.d-2, kind=fp)
  integer,       parameter          :: itmax = 20!2 + ceiling(-log10(EPSuv))
  integer,       parameter          :: itconv = 5

  type(type_curve),   intent(in)    :: curv
  type(type_surface), intent(in)    :: surf
  real(kind=fp),      intent(in)    :: lowerb(3)
  real(kind=fp),      intent(in)    :: upperb(3)
  integer,            intent(out)   :: stat
  real(kind=fp),      intent(inout) :: tuv(3)
  real(kind=fp),      intent(out)   :: xyz(3)
  real(kind=fp), dimension(3)       :: xyz_c, xyz_s, r
  real(kind=fp)                     :: resxyz, resconv
  real(kind=fp)                     :: jac(3,3), dtuv(3), cond, errtuv!, lambda
  integer                           :: rank
  integer                           :: it

  IF ( DEBUG ) THEN
     PRINT *,''; PRINT *,'';
     PRINT *,'NEWTON_CURVE_SURFACE'
     PRINT *,'LOWERB =',LOWERB
     PRINT *,'UPPERB =',UPPERB
  END IF

  stat = 1
  errtuv = 0._fp
  cond = 1._fp

  do it = 1,itmax
     !IF ( DEBUG ) PRINT *,'IT #',IT
     !! compute residual
     call eval(xyz_c, curv, tuv(1)  ) ! curve's position vector
     call eval(xyz_s, surf, tuv(2:3)) ! surface's position vector

     r = xyz_s - xyz_c

     !IF ( DEBUG ) PRINT *,'R =',R

     ! check signs of convergence
     resxyz = sum(r**2)
     IF ( DEBUG ) PRINT *,SQRT(RESXYZ), SQRT(ERRTUV), EPSUV*COND
     if ( it == 1 ) resconv = THRESHOLD * resxyz
     if ( it > itconv .and. resxyz > resconv ) then
        ! Newton sequence not likely to converge, presumably no solution
        IF ( DEBUG ) PRINT *,'NO SIGN OF CONVERGENCE'
        return
     end if
     
     !! termination criteria
     if ( errtuv < EPSuvsqr*cond**2 .and. it > 1 ) then
        if ( resxyz < EPSxyzsqr ) then
           ! converged to a curve-surface intersection point
           stat = 0
           xyz = 0.5_fp * (xyz_c + xyz_s)
           IF ( DEBUG ) PRINT *,'CONVERGED, TUV =',TUV,', XYZ =',XYZ
        else
           IF ( DEBUG ) PRINT *,'STAGNATION'
        end if
        return
     end if

     !if ( resxyz < EPSxyzsqr .and. errtuv < EPSuvsqr*cond**2 ) then ! <---+
     !   ! converged to a curve-surface intersection point                 !
     !   stat = 0                                                          !
     !   xyz = 0.5_fp * (xyz_c + xyz_s)                                    !
     !   IF ( DEBUG ) PRINT *,'CONVERGED, TUV =',TUV,', XYZ =',XYZ
     !   return                                                            !
     !end if ! <-----------------------------------------------------------+

     !! compute Jacobian matrix
     call evald1(jac(:,1), curv, tuv(1))
     jac(:,1) = -jac(:,1)
     call evald1(jac(:,2), surf, tuv(2:3), 1)
     call evald1(jac(:,3), surf, tuv(2:3), 2)

     !! solve for Newton step
     call linsolve_svd( &
          dtuv, &
          jac, &
          -r, &
          3, &
          3, &
          1, &
          cond, &
          rank )

     errtuv = sum( dtuv**2 )

     if ( rank < 3 ) then ! <------------------+
        ! singular Jacobian matrix             !
        if ( stat == 0 ) then ! <----+         !
           stat = 2                  !         !
           IF ( DEBUG ) PRINT *,'SINGULAR JACOBIAN AT SOLUTION'
           return                    !         !
        end if ! <-------------------+         !
     end if ! <--------------------------------+

     ! damp Newton step to keep the solution inside feasible region
     !call nd_box_constraint( &
     !     tuv, &
     !     lowerb, &
     !     upperb, &
     !     dtuv, &
     !     lambda )
     call nd_box_reflexions( &
          tuv, &
          lowerb, &
          upperb, &
          dtuv, &
          3 )

     !IF ( DEBUG ) PRINT *,'DAMPING =',lambda

     !if ( lambda < epsilon(1._fp) ) then ! <---+
     !   ! non-positive scaling factor          !
     !   IF ( DEBUG ) PRINT *,'NONPOSITIVE DAMPING FACTOR'
     !   return                                 !
     !end if ! <--------------------------------+

     ! update solution
     !dtuv = lambda * dtuv
     tuv = tuv + dtuv

  end do

end subroutine newton_curve_surface
