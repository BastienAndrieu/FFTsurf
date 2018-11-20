subroutine newton_curve_surface( &
     curv, &
     surf, &
     lowerb, &
     upperb, &
     stat, &
     tuv, &
     toltuv, &
     xyz )
  use mod_math
  use mod_linalg
  use mod_diffgeom
  use mod_tolerances    
  ! stat = 0 : converged
  !        1 : not converged
  !        2 : degeneracy
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .TRUE. )
  logical,       parameter          :: acceleration = .false.
  real(kind=fp), parameter          :: THRESHOLD = real(1.d-2, kind=fp)**2
  integer,       parameter          :: itmax = 30!2 + ceiling(-log10(EPSuv))
  integer,       parameter          :: itconv = 5
  type(type_curve),   intent(in)    :: curv
  type(type_surface), intent(in)    :: surf
  real(kind=fp),      intent(in)    :: lowerb(3)
  real(kind=fp),      intent(in)    :: upperb(3)
  integer,            intent(out)   :: stat
  real(kind=fp),      intent(inout) :: tuv(3)
  real(kind=fp),      intent(out)   :: toltuv
  real(kind=fp),      intent(out)   :: xyz(3)
  real(kind=fp), dimension(3)       :: xyz_c, xyz_s, r
  real(kind=fp)                     :: resxyz, resconv
  real(kind=fp)                     :: jac(3,3), dtuv(3), cond, errtuv, errtuvprev
  integer                           :: it, rank
  logical                           :: linear_conv
  real(kind=fp)                     :: fac
  integer                           :: stat_refl

  IF ( DEBUG ) THEN
     PRINT *,''; PRINT *,'';
     PRINT *,'NEWTON_CURVE_SURFACE'
     PRINT *,'LOWERB =',LOWERB
     PRINT *,'UPPERB =',UPPERB
  END IF

  stat = 1
  errtuv = 0._fp
  errtuvprev = 0._fp
  cond = 1._fp
  fac = 1._fp
  linear_conv = .false.
  !IF ( DEBUG ) PRINT *,'|XS - XC|, |DTUV|, EPS*COND(J)'
  do it = 1,itmax
     !! compute residual
     call eval(xyz_c, curv, tuv(1)  ) ! curve's position vector
     call eval(xyz_s, surf, tuv(2:3)) ! surface's position vector

     r = xyz_s - xyz_c

     ! check signs of convergence
     resxyz = sum(r**2)
     IF ( DEBUG ) PRINT *,SQRT(RESXYZ), SQRT(ERRTUV), EPSFP*COND
     if ( it == 1 ) resconv = THRESHOLD * resxyz
     if ( it > itconv ) then
        if ( resxyz > resconv ) then
           ! Newton sequence not likely to converge, presumably no solution
           IF ( DEBUG ) PRINT *,'NO SIGN OF CONVERGENCE'
           return
        elseif ( errtuv > THRESHOLD * errtuvprev ) then
           ! linear convergence
           IF ( DEBUG ) PRINT *,'LINEAR CONVERGENCE, MU =',sqrt(errtuv/errtuvprev)
           linear_conv = .true.
           fac = 2._fp
        else
           ! superlinear convergence
           fac = 1._fp
           linear_conv = .false.
        end if
     end if

     !! compute Jacobian matrix
     call evald1(jac(:,1), curv, tuv(1))
     jac(:,1) = -jac(:,1)
     call evald1(jac(:,2), surf, tuv(2:3), 1)
     call evald1(jac(:,3), surf, tuv(2:3), 2)
     IF ( DEBUG ) THEN
        PRINT *,'JAC ='
        CALL PRINT_MAT(JAC)
        PRINT *,'RHS =',-R
     END IF

     !! solve for Newton step
     call linsolve_svd( &
          dtuv, &
          jac, &
          -r, &
          3, &
          3, &
          1, &
          cond, &
          rank, &
          tol=EPSmath )
     !PRINT *,'RANK =',RANK
     if ( acceleration ) then
        !PRINT *,'FAC =',FAC
        if ( mod(it,3) == 0 ) dtuv = fac * dtuv
     end if
     errtuvprev = errtuv
     errtuv = sum(dtuv**2)

     ! correct Newton step to keep the iterate inside feasible region
     IF ( DEBUG ) THEN
        PRINT *,' TUV =',TUV
        PRINT *,'DTUV =',DTUV
     END IF
     if ( .true. ) then
        if ( errtuv > sum((upperb - lowerb)**2) ) then
           dtuv = sqrt(sum((upperb - lowerb)**2) / errtuv) * dtuv
           IF ( DEBUG ) PRINT *,'DTUV*=',DTUV
        end if
     end if
     call nd_box_reflexions( &
          tuv, &
          lowerb, &
          upperb, &
          dtuv, &
          3, &
          stat_refl )
     if ( stat_refl > 0 ) return
     IF ( .false. ) PRINT *,'...OK'

     ! update solution
     tuv = tuv + dtuv

     !! termination criterion
     if ( errtuv < max(EPSuvsqr, EPSfpsqr*cond**2) ) then
        if ( resxyz < EPSxyzsqr ) then
           if ( errtuv > EPSuvsqr ) then
              IF ( DEBUG ) PRINT *,'newton_curve_surface : /!\ toltuv > EPSuv'
              if ( linear_conv ) then
                 IF ( DEBUG ) THEN
                    PRINT *,'--> NEWTON_CS_DYNSHIFT'
                    PRINT *,'TUV =',TUV
                 END IF
                 call newton_curve_surface_dynshift( &
                      curv, &
                      surf, &
                      lowerb, &
                      upperb, &
                      stat, &
                      tuv, &
                      toltuv, &
                      xyz )
                 return
              end if
           end if
           ! converged to a curve-surface intersection point
           stat = 0
           toltuv = max(EPSmath, sqrt(errtuv))
           call eval(xyz_c, curv, tuv(1)  )
           call eval(xyz_s, surf, tuv(2:3))
           xyz = 0.5_fp * (xyz_c + xyz_s)
           IF ( DEBUG ) PRINT *,SQRT(RESXYZ), SQRT(ERRTUV), EPSFP*COND
           IF ( DEBUG ) PRINT *,'CONVERGED, TUV =',TUV,', XYZ =',XYZ
           !PAUSE
        else
           IF ( DEBUG ) PRINT *,'STAGNATION'
        end if
        return
     end if

  end do

end subroutine newton_curve_surface
