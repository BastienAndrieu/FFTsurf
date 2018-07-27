subroutine newton_curve_surface_dynshift( &
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
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  integer,       parameter          :: itmax = 20!2 + ceiling(-log10(EPSuv))
  real(kind=fp), parameter          :: shiftd = real(1.d-1, kind=fp)
  type(type_curve),   intent(in)    :: curv
  type(type_surface), intent(in)    :: surf
  real(kind=fp),      intent(in)    :: lowerb(3)
  real(kind=fp),      intent(in)    :: upperb(3)
  integer,            intent(out)   :: stat
  real(kind=fp),      intent(inout) :: tuv(3)
  real(kind=fp),      intent(out)   :: toltuv
  real(kind=fp),      intent(out)   :: xyz(3)
  real(kind=fp)                     :: shiftv(3), sgn
  real(kind=fp), dimension(3)       :: xyz_c, xyz_s, r
  real(kind=fp), dimension(3)       :: ct, ctt, su, sv, suu, suv, svv, n
  real(kind=fp)                     :: resxyz
  real(kind=fp)                     :: grad(3), hess(3,3), dtuv(3), cond, errtuv
  integer                           :: it
  
  call evald2(ctt, curv, tuv(1))

  stat = 1
  errtuv = 0._fp
  cond = 1._fp
  do it = 1,itmax
     call evald1(ct, curv, tuv(1))
     call evald1(su, surf, tuv(2:3), 1)
     call evald1(sv, surf, tuv(2:3), 2)
     n = cross(su, sv)
     
     sgn = sign(1._fp, dot_product(n, ctt))
     shiftv = sgn * shiftd * n
     
     !! compute residual
     call eval(xyz_c, curv, tuv(1)  ) ! curve's position vector
     call eval(xyz_s, surf, tuv(2:3)) ! surface's position vector
     r = xyz_s - xyz_c
     resxyz = sum(r**2)
     r = r - shiftv
     IF ( DEBUG ) PRINT *,SQRT(RESXYZ), SQRT(ERRTUV), EPSFP*COND
     
     !! compute gradient of distance
     grad(1) = -dot_product(ct, r)
     grad(2) =  dot_product(su, r)
     grad(3) =  dot_product(sv, r)

     !! compute Hessian matrix
     if ( it > 1 ) call evald2(ctt, curv, tuv(1))
     call evald2(suu, surf, tuv(2:3), 1)
     call evald2(suv, surf, tuv(2:3), 2)
     call evald2(svv, surf, tuv(2:3), 3)

     hess(1,1) = -dot_product(ctt,r) + dot_product(ct,ct)
     hess(1,2) = -dot_product(ct,su)
     hess(1,3) = -dot_product(ct,sv)
     hess(2,2) =  dot_product(suu,r) + dot_product(su,su)
     hess(2,3) =  dot_product(suv,r) + dot_product(su,sv)
     hess(3,3) =  dot_product(svv,r) + dot_product(sv,sv)
     hess(2,1) = hess(1,2) ! ...
     hess(3,1) = hess(1,3) ! symmetry
     hess(3,2) = hess(2,3) ! ...

     !! solve for Newton step
     call linsolve_svd( &
          dtuv, &
          hess, &
          -grad, &
          3, &
          3, &
          1, &
          cond )

     errtuv = sum(dtuv**2)

     ! correct Newton step to keep the iterate inside feasible region
     call nd_box_reflexions( &
          tuv, &
          lowerb, &
          upperb, &
          dtuv, &
          3 )

     ! update solution
     tuv = tuv + dtuv

     !! termination criterion
     if ( errtuv < max(EPSuvsqr, EPSfpsqr*cond**2) ) then
        if ( resxyz < EPSxyzsqr ) then
           if ( errtuv > EPSuvsqr ) then
              PRINT *,'newton_curve_surface_dynshift : /!\ toltuv > EPSuv'
           end if
           ! converged to a curve-surface intersection point
           stat = 0
           toltuv = max(EPSmath, sqrt(errtuv))
           call eval(xyz_c, curv, tuv(1)  )
           call eval(xyz_s, surf, tuv(2:3))
           xyz = 0.5_fp * (xyz_c + xyz_s)
           IF ( DEBUG ) PRINT *,SQRT(RESXYZ), SQRT(ERRTUV), EPSFP*COND
           IF ( DEBUG ) PRINT *,'CONVERGED, TUV =',TUV,', XYZ =',XYZ,', TOL =',TOLTUV
           !PAUSE
        else
           IF ( DEBUG ) PRINT *,'STAGNATION'
        end if
        return
     end if
  end do
end subroutine newton_curve_surface_dynshift
