subroutine newton_three_surfaces( &
     surf, &
     lowerb, &
     upperb, &
     stat, &
     uv, &
     xyz)
  use mod_math
  use mod_linalg
  use mod_diffgeom
  use mod_tolerances  
  implicit none
  LOGICAL, PARAMETER :: DEBUG = .false.!( GLOBALDEBUG .AND. .true. )
  integer,           parameter     :: itmax = 2 + ceiling(-log10(EPSuv))
  type(ptr_surface), intent(in)    :: surf(3)
  real(kind=fp),     intent(in)    :: lowerb(6)
  real(kind=fp),     intent(in)    :: upperb(6)
  integer,           intent(out)   :: stat
  real(kind=fp),     intent(inout) :: uv(2,3)
  real(kind=fp),     intent(out)   :: xyz(3)
  real(kind=fp)                    :: xyzs(3,3)
  real(kind=fp)                    :: r(6), jac(6,6), duv(6)
  real(kind=fp)                    :: cond, erruv
  integer                          :: it, isurf

  IF ( DEBUG ) THEN
     PRINT *,''; PRINT *,'';
     PRINT *,'NEWTON_THREE_SURFACES'
     PRINT *,'LOWERB =',LOWERB
     PRINT *,'UPPERB =',UPPERB
  END IF
  stat = 1
  erruv = 0._fp
  cond = 1._fp

  do it = 1,itmax
     !IF ( DEBUG ) PRINT *,'UV =',UV
     !! compute residual
     do isurf = 1,3
        call eval( &
             xyzs(:,isurf), &
             surf(isurf)%ptr, &
             uv(:,isurf) )
     end do
     r(1:3) = xyzs(:,1) - xyzs(:,2)
     r(4:6) = xyzs(:,1) - xyzs(:,3)

     !! compute Jacobian matrix
     call evald1(jac(1:3,1), surf(1)%ptr, uv(:,1), 1)
     call evald1(jac(1:3,2), surf(1)%ptr, uv(:,1), 2)
     call evald1(jac(1:3,3), surf(2)%ptr, uv(:,2), 1)
     call evald1(jac(1:3,4), surf(2)%ptr, uv(:,2), 2)
     call evald1(jac(4:6,5), surf(3)%ptr, uv(:,3), 1)
     call evald1(jac(4:6,6), surf(3)%ptr, uv(:,3), 2)
     jac(4:6,1:2) = jac(1:3,1:2)
     jac(1:6,3:6) = -jac(1:6,3:6)
     jac(1:3,5:6) = 0._fp
     jac(4:6,3:4) = 0._fp

     IF ( .false. ) THEN
        PRINT *,'JAC ='
        CALL PRINT_MAT(JAC)

        PRINT *,'RHS =',-r
     END IF

     !! solve for Newton step
     call linsolve_svd( &
          duv, &
          jac, &
          -r, &
          6, &
          6, &
          1, &
          cond, &
          tol=0.1_fp*EPSuv )

     IF ( COND > EPSUV/EPSFP ) THEN
        if ( max(sum(r(1:3)**2), sum(r(4:6)**2)) < EPSxyzsqr ) then
           stat = 0
           erruv = 10._fp*EPSuv !*******
           IF ( DEBUG ) PRINT *,'newton_three_surfaces : /!\ singular Jacobian at solution'
           xyz = sum(xyzs, 2)/3._fp
           return
        end if
     END IF

     erruv = maxval([sum(duv(1:2)**2), sum(duv(3:4)**2), sum(duv(5:6))**2])

     !! correct Newton step to keep the iterate inside feasible region
     call nd_box_reflexions( &
          reshape(uv, [6]), &
          lowerb, &
          upperb, &
          duv, &
          6 )

     !! update solution
     uv(:,1) = uv(:,1) + duv(1:2)
     uv(:,2) = uv(:,2) + duv(3:4)
     uv(:,3) = uv(:,3) + duv(5:6)

     IF ( DEBUG ) PRINT *,norm2(r(1:3)), norm2(r(4:6)), sqrt(erruv), EPSfp*cond
     !! termination criteria
     if ( erruv < max(EPSuvsqr, EPSfpsqr*cond**2) ) then
        if ( max(sum(r(1:3)**2), sum(r(4:6)**2)) < EPSxyzsqr ) then
           if ( erruv > EPSuvsqr ) then
              IF ( DEBUG ) PRINT *,'newton_three_surfaces : /!\ toluv > EPSuv'
              pause
           end if
           stat = 0
           do isurf = 1,3
              call eval( &
                   xyzs(:,isurf), &
                   surf(isurf)%ptr, &
                   uv(:,isurf) )
           end do
           xyz = sum(xyzs, 2)/3._fp
           IF ( DEBUG ) PRINT *,'CONVERGED, UV =',UV,', XYZ =',XYZ
        else
           IF ( DEBUG ) PRINT *,'STAGNATION'
        end if
        return
     end if

  end do

end subroutine newton_three_surfaces
