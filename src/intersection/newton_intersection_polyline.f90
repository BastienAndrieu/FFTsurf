subroutine newton_intersection_polyline( &
     surf, &
     lowerb, &
     upperb, &
     xyzp_prev, &
     htargetsqr, &
     stat, &
     uv, &
     xyzp )
  use mod_math
  use mod_linalg
  use mod_diffgeom
  use mod_tolerances  
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  integer,           parameter     :: itmax = 2 + ceiling(-log10(EPSuv))
  type(ptr_surface), intent(in)    :: surf(2)
  real(kind=fp),     intent(in)    :: lowerb(4)
  real(kind=fp),     intent(in)    :: upperb(4)
  real(kind=fp),     intent(in)    :: xyzp_prev(3)
  real(kind=fp),     intent(in)    :: htargetsqr
  integer,           intent(out)   :: stat
  real(kind=fp),     intent(inout) :: uv(2,2)
  real(kind=fp),     intent(out)   :: xyzp(3)
  real(kind=fp)                    :: xyz(3,2), resh
  real(kind=fp)                    :: r1(3), r2(3), jac(4,4), duv(4)
  real(kind=fp)                    :: cond, erruv
  integer                          :: it, isurf, ivar

  IF ( DEBUG ) THEN
     PRINT *,''; PRINT *,'';
     PRINT *,'NEWTON_INTERSECTION_POLYLINE'
     PRINT *,'HTARGETSQR =',htargetsqr
     PRINT *,'LOWERB =',LOWERB
     PRINT *,'UPPERB =',UPPERB
  END IF

  stat = 1
  erruv = 0._fp
  cond = 1._fp

  do it = 1,itmax
     IF ( DEBUG ) PRINT *,'NEWTON_INTERSECTION_POLYLINE, IT#',IT
     IF ( DEBUG ) PRINT *,'UV =',UV
     !! compute residual
     do isurf = 1,2
        call eval( &
             xyz(:,isurf), &
             surf(isurf)%ptr, &
             uv(:,isurf) )
     end do
     r1 = xyz(:,1) - xyz(:,2)
     r2 = xyz(:,1) - xyzp_prev
     resh = sum(r2**2) - htargetsqr

     !! compute Jacobian matrix
     do isurf = 1,2
        do ivar = 1,2
           call evald1( &
                jac(1:3,2*(isurf-1)+ivar), &
                surf(isurf)%ptr, &
                uv(:,isurf), &
                ivar)
        end do
     end do
     jac(1:3,3:4) = -jac(1:3,3:4)
     do ivar = 1,2
        jac(4,ivar) = 2._fp*dot_product(jac(1:3,ivar), r2)
     end do
     jac(4,3:4) = 0._fp

     !! solve for Newton step
     !PRINT *,'RHS=',-R1,-RESH
     !PRINT *,'JAC='
     !CALL PRINT_MAT(JAC)
     call linsolve_svd( &
          duv, &
          jac, &
          -[r1,resh], &
          4, &
          4, &
          1, &
          cond )
     erruv = max(sum(duv(1:2)**2), sum(duv(3:4)**2))

     !PRINT *,'DUV =',DUV

     !! correct Newton step to keep the iterate inside feasible region
     call nd_box_reflexions( &
          reshape(uv, [4]), &
          lowerb, &
          upperb, &
          duv, &
          4 )
     !PRINT *,'DUV* =',DUV
     
     !! update solution
     uv(:,1) = uv(:,1) + duv(1:2)
     uv(:,2) = uv(:,2) + duv(3:4)
     
     !! termination criteria
     IF ( DEBUG ) PRINT *,norm2(r1), sqrt(abs(resh)), sqrt(erruv), EPSfp*cond
     if ( erruv < max(EPSuvsqr, EPSfpsqr*cond**2) ) then ! <--------+
        if ( sum(r1**2) < EPSxyzsqr .and. &                         !
             abs(resh)  < tolhsqr ) then ! <---------------------+  !
           stat = 0                                              !  !
           do isurf = 1,2 ! <------------+                       !  !
              call eval( &               !                       !  !
                   xyz(:,isurf), &       !                       !  !
                   surf(isurf)%ptr, &    !                       !  !
                   uv(:,isurf) )         !                       !  !
           end do ! <--------------------+                       !  !
           xyzp = 0.5_fp * sum(xyz, 2)                           !  !
           IF ( DEBUG ) PRINT *,'CONVERGED, UV =',UV,', XYZ =',XYZP
        else ! --------------------------------------------------+  !
           IF ( DEBUG ) PRINT *,'STAGNATION'                     !  !
        end if ! ------------------------------------------------+  !
        return                                                      !
     end if ! <-----------------------------------------------------+

  end do

end subroutine newton_intersection_polyline
