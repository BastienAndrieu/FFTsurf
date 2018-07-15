subroutine newton_tracing( &
     surf, &
     lowerb, &
     upperb, &
     p, &
     wtarget, &
     tolw, &
     stat, &
     uv, &
     xyzp )
  use mod_math
  use mod_linalg
  use mod_diffgeom2
  use mod_tolerances  
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  integer,           parameter     :: itmax = 2 + ceiling(-log10(EPSuv))
  type(ptr_surface), intent(in)    :: surf(2)
  real(kind=fp),     intent(in)    :: lowerb(4)
  real(kind=fp),     intent(in)    :: upperb(4)
  real(kind=fp),     intent(in)    :: p(3)
  real(kind=fp),     intent(in)    :: wtarget
  real(kind=fp),     intent(in)    :: tolw
  integer,           intent(out)   :: stat
  real(kind=fp),     intent(inout) :: uv(2,2)
  real(kind=fp),     intent(out)   :: xyzp(3)
  real(kind=fp)                    :: xyz(3,2), dxyz_duv(3,2,2)
  real(kind=fp)                    :: resxyz, resw
  real(kind=fp)                    :: r(3), jac(4,4), duv(4)
  integer                          :: rank
  real(kind=fp)                    :: cond, erruv
  integer                          :: it, isurf, ivar
  
  !IF ( DEBUG ) THEN
  !   PRINT *,''; PRINT *,'';
  !   PRINT *,'NEWTON_INTERSECTION_POLYLINE'
  !   PRINT *,'LOWERB =',LOWERB
  !   PRINT *,'UPPERB =',UPPERB
  !END IF
  stat = 1
  erruv = 0._fp
  cond = 1._fp

  do it = 1,itmax
     !! compute residual
     do isurf = 1,2
        call eval( &
             xyz(:,isurf), &
             surf(isurf)%ptr, &
             uv(:,isurf) )
     end do
     r = xyz(:,1) - xyz(:,2)

     !! termination criteria
     resxyz = sum(r**2)
     resw   = dot_product(p, xyz(:,1)) - wtarget
     IF ( DEBUG ) PRINT *,sqrt(resxyz), abs(resw), sqrt(erruv), EPSuv*cond

     if ( erruv < EPSuvsqr*cond**2 .and. it > 1 ) then
        if ( resxyz < EPSxyzsqr .and. abs(resw) < tolw ) then
           stat = 0
           xyzp = 0.5_fp * sum(xyz, 2)
           IF ( DEBUG ) PRINT *,'CONVERGED, UV =',UV,', XYZ =',XYZP
        else
           IF ( DEBUG ) PRINT *,'STAGNATION'
        end if
        return
     end if

     !! compute Jacobian matrix
     do isurf = 1,2
        do ivar = 1,2
           call evald1( &
                dxyz_duv(:,ivar,isurf), &
                surf(isurf)%ptr, &
                uv(:,isurf), &
                ivar )
        end do
        jac(1:3,2*isurf-1:2*isurf) = real((-1)**(isurf+1), kind=fp) * dxyz_duv(1:3,1:2,isurf)
        if ( isurf == 1 ) then
           do ivar = 1,2
              jac(4,ivar) = dot_product(p, dxyz_duv(1:3,ivar,1))
           end do
        end if
     end do
     jac(4,3:4) = 0._fp

     !! solve for Newton step
     call linsolve_svd( &
          duv, &
          jac, &
          -[r,resw], &
          4, &
          4, &
          1, &
          cond, &
          rank )

     erruv = sum( duv**2 )

     if ( rank < 4 ) then ! <------------------+
        ! singular Jacobian matrix             !
        if ( stat == 0 ) then ! <----+         !
           stat = 2                  !         !
           IF ( DEBUG ) PRINT *,'SINGULAR JACOBIAN AT SOLUTION'
           return                    !         !
        end if ! <-------------------+         !
     end if ! <--------------------------------+

     !! correct Newton step to keep the iterate inside feasible region
     call nd_box_reflexions( &
          reshape(uv, [4]), &
          lowerb, &
          upperb, &
          duv, &
          4 )

     !! update solution
     uv(:,1) = uv(:,1) + duv(1:2)
     uv(:,2) = uv(:,2) + duv(3:4)

  end do


end subroutine newton_tracing
