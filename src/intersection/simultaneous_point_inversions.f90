subroutine simultaneous_point_inversions( &
     surf, &
     lowerb, &
     upperb, &
     stat, &
     uv, &
     xyz )
  use mod_math
  use mod_linalg
  use mod_diffgeom
  use mod_tolerances
  implicit none
  LOGICAL, PARAMETER :: DEBUG = .false. !( GLOBALDEBUG .AND. .true. )
  integer,           parameter     :: itmax = 2 + ceiling(-log10(EPSuv))
  type(ptr_surface), intent(in)    :: surf(2)
  real(kind=fp),     intent(in)    :: lowerb(4)
  real(kind=fp),     intent(in)    :: upperb(4)
  integer,           intent(out)   :: stat
  real(kind=fp),     intent(inout) :: uv(2,2)
  real(kind=fp),     intent(out)   :: xyz(3)
  real(kind=fp)                    :: xyzs(3,2), r(3)
  real(kind=fp)                    :: dxyz_duv(3,2), mat(2,2), rhs(2), duv(2)
  real(kind=fp), dimension(2)      :: cond, erruv
  integer                          :: it, isurf, ivar
  integer                          :: stat_refl
  
  !do isurf = 1,2
  !   uv(1:2,isurf) = max(lowerb(2*isurf-1:2*isurf), min(upperb(2*isurf-1:2*isurf), uv(1:2,isurf)))
  !end do
  
  IF ( DEBUG ) THEN
     PRINT *,''
     PRINT *,'SIMULTANEOUS_POINT_INVERSION'
     PRINT *,'UV =',uv
  END IF
  
  stat = 1
  erruv(1:2) = 0._fp
  cond(1:2) = 1._fp

  IF ( DEBUG ) PRINT *,'|X1 - X2| , |DUV1| , |DUV2| , eps*cond1 , eps*cond2'
  do it = 1,itmax ! <-------------------------------------------------------+
     !IF ( DEBUG ) PRINT *,'IT. #',IT
     !IF ( DEBUG ) THEN
     !   PRINT *,'UV='
     !   CALL PRINT_MAT(transpose(UV))
     !END IF
     do isurf = 1,2 ! <--------------+                                      !
        call eval( &                 !                                      !
             xyzs(:,isurf), &        !                                      !
             surf(isurf)%ptr, &      !                                      !
             uv(:,isurf) )           !                                      !
        !PRINT *,'xyz =',xyzs(:,isurf)
     end do ! <----------------------+                                      !
     r = 0.5_fp * (xyzs(:,1) - xyzs(:,2))                                   !
     !                                                                      !
     do isurf = 1,2 ! <--------------------------------------------------+  !
        do ivar = 1,2 ! <---------------+                                !  !
           call evald1( &               !                                !  !
                dxyz_duv(:,ivar), &     !                                !  !
                surf(isurf)%ptr, &      !                                !  !
                uv(:,isurf), &          !                                !  !
                ivar )                  !                                !  !
        end do ! <----------------------+                                !  !
        mat = matmul(transpose(dxyz_duv), dxyz_duv)                      !  !
        rhs = matmul(transpose(dxyz_duv), real((-1)**isurf, kind=fp)*r)  !  !
        call linsolve_svd( &                                             !  !
             duv, &                                                      !  !
             mat, &                                                      !  !
             rhs, &                                                      !  !
             2, &                                                        !  !
             2, &                                                        !  !
             1, &                                                        !  !
             cond(isurf), &                                              !  !
             tol=EPSmath )                                               !  !
        !                                                                !  !
        erruv(isurf) = sum(duv**2)                                       !  !
        !                                                                !  !
        ! modify Newton step => keep iterate inside feasible region      !  !
        call nd_box_reflexions( &                                        !  !
             uv(:,isurf), &                                              !  !
             lowerb(2*isurf-1:2*isurf), &                                !  !
             upperb(2*isurf-1:2*isurf), &                                !  !
             duv, &                                                      !  !
             2, &                                                        !  !
             stat_refl )                                                 !  !
        if ( stat_refl > 0 ) return                                      !  !
        !                                                                !  !
        ! update solution                                                !  !
        uv(:,isurf) = uv(:,isurf) + duv                                  !  !
     end do ! <----------------------------------------------------------+  !
     !                                                                      !
     IF ( DEBUG ) PRINT *,norm2(r), sqrt(erruv), EPSfp*cond
     ! termination criteria                                                 !
     if ( all(erruv < max(EPSuvsqr, EPSfpsqr*cond**2)) ) then ! <--------+  !
        if ( sum(r**2) < 0.25_fp*EPSxyzsqr ) then ! <----+               !  !
           stat = 0                                      !               !  !
           do isurf = 1,2 ! <--------------+             !               !  !
              call eval( &                 !             !               !  !
                   xyzs(:,isurf), &        !             !               !  !
                   surf(isurf)%ptr, &      !             !               !  !
                   uv(:,isurf) )           !             !               !  !
           end do ! <----------------------+             !               !  !
           xyz = 0.5_fp * sum(xyzs, 2)                   !               !  !
           IF ( DEBUG ) PRINT *,'CONVERGED, UV =',UV,', XYZ =',XYZ
        else ! ------------------------------------------+               !  !
           IF ( DEBUG ) PRINT *,'STAGNATION'
        end if ! <---------------------------------------+               !  !
        return                                                           !  !
     end if ! <----------------------------------------------------------+  !
  end do ! <----------------------------------------------------------------+

  if ( stat > 0 ) then
     if ( all(erruv < max(100._fp * EPSuvsqr, EPSfpsqr*cond**2)) ) then
        if ( sum(r**2) < 25_fp*EPSxyzsqr ) then
           stat = 0
           do isurf = 1,2 ! <--------------+
              call eval( &                 !
                   xyzs(:,isurf), &        !
                   surf(isurf)%ptr, &      !
                   uv(:,isurf) )           !
           end do ! <----------------------+
           xyz = 0.5_fp * sum(xyzs, 2)
           IF ( DEBUG ) PRINT *,'CONVERGED*, UV =',UV,', XYZ =',XYZ
        end if
     end if
  end if

end subroutine simultaneous_point_inversions
