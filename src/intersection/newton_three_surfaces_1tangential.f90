subroutine newton_three_surfaces_1tangential( &
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
  USE MOD_UTIL
  USE MOD_POLYNOMIAL
  implicit none
  LOGICAL, PARAMETER :: DEBUG = .false.!( GLOBALDEBUG .AND. .true. )
  integer,           parameter     :: itmax = 2 + ceiling(-log10(EPSuv))
  type(ptr_surface), intent(in)    :: surf(3)
  real(kind=fp),     intent(in)    :: lowerb(6)
  real(kind=fp),     intent(in)    :: upperb(6)
  integer,           intent(out)   :: stat
  real(kind=fp),     intent(inout) :: uv(2,3)
  real(kind=fp),     intent(out)   :: xyz(3)
  real(kind=fp)                    :: xyzs(3,3), r(6)
  type(ptr_surface)                :: surftng(2)
  real(kind=fp)                    :: uvtng(2,2)
  real(kind=fp)                    :: duv_ds(2,2,2), dxyz_ds(3,2), curvature(2)
  real(kind=fp)                    :: mat(3,3), dtuv(3), duv(6)
  integer                          :: stat_contactpoint
  real(kind=fp)                    :: cond, erruv
  integer                          :: it, isurf, ivar
  integer                          :: stat_refl
  
  IF ( DEBUG ) THEN
     PRINT *,''
     PRINT *,'NEWTON_THREE_SURFACES_1TANGENTIAL'
     PRINT *,'UV0 ='
     CALL PRINT_MAT(transpose(UV))
  END IF
  
  stat = 1
  erruv = 0._fp
  cond = 1._fp

  surftng(1)%ptr => surf(1)%ptr
  surftng(2)%ptr => surf(3)%ptr

  IF ( DEBUG ) PRINT *,'|X1 - X2| , |X1 - X3| , |DUV| , eps*cond'
  do it = 1,itmax
     IF ( DEBUG ) PRINT *,'IT. #',IT
     uvtng = uv(:,[1,3])
     !! relax tangential intersection
     call simultaneous_point_inversions( &
          surftng, &
          lowerb([1,2,5,6]) - 100._fp*EPSuv, &
          upperb([1,2,5,6]) + 100._fp*EPSuv, &
          stat, &
          uvtng, &
          xyz )
     if ( stat > 0 ) then
        CALL WRITE_POLYNOMIAL(surf(1)%ptr%x, '../debug/debug_n3s1t_surf1.cheb')
        CALL WRITE_POLYNOMIAL(surf(2)%ptr%x, '../debug/debug_n3s1t_surf2.cheb')
        CALL WRITE_POLYNOMIAL(surf(3)%ptr%x, '../debug/debug_n3s1t_surf3.cheb')
        CALL write_matrix(uv, 2, 3, '../debug/debug_n3s1t_uv.dat')
        STOP 'newton_three_surfaces_1tangential : failed to relax onto tangential intersection :('
     end if
     uv(:,[1,3]) = uvtng
     IF ( DEBUG ) THEN
        PRINT *,'UV='
        CALL PRINT_MAT(transpose(UV))
     END IF
     
     !! compute residual
     do isurf = 1,3
        call eval( &
             xyzs(:,isurf), &
             surf(isurf)%ptr, &
             uv(:,isurf) )
     end do
     IF ( DEBUG ) THEN
        PRINT *,'XYZ='
        CALL PRINT_MAT(transpose(XYZS))
     END IF
     
     r(1:3) = xyzs(:,1) - xyzs(:,2)
     r(4:6) = xyzs(:,1) - xyzs(:,3)
     IF ( .true. ) THEN
        if ( max(sum(r(1:3)**2), sum(r(4:6)**2)) < EPSxyzsqr ) then
           stat = 0
           xyz = sum(xyzs, 2)/3._fp
           IF ( DEBUG ) PRINT *,'CONVERGED, UV =',UV,', XYZ =',XYZ
           return
        end if
     END IF

     !! compute derivatives
     if ( sum(r(4:6)**2) > EPSxyzsqr ) then
        PRINT *,'|X1 - X3| = ',norm2(r(4:6))
        STOP 'newton_three_surfaces_1tangential : not an intersection point :('
     else
        ! tangent to tangential intersection curve
        call diffgeom_intersection( &
             surftng, &
             uv(:,[1,3]), &
             duv_ds, &
             dxyz_ds, &
             stat_contactpoint, &
             curvature )
        !IF ( DEBUG ) PRINT *,'STAT_CONTACTPOINT =',stat_contactpoint
        IF ( DEBUG ) PRINT *,'DUV_DS  =',duv_ds(:,1,:)
        !IF ( DEBUG ) PRINT *,'DXYZ_DS =',dxyz_ds(:,1)
        if ( stat_contactpoint /= 1 ) then
           PRINT *,'STAT_CONTACTPOINT =',STAT_CONTACTPOINT
           STOP 'newton_three_surfaces_1tangential : not a tangential intersection curve'
        end if
     end if
     mat(:,1) = -dxyz_ds(:,1)

     ! tangents to surface
     do ivar = 1,2
        call evald1( &
             mat(:,ivar+1), &
             surf(2)%ptr, &
             uv(:,2), &
             ivar )
     end do

     !! solve for Newton step
     IF ( DEBUG ) PRINT *,'MAT ='
     IF ( DEBUG ) CALL PRINT_MAT(MAT)
     IF ( DEBUG ) PRINT *,'RHS =',0.5_fp*(xyzs(:,1) + xyzs(:,3)) - xyzs(:,2)
     call linsolve_svd( &
          dtuv, &
          mat, &
          0.5_fp*(xyzs(:,1) + xyzs(:,3)) - xyzs(:,2), &
          3, &
          3, &
          1, &
          cond, &
          tol=EPSmath )
     PRINT *,'DTUV =',DTUV

     duv(1:2) = dtuv(1)*duv_ds(:,1,1)
     duv(3:4) = dtuv(2:3)
     duv(5:6) = dtuv(1)*duv_ds(:,1,2)
     
     erruv = sum(duv**2) / 3._fp

     !! correct Newton step to keep the iterate inside feasible region
     call nd_box_reflexions( &
          reshape(uv, [6]), &
          lowerb, &
          upperb, &
          duv, &
          6, &
          stat_refl )
     if ( stat_refl > 0 ) return
     !IF ( DEBUG ) PRINT *,'DUV =',DUV
     IF ( DEBUG ) THEN
        PRINT *,'DUV ='
        DO IVAR = 1,SIZE(DUV)
           PRINT *,DUV(IVAR)
        END DO
     END IF

     !! update solution
     uv(:,1) = uv(:,1) + duv(1:2)
     uv(:,2) = uv(:,2) + duv(3:4)
     uv(:,3) = uv(:,3) + duv(5:6)

     IF ( DEBUG ) PRINT *,norm2(r(1:3)), norm2(r(4:6)), sqrt(erruv), EPSfp*cond

     !! termination criteria
     if ( erruv < max(EPSuvsqr, EPSfpsqr*cond**2) ) then
        if ( max(sum(r(1:3)**2), sum(r(4:6)**2)) < EPSxyzsqr ) then
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


end subroutine newton_three_surfaces_1tangential
