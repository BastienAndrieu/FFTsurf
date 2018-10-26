subroutine trace_intersection_polyline( &
     surf, &
     uvbox, &
     param_vector, &
     uv_endpoints, &
     xyz_endpoints, &
     stat, &
     polyline, &
     w0, &
     tolchord, &
     hmin, &
     hmax )
  use mod_math
  use mod_diffgeom
  use mod_types_intersection
  use mod_tolerances
  ! Returns: stat =-1 if one surface is singular at a polyline point (normal = zero vector)
  !               = 0 if everything went OK :)
  !               = 1 if an error ocurred when reallocating polyline%uv
  !               = 2 if an error ocurred when reallocating polyline%xyz
  !               = 3 if an isolated tangential contact point was encountered
  !               = 4 if an high-order tangential contact point was encountered
  !               = 5 if backtracking did not terminate successfully
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  real(kind=fp), parameter                        :: FRACbacktrack = 0.5_fp
  real(kind=fp), parameter                        :: EPSbacktrack = real(1d-2, kind=fp)
  type(ptr_surface),                intent(in)    :: surf(2)
  real(kind=fp),                    intent(in)    :: uvbox(2,2,2)
  real(kind=fp),                    intent(inout) :: param_vector(3)
  real(kind=fp),                    intent(in)    :: uv_endpoints(2,2,2) ! u/v, #surf, #point
  real(kind=fp),                    intent(in)    :: xyz_endpoints(3,2)
  integer,                          intent(out)   :: stat
  type(type_intersection_polyline), intent(inout) :: polyline
  real(kind=fp),                    intent(out)   :: w0
  real(kind=fp),                    intent(in)    :: tolchord
  real(kind=fp),                    intent(in)    :: hmin, hmax
  real(kind=fp)                                   :: FRACcurvature_radius 
  integer                                         :: stat_tangent, stat_insertion, stat_newton
  real(kind=fp), dimension(4)                     :: lowerb, upperb
  real(kind=fp)                                   :: Dw, w, wprev, dist_from_end
  real(kind=fp)                                   :: duv_ds(2,2,2), dxyz_ds(3,2), curvature(2)
  real(kind=fp)                                   :: h_endpoints(2), h, EPSh, hnext
  real(kind=fp)                                   :: uv(2,2), xyz(3), lambda
  integer                                         :: ipt

  FRACcurvature_radius = 2._fp*sqrt(tolchord*(2._fp - tolchord))
  
  stat = 0
  lowerb = reshape(uvbox(1,1:2,1:2), [4]) - EPSuv
  upperb = reshape(uvbox(2,1:2,1:2), [4]) + EPSuv
  IF ( DEBUG ) THEN
     PRINT *,''; PRINT *,'';
     PRINT *,'TRACE_INTERSECTION_POLYLINE'
     PRINT *,' PARAM_VECTOR =', PARAM_VECTOR
     PRINT *,' UV_ENDPOINTS ='; PRINT *,uv_endpoints(:,:,1); PRINT *,uv_endpoints(:,:,2)
     PRINT *,'XYZ_ENDPOINTS ='; PRINT *,xyz_endpoints(:,1);  PRINT *,xyz_endpoints(:,2)
     PRINT *,'LOWERB =',LOWERB; PRINT *,'UPPERB =',UPPERB
  END IF
  w0 = dot_product(param_vector, xyz_endpoints(:,1))
  Dw = dot_product(param_vector, xyz_endpoints(:,2)) - w0
  param_vector = param_vector/Dw
  w0 = w0/Dw

  !PRINT *,'CURVATURE RADIUS AT ENDPOINTS ='
  do ipt = 2,1,-1
     ! compute tangent direction and curvature of the intersection curve at endpoints
     ! (in reverse order, so data at first endpoint is kept in memory)
     call diffgeom_intersection( &
          surf, &
          uv_endpoints(:,:,ipt), &
          duv_ds, &
          dxyz_ds, &
          stat_tangent, &
          curvature )
     !call diffgeom_intersection_curve( &
     !     surf, &
     !     uv_endpoints(:,:,ipt), &
     !     duv_ds, &
     !     dxyz_ds, &
     !     stat_tangent, &
     !     curvature )
     IF ( DEBUG ) PRINT *,'CURVATURE RADIUS = ',1._FP/CURVATURE(1)
     if ( stat_tangent < 0 .or. stat_tangent > 2 ) then ! <----+
        IF ( DEBUG ) PRINT *,'stat_tangent = ',stat_tangent    !
        stat = stat_tangent                                    !
        return                                                 !
     end if ! <------------------------------------------------+
     curvature(1) = max(EPSfp, curvature(1))
     h_endpoints(ipt) = FRACcurvature_radius / curvature(1)
  end do
  h_endpoints = max(h_endpoints, hmin) 
  h_endpoints = min(h_endpoints, hmax)

  
  ! first point
  call insert_polyline_point( &
       uv_endpoints(1:2,1:2,1), &
       xyz_endpoints(1:3,1), &
       stat_insertion, &
       polyline, &
       i=0 )
  if ( stat_insertion > 0 ) then ! <--------+
     stat = stat_insertion                  !
     return                                 !
  end if ! <--------------------------------+

  wprev = 0._fp
  h = h_endpoints(1)
  outer : do ! <-----------------------------------------------------------------+
     IF ( DEBUG ) PRINT *,'HTARGET =',H
     ! check if the current is close enough to the end of the polyline           !
     dist_from_end = sum((polyline%xyz(:,polyline%np) - xyz_endpoints(:,2))**2)  !
     if ( dist_from_end < h**2 ) then ! <-----------------------+                !
        if ( dist_from_end < h_endpoints(2)**2 ) then ! <----+  !                !
           if ( dist_from_end < (TOLh * min(h, h_endpoints(2)))**2 ) polyline%np = polyline%np - 1
           exit outer                                        !  !                !
        else ! ----------------------------------------------+  !                !
           IF ( DEBUG ) PRINT *,'SHORTEN H : H =',h,', DIST =',sqrt(dist_from_end),', HEND =',h_endpoints(2)
           h = sqrt(dist_from_end)                           !  !                !
        end if ! <-------------------------------------------+  !                !
     end if ! <-------------------------------------------------+                !
     !                                                                           !
     ! compute next point                                                        !
     EPSh = EPSbacktrack * h                                                     !
     inner : do ! <-----------------------------------------------------------+  !
        ! set initial iterate                                                 !  !
        call nd_box_constraint( &                                             !  !
             reshape(polyline%uv(:,:,polyline%np), [4]), &                    !  !
             lowerb, &                                                        !  !
             upperb, &                                                        !  !
             h*reshape(duv_ds(:,1,:), [4]), &                                 !  !
             lambda )                                                         !  !
        !                                                                     !  !
        IF (DEBUG .AND. LAMBDA < 1._FP) THEN
           PRINT *,' UV =',polyline%uv(:,:,polyline%np)
           PRINT *,'DUV =',h*duv_ds(:,1,:)
           PRINT *,'LAMBDA =',LAMBDA
        END IF
        uv = polyline%uv(:,:,polyline%np) + lambda*h*duv_ds(:,1,:)            !  !
        !                                                                     !  !
        ! refine using Newton-Raphson algorithm                               !  !
        IF ( DEBUG ) PRINT *,'UV0 =',UV
        call newton_intersection_polyline( &                                  !  !
             surf, &                                                          !  !
             lowerb, &                                                        !  !
             upperb, &                                                        !  !
             polyline%xyz(1:3,polyline%np), &                                 !  !
             (lambda*h)**2, &                                                 !  !
             stat_newton, &                                                   !  !
             uv, &                                                            !  !
             xyz )                                                            !  !
        !                                                                     !  !
        IF ( DEBUG ) PRINT *,'STAT_NEWTON =',STAT_NEWTON
        if ( stat_newton == 0 ) then ! <----------------------------------+   !  !
           ! Newton has converged, check if w is monotonic                !   !  !
           w = dot_product(param_vector, xyz) - w0                        !   !  !
           if ( is_in_open_interval(w, wprev, 1._fp) ) then ! <--------+  !   !  !
              ! compute tangent direction and curvature                !  !   !  !
              ! of the intersection curve at the current point         !  !   !  !
              call diffgeom_intersection( &                            !  !   !  !
                   surf, &                                             !  !   !  !
                   uv, &                                               !  !   !  !
                   duv_ds, &                                           !  !   !  !
                   dxyz_ds, &                                          !  !   !  !
                   stat_tangent, &                                     !  !   !  !
                   curvature )                                         !  !   !  !
              IF ( DEBUG ) PRINT *,'CURVATURE RADIUS = ',1._FP/CURVATURE(1)
              if ( stat_tangent < 0 ) then ! <--------------------+    !  !   !  !
                 stat = stat_tangent                              !    !  !   !  !
                 return ! singular surface (undefined normal)     !    !  !   !  !
              elseif ( stat > 2 ) then ! -------------------------+    !  !   !  !
                 stat = stat_tangent                              !    !  !   !  !
                 return ! undefined tangent to intersection curve !    !  !   !  !
              end if ! <------------------------------------------+    !  !   !  !
              !                                                        !  !   !  !
              curvature(1) = max(EPSfp, curvature(1))                  !  !   !  !
              hnext = FRACcurvature_radius / curvature(1)              !  !   !  !
              hnext = max(hnext, hmin)                                 !  !   !  !
              hnext = min(hnext, hmax)                                 !  !   !  !
              if ( lambda*h <= hnext ) then ! <---------+              !  !   !  !
                 h = hnext                              !              !  !   !  !
                 exit inner                             !              !  !   !  !
              end if ! <--------------------------------+              !  !   !  !
           end if ! <--------------------------------------------------+  !   !  !
        end if ! <--------------------------------------------------------+   !  !
        !                                                                     !  !
        ! Newton either failed to converge or produced                        !  !
        ! an unsatisfactory solution => backtrack                             !  !
        h = FRACbacktrack * h                                                 !  !
        if ( h < EPSh ) then ! <---------------+                              !  !
           ! h << h0 (indefinite backtracking) !                              !  !
           stat = 5                            !                              !  !
           return                              !                              !  !
        end if ! <-----------------------------+                              !  !
        !                                                                     !  !
     end do inner ! <---------------------------------------------------------+  !
     !                                                                           !
     ! insert the new point                                                      !
     call insert_polyline_point( &                                               !
          uv, &                                                                  !
          xyz, &                                                                 !
          stat_insertion, &                                                      !
          polyline, &                                                            !
          i=polyline%np )                                                        !
     if ( stat_insertion > 0 ) then ! <--------+                                 !
        stat = stat_insertion                  !                                 !
        return                                 !                                 !
     end if ! <--------------------------------+                                 !
     !                                                                           !
     wprev = w                                                                   !
     !                                                                           !
  end do outer ! <---------------------------------------------------------------+

  ! last point
  call insert_polyline_point( &
       uv_endpoints(1:2,1:2,2), &
       xyz_endpoints(1:3,2), &
       stat_insertion, &
       polyline, &
       i=polyline%np ) 
  if ( stat_insertion > 0 ) then ! <--------+
     stat = stat_insertion                  !
     return                                 !
  end if ! <--------------------------------+

  IF ( DEBUG ) PRINT *,'TRACING OK, NP =',polyline%np

end subroutine trace_intersection_polyline
