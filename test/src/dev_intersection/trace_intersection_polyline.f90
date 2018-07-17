subroutine trace_intersection_polyline( &
     surf, &
     uvbox, &
     param_vector, &
     uv_endpoints, &
     xyz_endpoints, &
     stat, &
     polyline, &
     hmin, &
     hmax )
  use mod_math
  use mod_diffgeom2
  use mod_types_intersection
  use mod_tolerances
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  real(kind=fp), parameter                        :: FRACbacktrack = 0.5_fp
  real(kind=fp), parameter                        :: EPSbacktrack = real(1d-2, kind=fp)
  type(ptr_surface),                intent(in)    :: surf(2)
  real(kind=fp),                    intent(in)    :: uvbox(2,2,2)
  real(kind=fp),                    intent(in)    :: param_vector(3)
  real(kind=fp),                    intent(in)    :: uv_endpoints(2,2,2) ! u/v, #surf, #point
  real(kind=fp),                    intent(in)    :: xyz_endpoints(3,2)
  integer,                          intent(out)   :: stat
  type(type_intersection_polyline), intent(inout) :: polyline
  real(kind=fp), optional,          intent(in)    :: hmin, hmax
  real(kind=fp), dimension(4)                     :: lowerb, upperb
  real(kind=fp)                                   :: w0, Dw, w, wprev
  real(kind=fp)                                   :: duv_ds(2,2,2), dxyz_ds(3,2), curvature(2)
  real(kind=fp)                                   :: h_endpoints(2), h, EPSh
  real(kind=fp)                                   :: uv(2,2), xyz(3)
  integer                                         :: ipt

  stat = 0
  lowerb = reshape(uvbox(1,1:2,1:2), [4])
  upperb = reshape(uvbox(2,1:2,1:2), [4])
  IF ( DEBUG ) THEN
     PRINT *,''; PRINT *,'';
     PRINT *,'TRACE_INTERSECTION_POLYLINE'
     PRINT *,' PARAM_VECTOR =', PARAM_VECTOR
     PRINT *,' UV_ENDPOINTS ='
     PRINT *,uv_endpoints(:,:,1)
     PRINT *,uv_endpoints(:,:,2)
     PRINT *,'XYZ_ENDPOINTS ='
     PRINT *,xyz_endpoints(:,1)
     PRINT *,xyz_endpoints(:,2)
     PRINT *,'LOWERB =',LOWERB
     PRINT *,'UPPERB =',UPPERB
  END IF

  w0 = dot_product( param_vector, xyz_endpoints(:,1) )
  Dw = dot_product( param_vector, xyz_endpoints(:,2) ) - w0

  do ipt = 2,1,-1
     ! get tangent direction and curvature of the intersection point at endpoints
     ! (in reverse order, so data at first endpoint is kept in memory)
     call diffgeom_intersection_curve( &
          surf, &
          uv_endpoints(:,:,ipt), &
          duv_ds, &
          dxyz_ds, &
          stat, &
          curvature )
     h_endpoints(ipt) = FRACcurvature_radius / curvature(1)
  end do
  if ( present(hmin) ) h_endpoints = max(h_endpoints, hmin) 
  if ( present(hmax) ) h_endpoints = min(h_endpoints, hmax)

  
  ! first point
  call insert_polyline_point( &
       uv_endpoints(1:2,1:2,1), &
       xyz_endpoints(1:3,1), &
       stat, &
       polyline, &
       i=0 )
  if ( stat > 0 ) then
     stat = stat + 2 
     return
  end if

  wprev = 0._fp
  h = h_endpoints(1)
  outer : do
     ! check if the current is close enough to the end of the polyline
     if ( sum( (polyline%xyz(:,polyline%np) - xyz_endpoints(:,2))**2 ) < &
          min(h, h_endpoints(2))**2 ) exit outer

     ! compute next point
     EPSh = EPSbacktrack * h
     inner : do
        ! initial iterate
        uv = polyline%uv(:,:,polyline%np) + h*duv_ds(:,1,:)

        call newton_intersection_polyline( &
             surf, &
             lowerb, &
             upperb, &
             polyline%xyz(1:3,polyline%np), &
             h**2, &
             stat, &
             uv, &
             xyz )

        if ( stat == 0 ) then
           ! Newton has converged, check whether w is monotonic along the polyline
           w = dot_product( param_vector, xyz )
           w = ( w - w0 ) / Dw
           if ( is_in_open_interval(w, wprev, 1._fp) ) then
              ! get tangent direction and curvature of the intersection point at the current point
              call diffgeom_intersection_curve( &
                   surf, &
                   uv, &
                   duv_ds, &
                   dxyz_ds, &
                   stat, &
                   curvature )
              if ( h <= FRACcurvature_radius / curvature(1) ) then
                 h = FRACcurvature_radius / curvature(1)
                 exit inner
              end if
           end if
        elseif ( stat == 2 ) then
           ! a singular Jacobian matrix has been encountered 
           return
        end if

        ! Newton either failed to converge or produced an unsatisfactory solution => backtrack
        h = FRACbacktrack * h
        if ( h < EPSh ) then
           ! h << h0 (indefinite backtracking)
           stat = 1
           return
        end if

     end do inner

     ! insert the new point
     call insert_polyline_point( &
          uv, &
          xyz, &
          stat, &
          polyline, &
          i=polyline%np ) 
     if ( stat > 0 ) then
        stat = stat + 2 
        return
     end if

     wprev = w

  end do outer

  ! last point
  call insert_polyline_point( &
       uv_endpoints(1:2,1:2,2), &
       xyz_endpoints(1:3,2), &
       stat, &
       polyline, &
       i=polyline%np ) 
  if ( stat > 0 ) then
     stat = stat + 2 
     return
  end if

end subroutine trace_intersection_polyline
