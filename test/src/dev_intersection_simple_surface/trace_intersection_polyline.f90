subroutine trace_intersection_polyline( &
     surf, &
     uv_endpoints, &
     xyz_endpoints, &
     param_vector, &
     polyline, &
     hmin, &
     hmax )
  use mod_math
  use mod_diffgeom2
  implicit none
  LOGICAL, PARAMETER :: DEBUG = .false.
  !real(kind=fp), parameter                        :: tolchord = real( 5.e-4, kind=fp )
  !real(kind=fp), parameter                        :: FRACcurvature_radius = 2._fp * sqrt( tolchord*(2._fp - tolchord ) )
  !real(kind=fp), parameter                        :: tolh = real( 1e-2, kind=fp )
  !real(kind=fp), parameter                        :: tolhsqr = tolh**2
  real(kind=fp), parameter                        :: tolw = tolh
  real(kind=fp), parameter                        :: FRACbacktrack = 0.5_fp
  real(kind=fp), parameter                        :: EPSbacktrack = real( 1e-2, kind=fp )

  type(ptr_surface),                intent(in)    :: surf(2)
  real(kind=fp),                    intent(in)    :: param_vector(3)
  real(kind=fp),                    intent(in)    :: uv_endpoints(2,2,2) ! u/v, #surf, #point
  real(kind=fp),                    intent(in)    :: xyz_endpoints(3,2)
  type(type_intersection_polyline), intent(inout) :: polyline
  real(kind=fp), optional,          intent(in)    :: hmin, hmax
  real(kind=fp)                                   :: w0, Dw
  real(kind=fp)                                   :: curvature(2), uv_s(2,2,2), xyz_s(3,2)
  real(kind=fp)                                   :: h, EPSh, uv(2,2), xyz(3), w, wprev
  integer                                         :: stat
  integer                                         :: ipt
  
  !PRINT *,' PARAM_VECTOR =', PARAM_VECTOR
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
  END IF

  w0 = dot_product( param_vector, xyz_endpoints(:,1) )
  Dw = dot_product( param_vector, xyz_endpoints(:,2) ) - w0

  ! first point
  call insert_polyline_point( &
       uv_endpoints(:,:,1), &
       xyz_endpoints(:,1), &
       polyline )

  wprev = 0._fp
  outer_while : do 
     ! get tangent direction and curvature of the intersection point at the current point
     call diffgeom_intersection_curve( &
          surf, &
          polyline%uv(:,:,polyline%np), &
          uv_s, &
          xyz_s, &
          stat, &
          curvature )
     if ( stat > 0 ) then
        SELECT CASE (STAT)
        CASE (1)
           PRINT *,'trace_intersection_polyline : BRANCH POINT'
        CASE (2)
           PRINT *,'trace_intersection_polyline : ISOLATED CONTACT POINT'
        CASE (3)
           PRINT *,'trace_intersection_polyline : HIGH-ORDER CONTACT POINT'
        END SELECT
        STOP
     end if

     IF ( DEBUG) PRINT *,'DXYZ_DS.P =',dot_product( xyz_s(:,1), param_vector )
     
     ! target segment length
     h = FRACcurvature_radius / curvature(1)
     
     ! enforce bounds on segment length
     if ( present(hmin) ) then
        h = max( h, hmin ) 
     else
        h = max( h, tiny(h) )
     end if
     if ( present(hmax) ) then
        h = min( h, hmax )
     else
         h = min( h, huge(h) )
     end if
     IF ( DEBUG ) PRINT *,'HTARGET =',H

     ! check if the current is close enough to the end of the polyline
     if ( sum( (polyline%xyz(:,polyline%np) - xyz_endpoints(:,2))**2 ) < ((1.0 + tolh)*h)**2 ) then
        IF ( DEBUG ) THEN
           PRINT *,'CLOSE TO THE END, D =',NORM2( polyline%xyz(:,polyline%np) - xyz_endpoints(:,2) )
           PRINT *,' W* =', wprev, ', 1+TOLW = ', 1._fp - tolw
        END IF
        !if ( wprev > 1._fp - tolw ) exit outer_while
        exit outer_while
     end if

     ! compute next point using a Newton-Raphson algorithm
     EPSh = EPSbacktrack * h
     inner_while : do 
        uv = polyline%uv(:,:,polyline%np) + h * uv_s(:,1,:)
        
        call newton_intersection_polyline( &
             surf, &
             polyline%xyz(:,polyline%np), &
             h**2, &
             uv, &
             xyz, &
             stat )    
        
        if ( stat == 0 ) then
           ! Newton has converged, check whether the parameter w is monotonic along the polyline
           w = dot_product( param_vector, xyz )
           w = ( w - w0 ) / Dw
           if ( w > wprev .and. w <= 1._fp ) then
              exit inner_while
           else
              IF ( DEBUG ) THEN
                 PRINT *,'XYZ =',XYZ
                 PRINT *,'   W*=', W
                 PRINT *,'BACKTRACK...'
              END IF
           end if
        elseif ( stat == 2 ) then
           ! a singular Jacobian matrix has been encountered 
           STOP 'trace_intersection_polyline : singular Jacobian matrix in newton_intersection_polyline'
        end if
        
        h = FRACbacktrack * h
        if ( h < EPSh ) then
           CALL WRITE_POLYNOMIAL( surf(1)%ptr%x, 'trace_intersection_polyline/c1.cheb' )
           CALL WRITE_POLYNOMIAL( surf(2)%ptr%x, 'trace_intersection_polyline/c2.cheb' )

           OPEN( UNIT=13, FILE='trace_intersection_polyline/data.dat', ACTION='WRITE' )
           WRITE (13,*) 'UV_ENDPOINTS'
           DO IPT = 1,2
              WRITE (13,*) uv_endpoints(:,1,ipt)
              WRITE (13,*) uv_endpoints(:,2,ipt)
           END DO
           WRITE (13,*) 'XYZ_ENDPOINTS'
           DO IPT = 1,2
              WRITE (13,*) xyz_endpoints(:,ipt)
           END DO
           WRITE (13,*) 'PARAM_VECTOR'
           WRITE (13,*) param_vector
           CLOSE(13)

           OPEN( UNIT=13, FILE='trace_intersection_polyline/xyz_polyline.dat', ACTION='WRITE' )
           DO IPT = 1,POLYLINE%NP
              WRITE (13,*) POLYLINE%XYZ(:,IPT)
           END DO
           CLOSE(13)

           OPEN( UNIT=13, FILE='trace_intersection_polyline/uv_polyline.dat', ACTION='WRITE' )
           DO IPT = 1,POLYLINE%NP
              WRITE (13,*) POLYLINE%UV(:,:,IPT)
           END DO
           CLOSE(13)

           OPEN( UNIT=13, FILE='trace_intersection_polyline/tangent_polyline.dat', ACTION='WRITE' )
           DO IPT = 1,POLYLINE%NP
              call diffgeom_intersection_curve( &
                   surf, &
                   polyline%uv(:,:,ipt), &
                   uv_s, &
                   xyz_s, &
                   stat, &
                   curvature )
              WRITE (13,*) xyz_s(:,1), uv_s(:,1,:)
           END DO
           CLOSE(13)

           
           call diffgeom_intersection_curve( &
                surf, &
                uv_endpoints(:,:,2), &
                uv_s, &
                xyz_s, &
                stat, &
                curvature )
           PRINT *,'AT 2ND ENDPOINT, XYZ_S.P =',dot_product( xyz_s(:,1), param_vector )

           STOP 'trace_intersection_polyline : h << h0, indefinite backtracking'
        end if

     end do inner_while
     
     ! insert the new point
     call insert_polyline_point( &
          uv, &
          xyz, &
          polyline ) 
     IF ( DEBUG ) THEN
        PRINT *,'+1 POINT :'
        PRINT *,'     UV =',UV
        PRINT *,'    XYZ =',XYZ
        PRINT *,'     W* =',W
        PRINT *,''
     END IF

     wprev = w

  end do outer_while


  ! last point
  call insert_polyline_point( &
       uv_endpoints(:,:,2), &
       xyz_endpoints(:,2), &
       polyline )  


  IF ( .false. ) THEN !( DEBUG ) THEN
     CALL WRITE_POLYNOMIAL( surf(1)%ptr%x, 'trace_intersection_polyline/c1.cheb' )
     CALL WRITE_POLYNOMIAL( surf(2)%ptr%x, 'trace_intersection_polyline/c2.cheb' )

     OPEN( UNIT=13, FILE='trace_intersection_polyline/data.dat', ACTION='WRITE' )
     WRITE (13,*) 'UV_ENDPOINTS'
     DO IPT = 1,2
        WRITE (13,*) uv_endpoints(:,1,ipt)
        WRITE (13,*) uv_endpoints(:,2,ipt)
     END DO
     WRITE (13,*) 'XYZ_ENDPOINTS'
     DO IPT = 1,2
        WRITE (13,*) xyz_endpoints(:,ipt)
     END DO
     WRITE (13,*) 'PARAM_VECTOR'
     WRITE (13,*) param_vector
     CLOSE(13)

     OPEN( UNIT=13, FILE='trace_intersection_polyline/xyz_polyline.dat', ACTION='WRITE' )
     DO IPT = 1,POLYLINE%NP
        WRITE (13,*) POLYLINE%XYZ(:,IPT)
     END DO
     CLOSE(13)

     OPEN( UNIT=13, FILE='trace_intersection_polyline/uv_polyline.dat', ACTION='WRITE' )
     DO IPT = 1,POLYLINE%NP
        WRITE (13,*) POLYLINE%UV(:,:,IPT)
     END DO
     CLOSE(13)

     OPEN( UNIT=13, FILE='trace_intersection_polyline/tangent_polyline.dat', ACTION='WRITE' )
     DO IPT = 1,POLYLINE%NP
        call diffgeom_intersection_curve( &
          surf, &
          polyline%uv(:,:,ipt), &
          uv_s, &
          xyz_s, &
          stat, &
          curvature )
        WRITE (13,*) xyz_s(:,1), uv_s(:,1,:)
     END DO
     CLOSE(13)
     
     STOP
  END IF


  IF ( DEBUG ) THEN
     PRINT *,''; PRINT *,''; PRINT *,'';
  END IF

end subroutine trace_intersection_polyline
