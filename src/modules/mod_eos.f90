module mod_eos

  use mod_math

  implicit none

contains

  !-------------------------------------------------------------
  function normal_speed( &
       x, &
       y, &
       z, &
       m, &
       n )
    implicit none
    integer,                       intent(in) :: m, n
    real(kind=fp), dimension(m,n), intent(in) :: x, y, z
    real(kind=fp)                             :: normal_speed(m,n)

    !normal_speed(1:m,1:n) = 0.15d0
    !normal_speed = 0.15d0*(1.d0+ 0.3d0*(x+1.d0) - 0.08d0*(y+1.d0) + 0.12d0*(z+1.d0))
    !normal_speed = 0.15d0*(1.d0 + 0.3d0*cos(6.d0*(x+y+z)))
    !normal_speed = 0.15d0*(1.d0 + 0.15d0*cos(5.d0*(x+y+z)))
    !normal_speed = 0.15d0*(1.d0 + 0.05d0*cos(5.d0*(x+y+z)))
    normal_speed = 0.15d0*(1.d0 + 0.1d0*sin(6.d0*(x - 0.4d0 + y - 0.4d0 + z - 0.4d0)))
  end function normal_speed
  !-------------------------------------------------------------



  !-------------------------------------------------------------
  subroutine trace_border_polyline( &
       surf, &
       iborder, &
       np, &
       uv, &
       xyz )
    use mod_diffgeom
    implicit none
    type(type_surface), intent(in)  :: surf
    integer,            intent(in)  :: iborder
    integer,            intent(in)  :: np
    real(kind=fp),      intent(out) :: uv(2,np)
    real(kind=fp),      intent(out) :: xyz(3,np)
    integer                         :: ip

    select case (iborder)
    case (1)
       uv(1,1:np) = linspace(-1._fp, 1._fp, np)
       uv(2,1:np) = -1._fp
    case (2)
       uv(1,1:np) = 1._fp  
       uv(2,1:np) = linspace(-1._fp, 1._fp, np)
    case (3)
       uv(1,1:np) = linspace(1._fp, -1._fp, np)
       uv(2,1:np) = 1._fp
    case (4)
       uv(1,1:np) = -1._fp  
       uv(2,1:np) = linspace(1._fp, -1._fp, np)
    end select

    do ip = 1,np
       call eval( &
            xyz(1:3,ip), &
            surf, &
            uv(1:2,ip) )
    end do

end subroutine trace_border_polyline
  !-------------------------------------------------------------


  !-------------------------------------------------------------
  subroutine trace_border_polyline_adaptive( &
       surf, &
       iborder, &
       tolchord, &
       hmin, &
       hmax, &
       stat, &
       uv, &
       xyz, &
       np )
    use mod_util
    use mod_polynomial
    use mod_diffgeom
    use mod_tolerances
    implicit none
    LOGICAL, PARAMETER :: DEBUG = .false.
    real(kind=fp), parameter                  :: FRACbacktrack = 0.5_fp
    type(type_surface),         intent(in)    :: surf
    integer,                    intent(in)    :: iborder
    real(kind=fp),              intent(in)    :: tolchord
    real(kind=fp),              intent(in)    :: hmin, hmax
    integer,                    intent(out)   :: stat
    real(kind=fp), allocatable, intent(inout) :: uv(:,:)
    real(kind=fp), allocatable, intent(inout) :: xyz(:,:)
    integer,                    intent(out)   :: np
    real(kind=fp)                             :: FRACcurvature_radius 
    integer                                   :: ivar, ival, jvar
    real(kind=fp)                             :: w_ends(2), sign_dw
    type(type_curve)                          :: curv
    real(kind=fp)                             :: curvature
    real(kind=fp)                             :: h_ends(2), xyz_ends(3,2)
    real(kind=fp)                             :: w, h
    real(kind=fp)                             :: dist_from_end
    real(kind=fp)                             :: dxyz_dw(3)
    real(kind=fp)                             :: norm_dxyz_dw
    real(kind=fp)                             :: w_next, xyz_next(3), uv_next(2), h_next
    integer                                   :: ipt, ntmp

    IF ( DEBUG ) THEN
       PRINT *,'-->> TRACE_BORDER_POLYLINE_ADAPTIVE'
       PRINT *,'IBORDER =',IBORDER
       PRINT *,'TOLCHORD =',TOLCHORD
       PRINT *,'HMIN =',HMIN
       PRINT *,'HMAX =',HMAX
    END IF

    stat = 0
    FRACcurvature_radius = 2._fp*sqrt(tolchord*(2._fp - tolchord))
    
    select case (iborder)
    case (1)
       ivar = 2
       ival = 1
       w_ends(1:2) = [-1._fp, 1._fp]
    case (2)
       ivar = 1
       ival = 2
       w_ends(1:2) = [-1._fp, 1._fp]
    case (3)
       ivar = 2
       ival = 2
       w_ends(1:2) = [1._fp, -1._fp]
    case (4)
       ivar = 1
       ival = 1
       w_ends(1:2) = [1._fp, -1._fp]
    end select
    ! uv(ivar) = (-1)**ival
    ! uv(jvar) = w
    jvar = 1 + mod(ivar,2) ! free parameter

    sign_dw = sign(1._fp, w_ends(2) - w_ends(1))

    call bivar2univar( &
         surf%x, &
         curv%x, &
         ivar, &
         ival )
    call compute_deriv1(curv)
    call compute_deriv2(curv)

    ! endpoints
    do ipt = 1,2
       call eval_curvature_curve( &
            curv, &
            w_ends(ipt), &
            curvature )
       IF ( DEBUG ) PRINT *,'  ENDPOINT #',IPT,', CURVATURE =', CURVATURE
       h_ends(ipt) = FRACcurvature_radius/max(EPSfp, curvature)
       
       call eval( &
            xyz_ends(1:3,ipt), &
            curv, &
            w_ends(ipt) )
    end do
    h_ends = max(h_ends, hmin)
    h_ends = min(h_ends, hmax)
    IF ( DEBUG ) PRINT *,'H_ENDS =',H_ENDS

    ! first point
    np = 0
    call append_vec( &
         xyz_ends(1:3,1), &
         3, &
         xyz, &
         np )
    uv_next(ivar) = real((-1)**ival, kind=fp)
    uv_next(jvar) = w_ends(1)
    ntmp = np - 1
    call append_vec( &
         uv_next, &
         2, &
         uv, &
         ntmp )
    
    w = w_ends(1)
    h = h_ends(1)

    outer : do
       ! check if the current is close enough to the end of the polyline
       dist_from_end = sum((xyz(1:3,np) - xyz_ends(1:3,2))**2)
       !IF ( DEBUG ) PRINT *,' DIST_FROM_END =',SQRT(DIST_FROM_END)

       if ( dist_from_end < h**2 ) then ! <--------------------------------------+
          if ( dist_from_end < h_ends(2)**2 ) then ! <------------------------+  !
             ! current point too close from endpoint, simply remove it        !  !
             if ( dist_from_end < (TOLh * min(h, h_ends(2)))**2 ) np = np - 1 !  !
             !IF ( DEBUG ) PRINT *,'TOO CLOSE'
             exit outer                                                       !  !
          else ! -------------------------------------------------------------+  !
             h = sqrt(dist_from_end)                                          !  !
          end if ! <----------------------------------------------------------+  !
       end if ! <----------------------------------------------------------------+

       ! compute next point
       call evald1( &
            dxyz_dw, &
            curv, &
            w )
       norm_dxyz_dw = norm2(dxyz_dw)

       backtracking : do
          w_next = w + sign_dw*h/norm_dxyz_dw
          w_next = min(1._fp + EPSuv, max(-(1._fp + EPSuv), w_next))

          call eval( &
               xyz_next, &
               curv, &
               w_next )

          if ( sum((xyz_next - xyz(1:3,np))**2) > (h*(1.0 + TOLh))**2 ) then
             h = FRACbacktrack*h
             if ( h < 0.001*hmin ) then
                stat = 1
                return
             end if
             cycle
          end if

          call eval_curvature_curve( &
            curv, &
            w_next, &
            curvature )
          h_next = FRACcurvature_radius/max(EPSfp, curvature)
          h_next = min(hmax, max(hmin, h_next))

          if ( h_next < h ) then
             h = h_next
          else      
             exit backtracking
          end if
       end do backtracking

       ! add point
       IF ( DEBUG ) PRINT *,' + 1 POINT'
       call append_vec( &
         xyz_next, &
         3, &
         xyz, &
         np )
       ntmp = np - 1
       uv_next(ivar) = uv(ivar,1)
       uv_next(jvar) = w_next
       call append_vec( &
         uv_next, &
         2, &
         uv, &
         ntmp )
       w = w_next
       h = h_next
            
    end do outer

    ! last point
    call append_vec( &
         xyz_ends(1:3,2), &
         3, &
         xyz, &
         np )
    ntmp = np - 1
    uv_next(ivar) = uv(ivar,1)
    uv_next(jvar) = w_ends(2)
    call append_vec( &
         uv_next, &
         2, &
         uv, &
         ntmp )

    !! get uv coordinates
    !if ( allocated(uv) ) then
    !   if ( size(uv,1) < 2 .or. size(uv,2) < np ) deallocate(uv)
    !end if
    !if ( .not.allocated(uv) ) allocate(uv(2,np))  
    
  end subroutine trace_border_polyline_adaptive
  !-------------------------------------------------------------



  !-------------------------------------------------------------
  subroutine eos_from_surface( &
       timestep, &
       surf, &
       surf_eos )
    use mod_chebyshev
    use mod_diffgeom
    use mod_propagation
    implicit none
    real(kind=fp),      intent(in)  :: timestep
    type(type_surface), intent(in)  :: surf
    type(type_surface), intent(out) :: surf_eos
    integer                         :: m, n, k
    real(kind=fp), allocatable      :: xyz(:,:,:), speed(:,:), s2e(:,:,:)

    m = size(surf%x%coef,1)
    n = size(surf%x%coef,2)

    allocate(xyz(m,n,3), speed(m,n), s2e(m,n,3))

    call ifcht2( &
         surf%x%coef(1:m,1:n,1:3), &
         xyz, &
         m, &
         n, &
         3 )

    speed = normal_speed( &
         xyz(1:m,1:n,1), &
         xyz(1:m,1:n,2), &
         xyz(1:m,1:n,3), &
         m, &
         n )

    call eos_surface( &
         xyz, &
         speed, &
         m, &
         n, &
         timestep, &
         1, &!sens_normal
         s2e )

    do k = 1,3
       xyz(1:m,1:n,k) = xyz(1:m,1:n,k) + timestep*speed(1:m,1:n)*s2e(1:m,1:n,k)
    end do

    call reset_polynomial(surf_eos%x, 2, 1, [m-1,n-1], 3)
    call fcht2( &
         xyz, &
         surf_eos%x%coef(1:m,1:n,1:3), &
         m, &
         n, &
         3, &
         EPSmath )

    deallocate(xyz, speed, s2e)

  end subroutine eos_from_surface
  !-------------------------------------------------------------



  !-------------------------------------------------------------
  subroutine eos_from_polyline( &
       g, &
       dg, &
       er, &
       el, &
       m, &
       v, &
       n, &
       e )
    use mod_geometry
    USE MOD_UTIL
    implicit none
    LOGICAL, PARAMETER :: DEBUG = .false.
    integer,                       intent(in)  :: m
    real(kind=fp), dimension(3,m), intent(in)  :: g, dg, er, el
    integer,                       intent(in)  :: n
    real(kind=fp),                 intent(in)  :: v(n)
    real(kind=fp),                 intent(out) :: e(m,n,3)
    real(kind=fp), dimension(3,m)              :: o, eor, eol
    !real(kind=fp)                              :: ei(3,n)
    integer                                    :: i, j
    real(kind=fp)                              :: angle, aj, eor90d(3)
    INTEGER :: FID, K
    
    o = g + dg*spread(0.5_fp*sum((er + el - 2._fp*g)*dg,1)/sum(dg**2,1), dim=1, ncopies=3)
    
    IF ( DEBUG ) CALL WRITE_MATRIX( &
         TRANSPOSE(O(1:3,1:M)), &
         M, 3, &
         '../debug/eos_from_polyline_center.dat' )
    
    eor = er - o
    eol = el - o

    IF ( DEBUG ) THEN
       CALL GET_FREE_UNIT(FID)
       OPEN(UNIT=FID, &
            FILE='../debug/eos_from_polyline_localframe.dat', &
            ACTION='WRITE')
    END IF

    do i = 1,m      
       !call slerp( &
       !     eor(1:3,i), &
       !     eol(1:3,i), &
       !     v, &
       !     n, &
       !     ei )

       eor90d = norm2(eol(1:3,i))*cross(dg(1:3,i), eor(1:3,i)) / &
            (norm2(eor(1:3,i))*norm2(dg(1:3,i)))
       angle = atan2( &
            dot_product(eor90d, eol(1:3,i)),&
            dot_product(eor(1:3,i), eol(1:3,i)) &
            )
       if ( angle > 0._fp ) angle = angle - 2._fp*CSTpi

       IF ( DEBUG ) WRITE (FID,*) EOR(1:3,I), EOR90D, ANGLE

       do j = 1,n
          !e(i,j,1:3) = o(1:3,i) + ei(1:3,j)
          aj = v(j)*angle
          e(i,j,1:3) = o(1:3,i) + &
               cos(aj)*eor(1:3,i) + sin(aj)*eor90d
       end do
    end do

    IF ( DEBUG ) CLOSE(FID)

    IF ( DEBUG ) THEN
       CALL GET_FREE_UNIT(FID)
       OPEN(UNIT=FID, &
            FILE='../debug/eos_from_polyline_cgl_tpgrid.dat', &
            ACTION='WRITE')
       WRITE (FID, *), M, N, 3
       DO K = 1,3
          DO J = 1,N
             DO I = 1,M
                WRITE (FID,*) E(I,J,K)
             END DO
          END DO
       END DO
       CLOSE(FID)
    END IF

  end subroutine eos_from_polyline
  !-------------------------------------------------------------





  !-------------------------------------------------------------
  subroutine eos_from_curve( &
       mode, &
       curve, &
       head, &
       tail, &
       enve_rl, &
       x, &
       xyzc, &
       m, &
       n )
    use mod_chebyshev
    use mod_diffgeom
    use mod_types_intersection
    implicit none
    integer,                               intent(in)    :: mode
    type(type_intersection_curve), target, intent(in)    :: curve
    integer,                               intent(in)    :: head
    integer,                               intent(in)    :: tail
    type(ptr_surface),                     intent(in)    :: enve_rl(2)
    real(kind=fp), allocatable,            intent(inout) :: x(:)
    real(kind=fp), allocatable,            intent(inout) :: xyzc(:,:,:)
    integer,                               intent(out)   :: m, n
    type(type_intersection_polyline), pointer            :: polyline => null()
    integer                                              :: step, np, degr
    real(kind=fp), allocatable                           :: y(:,:)
    real(kind=fp), allocatable                           :: c(:,:), d(:,:)
    real(kind=fp)                                        :: cond, errL2, errLinf
    real(kind=fp), allocatable                           :: u(:)
    real(kind=fp), allocatable                           :: g(:,:), dg(:,:)
    real(kind=fp), allocatable                           :: uv(:,:), erl(:,:,:)
    integer                                              :: irl, ipt

    polyline => curve%polyline
    
    !! fit Chebyshev polynomial to intersection curve
    step = sign(1,tail-head)
    np = 1 + (tail-head)/step
    degr = max(1, int(0.25*real(np)) - 1)
    m = degr+1
    allocate(c(m,7))

    ! curve parameter = chordal arclength approximation (normalized to [-1,1])
    !x(1) = 0._fp
    !x(2:np) = sqrt(sum((&
    !     polyline%xyz(1:3,head+step:tail:step) - &
    !     polyline%xyz(1:3,head:tail-step:step))**2, 1))
    !do ipt = 2,np
    !   x(ipt) = x(ipt-1) + x(ipt)
    !end do
    !x = 2._fp*x/x(np) - 1._fp
    call intersection_segment_parameter( &
         curve, &
         head, &
         tail, &
         mode, &
         x )

    allocate(y(np,7))
    y(1:np,1:2) = transpose(polyline%uv(1:2,1,head:tail:step))
    y(1:np,3:4) = transpose(polyline%uv(1:2,2,head:tail:step))
    y(1:np,5:7) = transpose(polyline%xyz(1:3,head:tail:step))
    call chebfit1( &
         x, &
         y, &
         np, &
         7, &
         c, &
         degr, &
         cond, &
         errL2, &
         errLinf )
    deallocate(y)
    PRINT *,'COND =', COND, 'ERR_(L2,LINF) =', ERRL2, ERRLINF

    !! compute the polynomial's derivative
    allocate(d(degr,3))
    call chebdiff( &
         c(1:m,5:7), &
         d, &
         degr, &
         3 )

    !! evaluate intersection curve at cgl points
    allocate(g(m,3), dg(m,3))
    call ifcht1( &
         c(1:m,5:7), &
         g, &
         m, &
         3 )
    call ifcht1( &
         d, &
         dg, &
         m, &
         3 )

    !! evaluate right and left boundary curves
    allocate(uv(m,2), erl(m,3,2))
    do irl = 1,2
       call ifcht1( &
            c(1:m,2*irl-1:2*irl), &
            uv, &
            m, &
            2 )

       do ipt = 1,m
          call eval( &
               erl(ipt,1:3,irl), &
               enve_rl(irl)%ptr, &
               uv(ipt,1:2) )
       end do
    end do

    !! contruct portion of canal surfaces between right/left boundary curves
    n = m
    allocate(u(n), xyzc(m,n,3))
    call cgl_nodes( &
         u, &
         degr, &
         1._fp, &
         0._fp) !-1._fp )

    call eos_from_polyline( &
         transpose(g), &
         transpose(dg), &
         transpose(erl(1:m,1:3,1)), &
         transpose(erl(1:m,1:3,2)), &
         m, &
         u, &
         n, &
         xyzc )

  end subroutine eos_from_curve
  !-------------------------------------------------------------





  !-------------------------------------------------------------
  subroutine long_lat_patch_from_points( &
       center, &
       xyz, &
       n, &
       degr, &
       stat, &
       surf, &
       uv, &
       longlat_range, &
       solidangle )
    use mod_diffgeom
    use mod_geometry
    use mod_chebyshev
    use mod_polynomial
    USE MOD_UTIL
    implicit none
    LOGICAL :: DEBUG = .false.
    real(kind=fp), parameter          :: mrg = 0.05_fp
    real(kind=fp),      intent(in)    :: center(3)
    integer,            intent(in)    :: n
    real(kind=fp),      intent(inout) :: xyz(3,n)
    integer,            intent(in)    :: degr
    integer,            intent(out)   :: stat
    type(type_surface), intent(inout) :: surf
    real(kind=fp),      intent(out)   :: uv(2,n)
    real(kind=fp),      intent(out)   :: longlat_range(2)
    real(kind=fp),      intent(out)   :: solidangle
    real(kind=fp)                     :: xyzavg(3), xyzavg_sqr, ravg, r
    real(kind=fp)                     :: axes(3,3)
    real(kind=fp)                     :: xyzproj(3,n), dot
    real(kind=fp)                     :: boxcenter(2), boxranges(2), boxaxes(3,3)
    real(kind=fp), dimension(2)       :: minuv, maxuv
    real(kind=fp), dimension(degr+1)  :: tcgl, cosucgl, sinucgl, cosvcgl, sinvcgl
    real(kind=fp)                     :: xyzcgl(3,(degr+1)**2)
    integer                           :: i, j, k, m
    real(kind=fp)                     :: longlat_center(2)
    INTEGER :: FID

    IF ( DEBUG ) CALL GET_FREE_UNIT(FID)  
    m = degr+1

    ! translate center to the origin and project onto the unit sphere
    ! and compute centroid of xyz points
    xyzavg(1:3) = 0._fp
    ravg = 0._fp
    do i = 1,n
       xyz(1:3,i) = xyz(1:3,i) - center
       r = norm2(xyz(1:3,i))
       !IF ( DEBUG ) PRINT *,'R =',R
       xyz(1:3,i) = xyz(1:3,i)/r
       xyzavg(1:3) = xyzavg(1:3) + xyz(1:3,i)
       ravg = ravg + r
    end do
    xyzavg_sqr = sum(xyzavg**2)
    if ( xyzavg_sqr < EPSmath**2 ) then
       stat = 2
       return
    end if
    
    xyzavg = xyzavg/sqrt(xyzavg_sqr)
    ravg = ravg/real(n, kind=fp)

    IF ( DEBUG ) CALL WRITE_MATRIX(TRANSPOSE(XYZ), N, 3, &
         '/d/bandrieu/GitHub/FFTsurf/test/demo_EoS_brep/debug/debuglonglat_xyz.dat')

    ! determine the sphere's orientation
    axes(1:3,1) = xyzavg
    call complete_orthonormal_basis(axes(1:3,1), axes(1:3,2), axes(1:3,3))
    axes(1:3,1:3) = axes(1:3,[2,3,1])
    IF ( DEBUG ) CALL WRITE_MATRIX(AXES, 3, 3, &
         '/d/bandrieu/GitHub/FFTsurf/test/demo_EoS_brep/debug/debuglonglat_axes1.dat')

    ! central projection onto the plane tangent to the unit sphere at xyzavg
    ! in the meantime, check if all points have positive dot product with xyzavg
    ! (i.e. all the points lie in the hemisphere with pole xyzavg)
    do i = 1,n
       dot = dot_product(xyz(1:3,i), xyzavg)
       if ( .FALSE. ) THEN!dot_product(xyz(1:3,i), xyzavg) < EPSfp ) then
          stat = 1 
          return
       else
          stat = 0
          xyzproj(1:3,i) = xyz(1:3,i)/dot
       end if
    end do
    IF ( DEBUG ) CALL WRITE_MATRIX(TRANSPOSE(XYZPROJ), N, 3, &
         '/d/bandrieu/GitHub/FFTsurf/test/demo_EoS_brep/debug/debuglonglat_xyzproj.dat')

    ! compute coordinates in the plane's basis
    xyzproj(1:2,1:n) = matmul(transpose(axes(1:3,1:2)), xyzproj(1:3,1:n))
    IF ( DEBUG ) CALL WRITE_MATRIX(TRANSPOSE(XYZPROJ(1:2,1:n)), N, 2, &
         '/d/bandrieu/GitHub/FFTsurf/test/demo_EoS_brep/debug/debuglonglat_xyplane.dat')

    ! oriented bounding box with minimum area
    call minimum_area_OBB_2d( &
         xyzproj(1:2,1:n), &
         n, &
         boxcenter, &
         boxranges, &
         boxaxes(1:2,1:2) )
    boxaxes(1:2,3) = 0._fp
    boxaxes(3,1:2) = 0._fp
    boxaxes(3,3) = 1._fp
    IF ( DEBUG ) THEN
       OPEN(UNIT=FID, &
            FILE='/d/bandrieu/GitHub/FFTsurf/test/demo_EoS_brep/debug/debuglonglat_obb.dat', &
            ACTION='WRITE')
       WRITE(FID,*) BOXCENTER
       WRITE(FID,*) BOXRANGES
       WRITE(FID,*) BOXAXES(1,1:2)
       WRITE(FID,*) BOXAXES(2,1:2)
       CLOSE(FID)
    END IF

    ! change sphere orientation to minimize the longitude-latitude rectangle's area
    axes = matmul(axes, boxaxes)
    axes = axes(1:3,[3,1,2])
    IF ( DEBUG ) CALL WRITE_MATRIX(AXES, 3, 3, &
         '/d/bandrieu/GitHub/FFTsurf/test/demo_EoS_brep/debug/debuglonglat_axes2.dat')

    ! compute longitude-latitude coordinates
    xyzproj(1:3,1:n) = matmul(transpose(axes), xyz(1:3,1:n))
    uv(1,1:n) = atan2(xyzproj(2,1:n), xyzproj(1,1:n))
    uv(2,1:n) = asin(max(-1._fp, min(1._fp, xyzproj(3,1:n))))
    IF ( DEBUG ) CALL WRITE_MATRIX(TRANSPOSE(UV), N, 2, &
         '/d/bandrieu/GitHub/FFTsurf/test/demo_EoS_brep/debug/debuglonglat_tl.dat')

    ! longitude-latitude bounding box
    minuv = minval(uv, dim=2)
    maxuv = maxval(uv, dim=2)
    if ( .false. ) then
      longlat_center(1:2) = 0._fp
      longlat_range = (1._fp + mrg)*max(abs(maxuv), abs(minuv))
    else
      longlat_center = 0.5_fp*(minuv + maxuv)
      longlat_range = (1._fp + mrg)*0.5_fp*(maxuv - minuv)
    end if
    do i = 1,n
      uv(1:2,i) = (uv(1:2,i) - longlat_center)/longlat_range
    end do
    IF ( DEBUG ) CALL WRITE_MATRIX(TRANSPOSE(UV), N, 2, &
         '/d/bandrieu/GitHub/FFTsurf/test/demo_EoS_brep/debug/debuglonglat_uv.dat')
    solidangle = longlat_range(1) * sin(longlat_range(2)) / CSTpi

    ! Chebyshev patch    
    call cgl_nodes(tcgl, degr)
    cosucgl = cos(longlat_center(1) + longlat_range(1)*tcgl)
    sinucgl = sin(longlat_center(1) + longlat_range(1)*tcgl)
    cosvcgl = cos(longlat_center(2) + longlat_range(2)*tcgl)
    sinvcgl = sin(longlat_center(2) + longlat_range(2)*tcgl)

    do j = 1,m
       do i = 1,m
          k = m*(j-1) + i
          xyzcgl(1,k) = cosucgl(i)*cosvcgl(j)
          xyzcgl(2,k) = sinucgl(i)*cosvcgl(j)
          xyzcgl(3,k) = sinvcgl(j)
       end do
    end do

    xyzcgl = matmul(axes, xyzcgl)

    do i = 1,m*m
       xyzcgl(1:3,i) = center(1:3) + ravg*xyzcgl(1:3,i) 
    end do

    call reset_polynomial(surf%x, 2, 1, [degr, degr], 3)
    call fcht2( &
         reshape(transpose(xyzcgl), [m,m,3]), &
         surf%x%coef(1:m,1:m,1:3), &
         m, &
         m, &
         3, &
         EPSmath )
    IF ( DEBUG ) CALL WRITE_POLYNOMIAL(surf%x, &
    '/d/bandrieu/GitHub/FFTsurf/test/demo_EoS_brep/debug/debuglonglat_surf.cheb')

    do i = 1,n
       xyz(1:3,i) = center(1:3) + ravg*xyz(1:3,i)
    end do

    IF ( DEBUG ) PAUSE

  end subroutine long_lat_patch_from_points
  !-------------------------------------------------------------










  !-------------------------------------------------------------
  subroutine eos_from_polyline2( &
       g, &
       dg, &
       er, &
       el, &
       m, &
       v, &
       n, &
       e )
    use mod_geometry
    USE MOD_UTIL
    implicit none
    integer,                       intent(in)  :: m
    real(kind=fp), dimension(3,m), intent(in)  :: g, dg, er, el
    integer,                       intent(in)  :: n
    real(kind=fp),                 intent(in)  :: v(n)
    real(kind=fp),                 intent(out) :: e(m,n,3)
    real(kind=fp), dimension(m)                :: r, ru
    !real(kind=fp), dimension(3,m)              :: o, psi1, psi2
    !real(kind=fp)                              :: ei(3,n)
    real(kind=fp)                              :: err(2,m,n)
    integer                                    :: i, j

    r = 0.5_fp*(sqrt(sum((er - g)**2,1)) + sqrt(sum((el - g)**2,1)))
    ru = -0.5_fp*sum((er + el - 2._fp*g)*dg,1) / r
    !o = g + dg*spread(-r*ru/sum(dg**2,1), dim=1, ncopies=3)

    !o = g + dg*spread(0.5_fp*sum((er + el - 2._fp*g)*dg,1)/sum(dg**2,1), dim=1, ncopies=3)
    !eor = er - o
    !eol = el - o
    !psi1 = er - o
    !psi2 = el - o

    !do i = 1,m
    !   cosa = max(-1._fp, min(1._fp, &
    !        dot_product(psi1(1:3,i), psi2(1:3,i)) / &
    !        ( norm2(psi1(1:3,i))*norm2(psi2(1:3,i)) ) &
    !        ))
    !   sina = sqrt(1._fp - cosa**2)

    !call slerp( &
    !     eor(1:3,i), &
    !     eol(1:3,i), &
    !     v, &
    !     n, &
    !     ei )
    !do j = 1,n
    !   e(i,j,1:3) = o(1:3,i) + ei(1:3,j)
    !end do
    !end do


    !ERR(1:2,1:M,1:N) = 0._FP
    !DO J = 1,N
    !   DO I = 1,M
    !      ERR(1,I,J) = ABS( NORM2(E(I,J,1:3) - G(1:3,I)) - R(I) )
    !      ERR(2,I,J) = ABS( DOT_PRODUCT(E(I,J,1:3) - G(1:3,I), DG(1:3,I)) + RU(I)*R(I) )
    !   END DO
    !END DO

    !CALL WRITE_MATRIX(ERR(1,1:M,1:N), M, N, 'demo_EoS_brep/debug/err_1_rad.dat')
    !CALL WRITE_MATRIX(ERR(2,1:M,1:N), M, N, 'demo_EoS_brep/debug/err_1_dot.dat')

    CALL eos_from_polyline( &
         g, &
         dg, &
         er, &
         el, &
         m, &
         v, &
         n, &
         e )

    ERR(1:2,1:M,1:N) = 0._FP
    DO J = 1,N
       DO I = 1,M
          ERR(1,I,J) = ABS( NORM2(E(I,J,1:3) - G(1:3,I)) - R(I) )
          ERR(2,I,J) = ABS( DOT_PRODUCT(E(I,J,1:3) - G(1:3,I), DG(1:3,I)) + RU(I)*R(I) )
       END DO
    END DO

    CALL WRITE_MATRIX(E(1:M,1:N,1), M, N, 'demo_EoS_brep/debug/eos_from_polyline_x.dat')
    CALL WRITE_MATRIX(E(1:M,1:N,2), M, N, 'demo_EoS_brep/debug/eos_from_polyline_y.dat')
    CALL WRITE_MATRIX(E(1:M,1:N,3), M, N, 'demo_EoS_brep/debug/eos_from_polyline_z.dat')
    CALL WRITE_MATRIX(ERR(1,1:M,1:N), M, N, 'demo_EoS_brep/debug/eos_from_polyline_err_rad.dat')
    CALL WRITE_MATRIX(ERR(2,1:M,1:N), M, N, 'demo_EoS_brep/debug/eos_from_polyline_err_dot.dat')

    PAUSE

  end subroutine eos_from_polyline2
  !-------------------------------------------------------------









  !-------------------------------------------------------------
  subroutine make_partial_eos( &
       parameterization_mode, &
       surf, &
       nsurf, &
       interdata, &
       curvetype, &
       brep, &
       timestep, &
       surf_new, &
       nsurf_new, &
       interdata_new )
    use mod_diffgeom
    use mod_types_intersection
    use mod_types_brep
    use mod_halfedge
    use mod_tolerances
    use mod_propagation
    use mod_brep
    use mod_intersection
    use mod_util
    implicit none
    integer,                                 intent(in)    :: parameterization_mode
    integer,                                 intent(in)    :: nsurf
    type(type_surface), target,              intent(in)    :: surf(nsurf)
    type(type_intersection_data), target,    intent(in)    :: interdata
    integer,                                 intent(in)    :: curvetype(interdata%nc)
    type(type_brep),                         intent(in)    :: brep
    real(kind=fp),                           intent(in)    :: timestep
    type(type_surface), allocatable, target, intent(inout) :: surf_new(:)
    integer,                                 intent(out)   :: nsurf_new
    type(type_intersection_data), target,    intent(inout) :: interdata_new

    integer                                                :: isurf, jsurf
    type(ptr_surface), allocatable                         :: surf_eos(:), edge_eos(:), vert_eos(:)
    integer                                                :: m, n
    real(kind=fp), allocatable                             :: xyz(:,:,:), speed(:,:), direction_eos(:,:,:)

    integer                                                :: iedge
    integer                                                :: edgetype(brep%ne)
    type(type_intersection_curve), pointer                 :: curve => null(), curve_new => null()
    type(ptr_surface)                                      :: enve_rl(2), surfpair(2)
    integer                                                :: head, tail, sens, np
    real(kind=fp)                                          :: uvendpoint(2,2), xyzendpoint(3)
    integer                                                :: idendpoint(2)
    integer                                                :: ipoint, jpoint, icurv, jcurv
    real(kind=fp), allocatable                             :: abscissa(:)

    integer                                                :: ivert
    integer                                                :: ihedg(2), iface, jface
    integer                                                :: stat, nptot, nconvex, iborder
    integer, allocatable                                   :: iconvex(:,:) ![iedge;iborder]
    real(kind=fp), dimension(:,:), allocatable             :: xyzlonglat, uveosedg, tmp, uvlonglat
    real(kind=fp)                                          :: longlat_range(2), solidangle
    integer                                                :: i, j
    CHARACTER(3) :: STRNUM3

    ! >>> ---------- EoS from surfaces
    allocate(surf_new(nsurf+brep%ne+brep%nv))
    allocate(surf_eos(nsurf))
    nsurf_new = 0
    do isurf = 1,nsurf
       nsurf_new = nsurf_new + 1
       m = 2*size(surf(isurf)%x%coef,1)
       n = 2*size(surf(isurf)%x%coef,2)

       allocate(xyz(m,n,3), speed(m,n), direction_eos(m,n,3))
       call ifcht2( &
            surf(isurf)%x%coef, &
            xyz, &
            m, &
            n, &
            3 )

       speed = normal_speed( &
            xyz(1:m,1:n,1), &
            xyz(1:m,1:n,2), &
            xyz(1:m,1:n,3), &
            m, &
            n )

       call eos_surface2( &
            xyz, &
            speed, &
            m, &
            n, &
            timestep, &
            1, &
            direction_eos )

       do i = 1,3
          xyz(1:m,1:n,i) = xyz(1:m,1:n,i) + timestep*speed(1:m,1:n)*direction_eos(1:m,1:n,i)
       end do

       call reset_polynomial(surf_new(nsurf_new)%x, 2, 1, [m-1,n-1], 3)
       call fcht2( &
            xyz, &
            surf_new(nsurf_new)%x%coef(1:m,1:n,1:3), &
            m, &
            n, &
            3, &
            EPSmath )
       deallocate(xyz, speed, direction_eos)

       surf_eos(isurf)%ptr => surf_new(nsurf_new)
       call economize2( &
            surf_eos(isurf)%ptr%x, &
            EPSmath )
       call compute_deriv1(surf_eos(isurf)%ptr)
       call compute_deriv2(surf_eos(isurf)%ptr)

       WRITE (STRNUM3, '(I3.3)') ISURF
       CALL WRITE_POLYNOMIAL( &
            SURF_EOS(ISURF)%PTR%X, &
            '../debug/eos_surf_c_'//strnum3 //'.cheb' )
    end do
    ! ----------<<<


    ! >>> ---------- Mark edges as smooth, concave or convex
    do iedge = 1,brep%ne
       do icurv = 1,interdata%nc
          if ( associated(brep%edges(iedge)%curve, interdata%curves(icurv)) ) then
             edgetype(iedge) = curvetype(icurv)
             exit
          end if
       end do
    end do
    ! ----------<<<


    ! >>> ---------- EoS from convex edges
    allocate(edge_eos(brep%ne))
    do iedge = 1,brep%ne
       PRINT *,'IEDGE =',IEDGE,', TYPE =',EDGETYPE(IEDGE)
       if ( edgetype(iedge) == 1 ) cycle ! concave edge

       curve => brep%edges(iedge)%curve
       do jsurf = 1,2
          do isurf = 1,nsurf
             if ( associated(curve%surf(jsurf)%ptr, surf(isurf)) ) then
                enve_rl(jsurf)%ptr => surf_eos(isurf)%ptr
                exit
             end if
          end do
       end do

       if ( edgetype(iedge) == 0 ) then ! smooth edge
          ! preserve smooth intersection curve
          np = curve%polyline%np
          do jpoint = 1,2
             ipoint = 1 + (jpoint-1)*(np-1)
             uvendpoint(1:2,1:2) = curve%polyline%uv(1:2,1:2,ipoint)
             call eval( &
                  xyzendpoint(1:3), &
                  enve_rl(1)%ptr, &
                  uvendpoint(1:2,1) )
             call add_intersection_point( &
                  uvendpoint, &
                  xyzendpoint, &
                  enve_rl, &
                  2, &
                  interdata_new, &
                  idendpoint(jpoint) )
          end do

          call add_intersection_curve( &
               interdata_new, &
               [0._fp, 0._fp, 0._fp], &
               idendpoint, &
               spread(spread([-1._fp-EPSuv, 1._fp+EPSuv], 2, 2), 3, 2) )
          curve_new => interdata_new%curves(interdata_new%nc)
          curve_new%smooth = .true.
          curve_new%surf = enve_rl
          curve_new%isplit(2,1:2) = [1,np]
          allocate(curve_new%iedge(1))
          curve_new%iedge(:) = 0
          allocate(curve_new%polyline)
          curve_new%polyline%np = np
          allocate(curve_new%polyline%uv(2,2,np))
          curve_new%polyline%uv(1:2,1:2,1:np) = curve%polyline%uv(1:2,1:2,1:np)
          ! update xyz
          allocate(curve_new%polyline%xyz(3,np))
          do ipoint = 1,np
             call eval( &
                  curve_new%polyline%xyz(1:3,ipoint), &
                  curve_new%surf(1)%ptr, &
                  curve_new%polyline%uv(1:2,1,ipoint) )
          end do

       else ! convex edge
          nsurf_new = nsurf_new + 1

          call get_polyline_endpoints( &
               brep, &
               [iedge,1], &
               head, &
               tail, &
               sens, &
               np )

          !allocate(abscissa(np))
          if ( parameterization_mode < 0 ) then
             call eos_from_curve2( &
                  32, &
                  curve, &
                  head, &
                  tail, &
                  enve_rl, &
                  abscissa, &
                  xyz, &
                  m, &
                  n )
          else
             call eos_from_curve( &
                  parameterization_mode, &
                  curve, &
                  head, &
                  tail, &
                  enve_rl, &
                  abscissa, &
                  xyz, &
                  m, &
                  n )
          end if

          !PRINT *,'EOS_FROM_CURVE OK'
          
          call reset_polynomial(surf_new(nsurf_new)%x, 2, 1, [m-1,n-1], 3)
          call fcht2( &
               xyz, &
               surf_new(nsurf_new)%x%coef(1:m,1:n,1:3), &
               m, &
               n, &
               3, &
               EPSmath )

          edge_eos(iedge)%ptr => surf_new(nsurf_new)
          call economize2( &
               edge_eos(iedge)%ptr%x, &
               EPSmath )
          call compute_deriv1(edge_eos(iedge)%ptr)
          call compute_deriv2(edge_eos(iedge)%ptr)

          WRITE (STRNUM3, '(I3.3)') IEDGE
          CALL WRITE_POLYNOMIAL( &
               EDGE_EOS(IEDGE)%PTR%X, &
               '../debug/eos_edge_c_'//strnum3 //'.cheb' )

          ! 1st curve (R-E), 2nd curve (E-L)
          do jcurv = 1,2
             surfpair(jcurv)%ptr => enve_rl(jcurv)%ptr
             surfpair(1+mod(jcurv,2))%ptr => edge_eos(iedge)%ptr
             do jpoint = 1,2
                ipoint = 1 + (jpoint-1)*(np-1)
                uvendpoint(1:2,jcurv) = curve%polyline%uv(1:2,jcurv,ipoint)
                uvendpoint(1,1+mod(jcurv,2)) = real((-1)**jpoint, kind=fp)
                uvendpoint(2,1+mod(jcurv,2)) = real((-1)**jcurv, kind=fp)
                call eval( &
                     xyzendpoint(1:3), &
                     surfpair(1)%ptr, &
                     uvendpoint(1:2,1) )
                call add_intersection_point( &
                     uvendpoint, &
                     xyzendpoint, &
                     surfpair, &
                     2, &
                     interdata_new, &
                     idendpoint(jpoint) )
             end do

             call add_intersection_curve( &
                  interdata_new, &
                  [0._fp, 0._fp, 0._fp], &
                  idendpoint, &
                  spread(spread([-1._fp-EPSuv, 1._fp+EPSuv], 2, 2), 3, 2) )
             curve_new => interdata_new%curves(interdata_new%nc)
             curve_new%smooth = .true.
             curve_new%surf = surfpair
             curve_new%isplit(2,1:2) = [1,np]
             allocate(curve_new%iedge(1))
             curve_new%iedge(:) = 0
             allocate(curve_new%polyline)
             curve_new%polyline%np = np
             allocate(curve_new%polyline%uv(2,2,np))
             curve_new%polyline%uv(1:2,jcurv,1:np) = curve%polyline%uv(1:2,jcurv,1:np)
             curve_new%polyline%uv(1,1+mod(jcurv,2),1:np) = abscissa(1:np)
             curve_new%polyline%uv(2,1+mod(jcurv,2),1:np) = real((-1)**jcurv, kind=fp)
             ! update xyz
             allocate(curve_new%polyline%xyz(3,np))
             do ipoint = 1,np
                call eval( &
                     curve_new%polyline%xyz(1:3,ipoint), &
                     curve_new%surf(1)%ptr, &
                     curve_new%polyline%uv(1:2,1,ipoint) )
             end do
          end do

          deallocate(abscissa, xyz)
       end if
    end do
    ! ----------<<<



    ! >>> ---------- EoS from convex vertices
    allocate(vert_eos(brep%nv))
    np = 30 ! HARD-CODED ---> in the future, use adaptive sampling
    allocate(xyzlonglat(3,np), uveosedg(2,np))
    do ivert = 1,brep%nv
       PRINT *, 'IVERT =',IVERT
       ihedg = brep%verts(ivert)%halfedge
       iface = get_face(brep, ihedg)
       nconvex = 0
       nptot = 0
       ! cycle incident halfedges:
       !    for each convex edge, store:
       !       - xyz polyline vertex coordinates
       !       - uv polyline vertex coordinates in incident eos_surf
       do
          PRINT *,ihedg(1), edgetype(ihedg(1)), BREP%EDGES(IHEDG(1))%HYPEREDGE
          if ( edgetype(ihedg(1)) == 2 ) then ! current incident edge is convex
             if ( ihedg(2) == 1 ) then
                iborder = 4
             else
                iborder = 2
             end if

             if ( .not.allocated(xyzlonglat) ) allocate(xyzlonglat(3,nptot+np))
             if ( .not.allocated(uveosedg) ) allocate(uveosedg(2,nptot+np))

             if ( size(xyzlonglat,2) < nptot+np ) then
                call move_alloc(from=xyzlonglat, to=tmp)
                allocate(xyzlonglat(3,nptot+np))
                xyzlonglat(1:3,1:nptot) = tmp(1:3,1:nptot)
                deallocate(tmp)
             end if
             if ( size(uveosedg,2) < nptot+np ) then
                call move_alloc(from=uveosedg, to=tmp)
                allocate(uveosedg(2,nptot+np))
                uveosedg(1:2,1:nptot) = tmp(1:2,1:nptot)
                deallocate(tmp)
             end if

             allocate(tmp(5,np))
             call trace_border_polyline( &
                  edge_eos(ihedg(1))%ptr, &
                  iborder, &
                  np, &
                  tmp(1:2,1:np), &
                  tmp(3:5,1:np) )
             xyzlonglat(1:3,nptot+1:nptot+np) = tmp(3:5,1:np)
             uveosedg(1:2,nptot+1:nptot+np) = tmp(1:2,1:np)

             call append_vec_i( &
                  [ihedg(1),iborder,nptot+1,nptot+np], &
                  4, &
                  iconvex, &
                  nconvex )

             nptot = nptot + np
             deallocate(tmp)
          end if

          ! traverse halfedges counter-clockwise
          ihedg = get_prev(brep, ihedg) ! previous halfedge
          ihedg = get_twin(ihedg) ! twin halfedge
          jface = get_face(brep, ihedg)

          if ( jface == iface ) exit
       end do

       if ( nconvex >= 2 ) then
          !PRINT *, IVERT, brep%edges(iconvex(1,1))%hyperedge, brep%edges(iconvex(2,1))%hyperedge
          if ( nconvex == 2 .and. &
               brep%edges(iconvex(1,1))%hyperedge == brep%edges(iconvex(1,2))%hyperedge ) then
             ! both incident convex edges are on the same hyperedge,
             ! we just add a tangential instersection curve between
             ! the new surfaces associated to those edges
             do iedge = 1,2
                surfpair(iedge)%ptr => edge_eos(iconvex(1,1+mod(iedge,2)))%ptr
             end do
             ! for now we assume np is equal for both incident curves' polylines
             ! (=> in particular ntot = 2*np)
             head = iconvex(3,1)
             tail = iconvex(4,1)
             np = tail - head + 1
             do jpoint = 1,2
                ipoint = 1 + (jpoint-1)*(np-1)
                uvendpoint(1:2,1) = uveosedg(1:2,2*np+1-ipoint)
                uvendpoint(1:2,2) = uveosedg(1:2,ipoint)
                call add_intersection_point( &
                     uvendpoint, &
                     xyzlonglat(:,ipoint), &
                     surfpair, &
                     2, &
                     interdata_new, &
                     idendpoint(jpoint) )
             end do

             call add_intersection_curve( &
                  interdata_new, &
                  [0._fp, 0._fp, 0._fp], &
                  idendpoint, &
                  spread(spread([-1._fp-EPSuv, 1._fp+EPSuv], 2, 2), 3, 2) )
             curve_new => interdata_new%curves(interdata_new%nc)
             curve_new%smooth = .true.
             curve_new%surf = surfpair
             curve_new%isplit(2,1:2) = [1,np]
             allocate(curve_new%iedge(1))
             curve_new%iedge(:) = 0
             allocate(curve_new%polyline)
             curve_new%polyline%np = np
             allocate(curve_new%polyline%uv(2,2,np), curve_new%polyline%xyz(3,np))
             curve_new%polyline%uv(1:2,1,1:np) = uveosedg(1:2,2*np:np+1:-1)
             curve_new%polyline%uv(1:2,2,1:np) = uveosedg(1:2,1:np)
             curve_new%polyline%xyz(1:3,1:np) = &
                  0.5_fp*(xyzlonglat(1:3,1:np) + xyzlonglat(1:3,2*np:np+1:-1))

          else
             nsurf_new = nsurf_new + 1
             allocate(uvlonglat(2,nptot))
             call long_lat_patch_from_points( &
                  brep%verts(ivert)%point%xyz, &
                  xyzlonglat(1:3,1:nptot), &
                  nptot, &
                  16, & ! HARD-CODED
                  stat, &
                  surf_new(nsurf_new), &
                  uvlonglat, &
                  longlat_range, &
                  solidangle )
             IF ( STAT > 0 ) STOP 'make_partial_eos: failed to construct long-lat patch'
             vert_eos(ivert)%ptr => surf_new(nsurf_new)
             call economize2( &
                  vert_eos(ivert)%ptr%x, &
                  EPSmath )
             call compute_deriv1(vert_eos(ivert)%ptr)
             call compute_deriv2(vert_eos(ivert)%ptr)

             WRITE (STRNUM3, '(I3.3)') IVERT
             CALL WRITE_POLYNOMIAL( &
                  VERT_EOS(IVERT)%PTR%X, &
                  '../debug/eos_vert_c_'//strnum3 //'.cheb' )

             surfpair(1)%ptr => vert_eos(ivert)%ptr
             do iedge = 1,nconvex
                surfpair(2)%ptr => edge_eos(iconvex(1,iedge))%ptr
                do i = 1,2
                   j = iconvex(2+i,iedge)
                   uvendpoint(1:2,1) = uvlonglat(1:2,j)
                   uvendpoint(1:2,2) = uveosedg(1:2,j)
                   call add_intersection_point( &
                        uvendpoint, &
                        xyzlonglat(1:3,j), &
                        surfpair, &
                        2, &
                        interdata_new, &
                        idendpoint(i) ) 
                end do

                call add_intersection_curve( &
                     interdata_new, &
                     [0._fp, 0._fp, 0._fp], &
                     idendpoint, &
                     spread(spread([-1._fp-EPSuv, 1._fp+EPSuv], 2, 2), 3, 2) )

                curve_new => interdata_new%curves(interdata_new%nc)
                curve_new%smooth = .true.
                curve_new%surf = surfpair
                head = iconvex(3,iedge)
                tail = iconvex(4,iedge)
                np = tail - head + 1
                curve_new%isplit(2,1:2) = [1,np]
                allocate(curve_new%iedge(1))
                curve_new%iedge(:) = 0
                allocate(curve_new%polyline)
                curve_new%polyline%np = np
                allocate(curve_new%polyline%uv(2,2,np), curve_new%polyline%xyz(3,np))
                curve_new%polyline%uv(1:2,1,1:np) = uvlonglat(1:2,head:tail)
                curve_new%polyline%uv(1:2,2,1:np) = uveosedg(1:2,head:tail)
                curve_new%polyline%xyz(1:3,1:np) = xyzlonglat(1:3,head:tail)
             end do
             deallocate(uvlonglat)
          end if
       end if
       if ( allocated(xyzlonglat) ) deallocate(xyzlonglat)
       if ( allocated(uveosedg) ) deallocate(uveosedg)
    end do
    ! ----------<<<

  end subroutine make_partial_eos
  !-------------------------------------------------------------



  !-------------------------------------------------------------
  subroutine intersection_segment_parameter( &
       curve, &
       head, &
       tail, &
       mode, &
       w )
    use mod_types_intersection
    implicit none
    type(type_intersection_curve), intent(in), target :: curve
    integer,                       intent(in)         :: head
    integer,                       intent(in)         :: tail
    integer,                       intent(in)         :: mode
    real(kind=fp), allocatable,    intent(inout)      :: w(:)
    type(type_intersection_polyline), pointer         :: polyline => null()
    real(kind=fp)                                     :: dp(3), dpsqr, invdpsqr
    integer                                           :: step, np, i, j

    if ( .not.associated(curve%polyline) ) STOP 'intersection_segment_parameter: polyline not associated'

    polyline => curve%polyline
    if ( .not.allocated(polyline%xyz) ) STOP 'intersection_segment_parameter: xyz not allocated'
    step = sign(1,tail-head)
    np = 1 + (tail-head)/step

    if ( np < 1 ) STOP 'intersection_segment_parameter: np < 1'
    allocate(w(np))

    if ( mode > 2 ) then
       dp = polyline%xyz(1:3,tail) - polyline%xyz(1:3,head)
       dpsqr = sum(dp**2)
       if ( abs(dpsqr) < EPSmath**2 ) STOP 'intersection_segment_parameter: |dp| << 1'
       invdpsqr = 1._fp / dpsqr
    end if

    select case (mode)
    case (1) ! Hohmeyer's parameter (injectivity guaranteed)
       do i = 1,np
          j = head + (i-1)*step
          w(i) = dot_product(polyline%xyz(1:3,j), curve%param_vector)
       end do

    case (2) ! chordal approximation of arclength (injectivity guaranteed)
       w(1) = 0._fp
       w(2:np) = sqrt(sum((&
            polyline%xyz(1:3,head+step:tail:step) - &
            polyline%xyz(1:3,head:tail-step:step))**2, 1))
       do i = 2,np
          w(i) = w(i-1) + w(i)
       end do

    case (3) ! chord fraction (injectivity NOT guaranteed)
       w(1) = 0._fp
       do i = 2,np-1
          j = head + (i-1)*step
          w(i) = sqrt(invdpsqr*sum((polyline%xyz(1:3,j) - polyline%xyz(1:3,head))**2))
       end do
       w(np) = 1._fp

    case (4) ! ratio of distances from endpoints (injectivity NOT guaranteed)
       w(1) = 0._fp
       do i = 2,np-1
          j = head + (i-1)*step
          w(i) = invdpsqr * dot_product(polyline%xyz(1:3,j) - polyline%xyz(1:3,head), dp)
       end do
       w(np) = 1._fp

    end select

    ! normalize to [-1,1]
    w(1:np) = -1._fp + 2._fp*(w(1:np) - w(1))/(w(np) - w(1))

  end subroutine intersection_segment_parameter
  !-------------------------------------------------------------







  subroutine eos_from_curve2( &
       degr, &
       curve, &
       head, &
       tail, &
       enve_rl, &
       w, &
       xyzc, &
       m, &
       n )
    use mod_chebyshev
    use mod_diffgeom
    use mod_types_intersection
    use mod_intersection
    use mod_tolerances
    USE MOD_UTIL
    implicit none
    LOGICAL, PARAMETER :: DEBUG = .false.
    integer,                               intent(in)    :: degr
    type(type_intersection_curve), target, intent(in)    :: curve
    integer,                               intent(in)    :: head
    integer,                               intent(in)    :: tail
    type(ptr_surface),                     intent(in)    :: enve_rl(2)
    real(kind=fp), allocatable,            intent(inout) :: w(:)
    real(kind=fp), allocatable,            intent(inout) :: xyzc(:,:,:)
    integer,                               intent(out)   :: m, n
    type(type_intersection_polyline), pointer            :: polyline => null()
    real(kind=fp)                                        :: whead, wtail, sign_dw, wj
    real(kind=fp), allocatable                           :: wcgl(:)
    integer                                              :: step, np
    real(kind=fp), allocatable                           :: uv(:,:,:)
    real(kind=fp), allocatable                           :: xyz(:,:)
    real(kind=fp), dimension(4)                          :: lowerb, upperb
    real(kind=fp), allocatable                           :: dxyz_dw(:,:)
    real(kind=fp)                                        :: duv_ds(2,2,2), dxyz_ds(3,2)
    integer                                              :: stat
    real(kind=fp), allocatable                           :: erl(:,:,:)
    real(kind=fp)                                        :: v(degr+1)
    integer                                              :: ipt, jpt, irl

    IF ( DEBUG) PRINT *,'HEAD, TAIL =', HEAD, TAIL
    
    m = degr + 1
    
    ! Choose polynomial degree for new surface
    polyline => curve%polyline
    whead = dot_product(curve%param_vector, polyline%xyz(1:3,head))
    wtail = dot_product(curve%param_vector, polyline%xyz(1:3,tail))
    if ( allocated(w) ) deallocate(w)
    allocate(wcgl(m))
    call cgl_nodes( &
         wcgl, &
         degr, &
         whead, &
         wtail )
    IF ( DEBUG) PRINT *,'WHEAD, WTAIL =', WHEAD, WTAIL
    sign_dw = sign(1._fp, wtail - whead)
    
    ! sample intersection curve at CGL nodes
    ! (w is chosen to be Hohmeyer's parameter)
    allocate(xyz(3,m), uv(2,2,m))
    ! endpoints
    uv(1:2,1:2,1) = polyline%uv(1:2,1:2,head)
    xyz(1:3,1) = polyline%xyz(1:3,head)
    
    uv(1:2,1:2,m) = polyline%uv(1:2,1:2,tail)
    xyz(1:3,m) = polyline%xyz(1:3,tail)

    ! interior points
    step = sign(1, tail - head)
    jpt = head
    do ipt = 2,degr
       IF ( DEBUG) PRINT *,'IPT =', IPT
       do 
          wj = dot_product(curve%param_vector, polyline%xyz(1:3,jpt+step))
          IF ( DEBUG) PRINT *,'  JPT =', JPT
          IF ( DEBUG) PRINT *,'  WJ - W(IPT) =', wj - wcgl(ipt)
          if ((wj - wcgl(ipt))*sign_dw > 0._fp) then
             exit
          else
             jpt = jpt + step
             cycle
          end if
       end do
       ! first iterate
       uv(1:2,1:2,ipt) = polyline%uv(1:2,1:2,jpt)
       xyz(1:3,ipt) = polyline%xyz(1:3,jpt)

       ! refine
       lowerb = reshape(curve%uvbox(1,1:2,1:2), [4]) - EPSuv
       upperb = reshape(curve%uvbox(2,1:2,1:2), [4]) + EPSuv
       call newton_intersection_polyline2( &
            curve%surf, &
            curve%param_vector, &
            lowerb, &
            upperb, &
            wcgl(ipt), &
            EPSuv*abs(wtail - whead), &
            stat, &
            uv(1:2,1:2,ipt), &
            xyz(1:3,ipt) )
       
    end do

    !! compute tangents
    allocate(dxyz_dw(3,m))
    do ipt = 1,m
       call diffgeom_intersection( &
            curve%surf, &
            uv(1:2,1:2,ipt), &
            duv_ds, &
            dxyz_ds, &
            stat )
       IF ( STAT /= 0 ) PRINT *,'eos_from_curve2: diffgeom_intersection: STAT =',STAT
       dxyz_dw(1:3,ipt) = -dxyz_ds(1:3,1)*sign_dw
    end do

    !! evaluate right and left boundary curves
    allocate(erl(3,m,2))
    do irl = 1,2
       do ipt = 1,m
          call eval( &
               erl(1:3,ipt,irl), &
               enve_rl(irl)%ptr, &
               uv(1:2,irl,ipt) )
       end do
    end do

    !! contruct portion of canal surfaces between right/left boundary curves
    np = 1 + (tail - head)/step
    allocate(w(np))
    do ipt = head,tail,step
       w(ipt) = dot_product(curve%param_vector, polyline%xyz(1:3,ipt))
    end do
    n = m
    allocate(xyzc(m,n,3))

    call cgl_nodes( &
         v, &
         degr, &
         1._fp, &
         0._fp)
    call eos_from_polyline( &
         xyz, &
         dxyz_dw, &
         erl(1:3,1:m,1), &
         erl(1:3,1:m,2), &
         m, &
         v, &
         n, &
         xyzc )

    xyzc(1:m,1:n,1:3) = xyzc(m:1:-1,1:n,1:3)

    IF ( DEBUG ) THEN
       CALL WRITE_MATRIX( &
            TRANSPOSE(ERL(1:3,1:M,1)), &
            M, 3, &
            '../debug/eos_from_curve2_erl_1.dat' )
       CALL WRITE_MATRIX( &
            TRANSPOSE(ERL(1:3,1:M,2)), &
            M, 3, &
            '../debug/eos_from_curve2_erl_2.dat' )
       CALL WRITE_MATRIX( &
            TRANSPOSE(XYZ(1:3,1:M)), &
            M, 3, &
            '../debug/eos_from_curve2_xyz.dat' )
       CALL WRITE_MATRIX( &
            TRANSPOSE(DXYZ_DW(1:3,1:M)), &
            M, 3, &
            '../debug/eos_from_curve2_dxyz_dw.dat' )
       CALL WRITE_MATRIX( &
            [W], &
            M, 1, &
            '../debug/eos_from_curve2_w.dat' )
       PAUSE
    END IF

    ! map w to [-1,1]
    w = -1._fp + 2._fp*(w - whead)/(wtail - whead)

    deallocate(wcgl, uv, xyz, dxyz_dw)
    
  end subroutine eos_from_curve2




  subroutine newton_intersection_polyline2( &
       surf, &
       param_vector, &
       lowerb, &
       upperb, &
       wtarget, &
       tolw, &
       stat, &
       uv, &
       xyz )
    use mod_math
    use mod_linalg
    use mod_diffgeom
    use mod_tolerances
    implicit none
    LOGICAL, PARAMETER :: DEBUG = .false.
    integer,           parameter     :: itmax = 2 + ceiling(-log10(EPSuv))
    type(ptr_surface), intent(in)    :: surf(2)
    real(kind=fp),     intent(in)    :: param_vector(3)
    real(kind=fp),     intent(in)    :: lowerb(4)
    real(kind=fp),     intent(in)    :: upperb(4)
    real(kind=fp),     intent(in)    :: wtarget
    real(kind=fp),     intent(in)    :: tolw
    integer,           intent(out)   :: stat
    real(kind=fp),     intent(inout) :: uv(2,2)
    real(kind=fp),     intent(out)   :: xyz(3)
    real(kind=fp)                    :: xyz_tmp(3,2)
    real(kind=fp)                    :: resw, resxyz(3)
    real(kind=fp)                    :: jac(4,4), duv(4)
    real(kind=fp)                    :: cond, erruv
    integer                          :: it, isurf, ivar
    integer                          :: stat_refl

    IF ( DEBUG ) THEN
       PRINT *,''; PRINT *,'';
       PRINT *,'NEWTON_INTERSECTION_POLYLINE2'
       PRINT *,'WTARGET =',wtarget
       PRINT *,'LOWERB =',LOWERB
       PRINT *,'UPPERB =',UPPERB
    END IF
    
    stat = 1
    erruv = 0._fp
    cond = 1._fp

    do it = 1,itmax
       !! compute residual
       do isurf = 1,2
          call eval( &
               xyz_tmp(:,isurf), &
               surf(isurf)%ptr, &
               uv(:,isurf) )
       end do

       resxyz = xyz_tmp(1:3,2) - xyz_tmp(1:3,1)
       resw = wtarget - dot_product(param_vector, xyz_tmp(1:3,1))

       IF ( DEBUG ) THEN
          PRINT *,'RESXYZ =',NORM2(RESXYZ), ', RESW =',ABS(RESW)
       END IF

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
        jac(4,ivar) = dot_product(param_vector, jac(1:3,ivar))
     end do
     jac(4,3:4) = 0._fp

     !! solve for Newton step
     call linsolve_svd( &
          duv, &
          jac, &
          [resxyz, resw], &
          4, &
          4, &
          1, &
          cond )
     erruv = max(sum(duv(1:2)**2), sum(duv(3:4)**2))

     !! correct Newton step to keep the iterate inside feasible region
     call nd_box_reflexions( &
          reshape(uv, [4]), &
          lowerb, &
          upperb, &
          duv, &
          4, &
          stat_refl )
     if ( stat_refl > 0 ) return

     !! update solution
     uv(:,1) = uv(:,1) + duv(1:2)
     uv(:,2) = uv(:,2) + duv(3:4)

     !! termination criteria
     if ( erruv < max(EPSuvsqr, EPSfpsqr*cond**2) ) then ! <--------+
        if ( sum(resxyz**2) < EPSxyzsqr .and. &                     !
             abs(resw) < tolw ) then ! <-------------------------+  !
           stat = 0                                              !  !
           do isurf = 1,2 ! <------------+                       !  !
              call eval( &               !                       !  !
                   xyz_tmp(:,isurf), &   !                       !  !
                   surf(isurf)%ptr, &    !                       !  !
                   uv(:,isurf) )         !                       !  !
           end do ! <--------------------+                       !  !
           xyz = 0.5_fp * sum(xyz_tmp, 2)                        !  !
        else ! --------------------------------------------------+  !
        end if ! ------------------------------------------------+  !
        return                                                      !
     end if ! <-----------------------------------------------------+

  end do
  
  end subroutine newton_intersection_polyline2







   !subroutine long_lat_patch_from_arcs( &
   !   vxyz, &
   !   m, &
   !   cxyz, &
   !   occ, &
   !   tng, &
   !   surf, &
   !   longlat_ctr, &
   !   longlat_rng )
   !   use mod_diffgeom
   !   use mod_geometry
   !   implicit none
   !   real(kind=fp),      intent(in)    :: vxyz(3)
   !   integer,            intent(in)    :: m
   !   real(kind=fp),      intent(in)    :: cxyz(3,m)
   !   real(kind=fp),      intent(in)    :: occ(3,m)
   !   real(kind=fp),      intent(in)    :: tng(3,m)
   !   type(type_surface), intent(inout) :: surf
   !   real(kind=fp),      intent(out)   :: uv(2,n)
   !   real(kind=fp),      intent(out)   :: longlat_ctr(2)
   !   real(kind=fp),      intent(out)   :: longlat_rng(2)
   !   real(kind=fp)                     :: ctr(3), radsqr
   !   real(kind=fp)                     :: rot(3,3)
   !  
   !   ! 1) get axis of smallest cone bounding the points cxyz
   !   !   a) get smallest sphere bounding the points cxyz
   !   call smallest_enclosing_ball( &
   !      m, &
   !      3, &
   !      cxyz - spread(vxz, dim=2, ncopies=m), &
   !      ctr, &
   !      radsqr )
   !   !   b) get cone axis
   !   rot(1:3,1) = ctr/norm2(ctr)  
   !  
   !end subroutine long_lat_patch_from_arcs
end module mod_eos
