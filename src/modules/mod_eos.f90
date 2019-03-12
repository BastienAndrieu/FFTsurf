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
    normal_speed = 0.15d0*(1.d0 + 0.05d0*cos(5.d0*(x+y+z)))

  end function normal_speed
  !-------------------------------------------------------------


  
  
  !-------------------------------------------------------------
  subroutine make_all_eos( &
       timestep, &
       brep, &
       surf_old, &
       nsurf_old, &
       interdata_old, &
       curve_type, &
       surf_new, &
       nsurf_new, &
       interdata_new )
    use mod_diffgeom
    use mod_types_brep
    use mod_types_intersection
    use mod_intersection
    implicit none
    real(kind=fp),                   intent(in)            :: timestep
    type(type_brep),                 intent(in)            :: brep
    integer,                         intent(in)            :: nsurf_old
    type(type_surface),              intent(in), target    :: surf_old(nsurf_old)
    type(type_intersection_data),    intent(in), target    :: interdata_old
    integer,                         intent(in)            :: curve_type(interdata_old%nc)
    integer,                         intent(out)           :: nsurf_new
    type(type_surface), allocatable, intent(inout), target :: surf_new(:)
    type(type_intersection_data),    intent(out), target   :: interdata_new
    integer                                                :: icurv
    integer                                                :: isurf, jsurf
    type(type_intersection_curve), pointer                 :: curve => null()
    integer                                                :: curve2surfid(2,interdata_old%nc) 

    ! get INDICES of surfaces incident to each intersection curve
    curve2surfid(:,:) = 0
    do icurv = 1,interdata_old%nc
       curve => interdata_old%curves(icurv)
       do jsurf = 1,2
          do isurf = 1,nsurf_old
             if ( associated(curve%surf(jsurf)%ptr, surf_old(isurf)) ) then
                curve2surfid(jsurf,icurv) = isurf
                exit
             end if
          end do
       end do
    end do
    ! - - - - - - - - - - - - - - - - - - - - - - - -

  end subroutine make_all_eos
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
       uvxyz, &
       np )
    use mod_diffgeom
    implicit none
    type(type_surface),         intent(in)    :: surf
    integer,                    intent(in)    :: iborder
    real(kind=fp),              intent(in)    :: tolchord
    real(kind=fp),              intent(in)    :: hmin, hmax
    integer,                    intent(out)   :: stat
    real(kind=fp), allocatable, intent(inout) :: uvxyz(:,:)
    integer,                    intent(out)   :: np
    real(kind=fp)                             :: uvends(2,2)
    real(kind=fp)                             :: h_endpoints(2)
    
    select case (iborder)
    case (1)
       uvends(1,1:2) = [-1._fp, 1._fp]
       uvends(2,1:2) = -1._fp
    case (2)
       uvends(1,1:2) = 1._fp  
       uvends(2,1:2) = [-1._fp, 1._fp]
    case (3)
       uvends(1,1:2) = [1._fp, -1._fp]
       uvends(2,1:2) = 1._fp
    case (4)
       uvends(1,1:2) = -1._fp  
       uvends(2,1:2) = [1._fp, -1._fp]
    end select

    !(...)
    
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
    implicit none
    integer,                       intent(in)  :: m
    real(kind=fp), dimension(3,m), intent(in)  :: g, dg, er, el
    integer,                       intent(in)  :: n
    real(kind=fp),                 intent(in)  :: v(n)
    real(kind=fp),                 intent(out) :: e(m,n,3)
    real(kind=fp), dimension(3,m)              :: o, eor, eol
    real(kind=fp)                              :: ei(3,n)
    integer                                    :: i, j

    o = g + dg*spread(0.5_fp*sum((er + el - 2._fp*g)*dg,1)/sum(dg**2,1), dim=1, ncopies=3)

    eor = er - o
    eol = el - o

    do i = 1,m      
       call slerp( &
            eor(1:3,i), &
            eol(1:3,i), &
            v, &
            n, &
            ei )

       do j = 1,n
          e(i,j,1:3) = o(1:3,i) + ei(1:3,j)
       end do
    end do

  end subroutine eos_from_polyline
  !-------------------------------------------------------------


  
  

  !-------------------------------------------------------------
  subroutine eos_from_curve( &
       polyline, &
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
    !real(kind=fp), parameter                        :: FITTOL_uv = 1.d-7
    !real(kind=fp), parameter                        :: FITTOL_xyz = 1.d-7
    type(type_intersection_polyline), intent(in)    :: polyline
    integer,                          intent(in)    :: head
    integer,                          intent(in)    :: tail
    type(ptr_surface),                intent(in)    :: enve_rl(2)
    real(kind=fp),                    intent(out)   :: x(polyline%np)
    real(kind=fp), allocatable,       intent(inout) :: xyzc(:,:,:)
    integer,                          intent(out)   :: m, n
    integer                                         :: step, np, degr
    real(kind=fp)                                   :: y(polyline%np,7)
    real(kind=fp), allocatable                      :: c(:,:), d(:,:)
    real(kind=fp)                                   :: cond, errL2, errLinf
    real(kind=fp), allocatable                      :: u(:)
    real(kind=fp), allocatable                      :: g(:,:), dg(:,:)
    real(kind=fp), allocatable                      :: uv(:,:), erl(:,:,:)
    integer                                         :: irl, ipt

    !! fit Chebyshev polynomial to intersection curve
    step = sign(1,tail-head)
    np = 1 + (tail-head)/step
    degr = max(1, int(0.25*real(np)) - 1)
    m = degr+1
    allocate(c(m,7))

    x(1) = 0._fp
    x(2:np) = sqrt(sum((&
         polyline%xyz(1:3,head+step:tail:step) - &
         polyline%xyz(1:3,head:tail-step:step))**2, 1))
    do ipt = 2,np
       x(ipt) = x(ipt-1) + x(ipt)
    end do
    x = 2._fp*x/x(np) - 1._fp

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
    LOGICAL :: DEBUG = .FALSE.
    real(kind=fp), parameter          :: mrg = 0.1_fp
    real(kind=fp),      intent(in)    :: center(3)
    integer,            intent(in)    :: n
    real(kind=fp),      intent(inout) :: xyz(3,n)
    integer,            intent(in)    :: degr
    integer,            intent(out)   :: stat
    type(type_surface), intent(inout) :: surf
    real(kind=fp),      intent(out)   :: uv(2,n)
    real(kind=fp),      intent(out)   :: longlat_range(2)
    real(kind=fp),      intent(out)   :: solidangle
    real(kind=fp)                     :: xyzavg(3), ravg, r
    real(kind=fp)                     :: axes(3,3)
    real(kind=fp)                     :: xyzproj(3,n), dot
    real(kind=fp)                     :: boxcenter(2), boxranges(2), boxaxes(3,3)
    real(kind=fp), dimension(2)       :: minuv, maxuv
    real(kind=fp), dimension(degr+1)  :: tcgl, cosucgl, sinucgl, cosvcgl, sinvcgl
    real(kind=fp)                     :: xyzcgl(3,(degr+1)**2)
    integer                           :: i, j, k, m
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
       IF ( DEBUG ) PRINT *,'R =',R
       xyz(1:3,i) = xyz(1:3,i)/r
       xyzavg(1:3) = xyzavg(1:3) + xyz(1:3,i)
       ravg = ravg + r
    end do
    xyzavg = xyzavg/norm2(xyzavg)!real(n, kind=fp)
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
       if ( dot_product(xyz(1:3,i), xyzavg) < EPSfp ) then
          stat = 1 
          return
       else
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
         '/d/bandrieu/GitHub/FFTsurf/test/demo_EoS_brep/debug/debuglonglat_uv.dat')

    ! longitude-latitude bounding box
    minuv = minval(uv, dim=2)
    maxuv = maxval(uv, dim=2)
    longlat_range = (1._fp + mrg)*max(abs(maxuv), abs(minuv))
    longlat_range = longlat_range
    do i = 1,n
       uv(1:2,i) = uv(1:2,i)/longlat_range
    end do
    solidangle = longlat_range(1) * sin(longlat_range(2)) / CSTpi

    ! Chebyshev patch    
    call cgl_nodes(tcgl, degr)
    cosucgl = cos(longlat_range(1)*tcgl)
    sinucgl = sin(longlat_range(1)*tcgl)
    cosvcgl = cos(longlat_range(2)*tcgl)
    sinvcgl = sin(longlat_range(2)*tcgl)
    
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

    do i = 1,n
       xyz(1:3,i) = center(1:3) + ravg*xyz(1:3,i)
    end do
        
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

  
end module mod_eos
