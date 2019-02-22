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
    normal_speed = 0.15d0*(1.d0 + 0.15d0*cos(5.d0*(x+y+z)))

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
    integer                                                :: isurf, jsurf
    integer                                                :: icurv, jcurv
    integer                                                :: ipoint, jpoint
    integer                                                :: ieos
    integer                                                :: endpts(2)
    type(type_intersection_curve), pointer                 :: curve => null(), curve_new => null()
    type(ptr_surface)                                      :: enve_rl(2), surfpair(2)
    integer                                                :: m, n, np
    real(kind=fp), allocatable                             :: xyz(:,:,:)
    real(kind=fp)                                          :: uvend(2,2)
    

    nsurf_new = nsurf_old
    do icurv = 1,interdata_old%nc
       if ( curve_type(icurv) == 2 ) nsurf_new = nsurf_new + 1
    end do

    allocate(surf_new(nsurf_new))

    ! EoS from surfaces - - -
    do isurf = 1,nsurf_old
       call eos_from_surface( &
            timestep, &
            surf_old(isurf), &
            surf_new(isurf) )
    end do
    ! - - - - - - - -

    
    ! EoS from (convex) curves - - -
    ieos = nsurf_old
    do icurv = 1,interdata_old%nc
       curve => interdata_old%curves(icurv)

       do jsurf = 1,2
          do isurf = 1,nsurf_old
             if ( associated(curve%surf(jsurf)%ptr, surf_old(isurf)) ) then
                enve_rl(jsurf)%ptr => surf_new(isurf)
                exit
             end if
          end do
       end do
       
       if ( curve_type(icurv) == 0 ) then
          do jpoint = 1,2
             ipoint = 1 + (jpoint-1)*(curve%polyline%np-1)
             call add_intersection_point( &
                  curve%polyline%uv(:,:,ipoint), &
                  curve%polyline%xyz(:,ipoint), &
                  curve%surf, &
                  2, &
                  interdata_new, &
                  endpts(jpoint) )
          end do
          
          call add_intersection_curve( &
               interdata_new, &
               [0._fp, 0._fp, 0._fp], &
               endpts, &
               interdata_old%curves(icurv)%uvbox )

          curve_new => interdata_new%curves(interdata_new%nc)
          allocate(curve_new%polyline)
          np = curve%polyline%np
          curve_new%polyline%np = np
          allocate(curve_new%polyline%uv(2,2,np))
          curve_new%polyline%uv(1:2,1:2,1:np) = curve%polyline%uv(1:2,1:2,1:np)
          allocate(curve_new%polyline%xyz(3,np))
          do ipoint = 1,np
             call eval( &
                  curve_new%polyline%xyz(1:3,ipoint), &
                  surfpair(1)%ptr, &
                  curve_new%polyline%uv(1:2,2,ipoint) )
          end do

          curve_new%smooth = curve%smooth
          allocate(curve_new%iedge(curve%nsplit-1))
          curve_new%iedge = curve%iedge
          cycle
       elseif ( curve_type(icurv) == 1 ) then
          cycle
       end if

       ieos = ieos + 1

       call eos_from_curve( &
          curve%polyline, &
          1, &
          curve%polyline%np, &
          enve_rl, &
          xyz, &
          m, &
          n )

       call reset_polynomial(surf_new(ieos)%x, 2, 1, [m-1,n-1], 3)
       call fcht2( &
          xyz, &
          surf_new(ieos)%x%coef(1:m,1:n,1:3), &
          m, &
          n, &
          3, &
          EPSmath )


       ! 1st curve (R-E), 2nd curve (E-L)
       do jcurv = 1,2
          surfpair(jcurv)%ptr => curve%surf(jcurv)%ptr
          surfpair(1+mod(jcurv,2))%ptr => surf_new(ieos)

          do jpoint = 1,2
             ipoint = 1 + (jpoint-1)*(curve%polyline%np-1)
             uvend(1:2,jcurv) = curve%polyline%uv(:,1,ipoint)
             uvend(1,1+mod(jcurv,2)) = real((-1)**(ipoint+jcurv), kind=fp)
             uvend(2,1+mod(jcurv,2)) = real((-1)**jcurv, kind=fp)
             call add_intersection_point( &
                  uvend, &
                  curve%polyline%xyz(:,ipoint), &
                  surfpair, &
                  2, &
                  interdata_new, &
                  endpts(jpoint) )
          end do

          call add_intersection_curve( &
               interdata_new, &
               [0._fp, 0._fp, 0._fp], &
               endpts, &
               interdata_old%curves(icurv)%uvbox )
          curve_new => interdata_new%curves(interdata_new%nc)
          curve_new%uvbox(1,1:2,1+mod(jcurv,2)) = -1._fp
          curve_new%uvbox(2,1:2,1+mod(jcurv,2)) = 1._fp
          allocate(curve_new%polyline)
          np = curve%polyline%np
          curve_new%polyline%np = np
          allocate(curve_new%polyline%uv(2,2,np))
          curve_new%polyline%uv(1:2,jcurv,1:np) = curve%polyline%uv(1:2,jcurv,1:np)
          curve_new%polyline%uv(1,1+mod(jcurv,2),1:np) = &
               curve_new%polyline%s(1:np)*real((-1)**jcurv, kind=fp)
          curve_new%polyline%uv(2,1+mod(jcurv,2),1:np) = real((-1)**jcurv, kind=fp)
          ! update xyz
          allocate(curve_new%polyline%xyz(3,np))
          do ipoint = 1,np
             call eval( &
                  curve_new%polyline%xyz(1:3,ipoint), &
                  surfpair(1)%ptr, &
                  curve_new%polyline%uv(1:2,2,ipoint) )
          end do
       end do
       ! ---------------

       
    end do
    ! - - - - - - - -
    
    
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
    real(kind=fp), dimension(3)                :: bw, ewr, ewl, uw, vw, gw
    real(kind=fp)                              :: wr, wl, w, a, rwr, rwl, rw, aj, tj
    integer                                    :: i, j

    do i = 1,m
       ewr = er(1:3,i) - g(1:3,i)
       ewl = el(1:3,i) - g(1:3,i)

       bw = dg(1:3,i)/norm2(dg(1:3,i))
       wr = dot_product(ewr, bw)
       wl = dot_product(ewl, bw)
       w = 0.5_fp * (wr + wl)

       ewr = ewr + wr*bw
       ewl = ewl + wl*bw

       a = angle_between_vectors_3d(ewr, ewl)

       rwr = norm2(ewr)
       rwl = norm2(ewl)

       uw = ewr/rwr
       vw = ewl - dot_product(ewl,uw)*uw
       vw = vw/norm2(vw)

       gw = g(1:3,i) + w*bw
       do j = 1,n
          tj = 0.5_fp*(v(j) + 1._fp)
          aj = tj*a
          rw = (1._fp - tj)*rwr + tj*rwl
          w = (1._fp - tj)*wr + tj*wl
          gw = g(1:3,i) + w*bw
          e(i,j,1:3) = gw + rw*(uw*cos(aj) + vw*sin(aj))
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
    real(kind=fp), allocatable,       intent(inout) :: xyzc(:,:,:)
    integer,                          intent(out)   :: m, n
    integer                                         :: step, np, degr
    real(kind=fp), allocatable                      :: x(:)
    real(kind=fp), allocatable                      :: y(:,:)
    real(kind=fp), allocatable                      :: c(:,:), d(:,:)
    real(kind=fp)                                   :: cond, errL2, errLinf
    real(kind=fp), allocatable                      :: g(:,:), dg(:,:)
    real(kind=fp), allocatable                      :: uv(:,:), erl(:,:,:)
    integer                                         :: irl, ipt

    !! fit Chebyshev polynomial to intersection curve
    step = sign(1,tail-head)
    np = 1 + (tail-head)/step
    degr = max(1, int(0.25*real(np)) - 1)
    m = degr+1
    allocate(c(m,7))

    allocate(x(np))
    x(1) = 0._fp
    x(2:np) = sqrt(sum((&
         polyline%xyz(1:3,head+step:tail:step) - &
         polyline%xyz(1:3,head:tail-step:step))**2, 1))
    do ipt = 2,np
       x(ipt) = x(ipt-1) + x(ipt)
    end do
    x = 2._fp*x/x(np) - 1._fp

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

    PRINT *,'COND =', COND, 'ERR_(L2,LINF) =', ERRL2, ERRLINF
    deallocate(x,y)

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
    allocate(x(n), xyzc(m,n,3))
    call cgl_nodes( &
         x, &
         degr, &
         1._fp, &
         -1._fp )

    call eos_from_polyline( &
         transpose(g), &
         transpose(dg), &
         transpose(erl(1:m,1:3,1)), &
         transpose(erl(1:m,1:3,2)), &
         m, &
         x, &
         n, &
         xyzc )

  end subroutine eos_from_curve
  !-------------------------------------------------------------




  
end module mod_eos
