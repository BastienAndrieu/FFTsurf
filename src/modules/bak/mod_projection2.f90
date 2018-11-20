module mod_projection

  use mod_math

  implicit none

contains


subroutine projection_hyperface( &
     brep, &
     iface, &
     uv, &
     xyz, &
     duv, &
     dxyz, &
     ifacetmp, &
     uvtmp, &
     debug, &
     stat_proj )
  use mod_types_brep
  use mod_types_intersection
  use mod_brep
  implicit none
  logical, intent(in) :: debug
  real(kind=fp), parameter                  :: lambmax = 1.0_fp
  type(type_brep), intent(in)               :: brep
  integer,         intent(in)               :: iface
  real(kind=fp),   intent(in)               :: uv(2)
  real(kind=fp),   intent(in)               :: xyz(3)
  real(kind=fp),   intent(in)               :: duv(2)
  real(kind=fp),   intent(in)               :: dxyz(3)
  integer,         intent(out)              :: ifacetmp
  real(kind=fp),   intent(out)              :: uvtmp(2)
  integer,         intent(out)              :: stat_proj
  real(kind=fp), dimension(4)               :: lowerb, upperb
  real(kind=fp)                             :: duvtmp(2), xyztmp(3), dxyztmp(3)
  real(kind=fp)                             :: xyztarget(3)
  logical                                   :: change_face
  integer                                   :: iwire, ihedg(2), istart
  type(type_intersection_curve), pointer    :: curve => null()
  type(type_matrix)                         :: uvpoly(3)
  real(kind=fp), allocatable                :: xyzpoly(:,:)
  integer                                   :: ihead, itail, np, sens
  integer                                   :: ninter
  integer, allocatable                      :: ipls(:,:)
  real(kind=fp), allocatable                :: lamb(:,:)
    
  IF ( DEBUG ) PRINT *,'PROJECTION HYPERFACE'

  upperb(1:4) = 1._fp + EPSuv
  lowerb = -upperb

  ifacetmp = iface
  uvtmp = uv
  stat_proj = 1
  
  duvtmp = duv
  xyztmp = xyz
  dxyztmp = dxyz
  
  xyztarget = xyz + dxyz

  allocate(uvpoly(1)%mat(2,2))

  faces : do
     IF ( DEBUG ) PRINT *,'FACE',ifacetmp
     ! the provisional point uvtmp is inside the parametric domain of
     ! the face ifacetmp (possibly on its boundary)
     ! check whether the point uvtmp + duvtmp is inside/outside this domain
     uvpoly(1)%mat(:,1) = uvtmp
     uvpoly(1)%mat(:,2) = uvtmp + lambmax*duvtmp

     change_face = .false.

     ! to do so, intersect the segment [uvtmp, uvtmp+duvtmp] with all the curve segments
     ! that form the boundary of the domain (outer and inner wires of the face)
     wires : do iwire = 0,brep%faces(ifacetmp)%ninner
          ! get index of first halfedge on the wire
          if ( iwire == 0 ) then ! <---------------------+
             ihedg = brep%faces(ifacetmp)%outer          !
          else ! ----------------------------------------+
             ihedg = brep%faces(ifacetmp)%inner(:,iwire) !
          end if ! <-------------------------------------+
          !
          istart = get_orig(brep, ihedg) ! first brep vertex on the wire
          halfedges : do
             IF ( DEBUG ) PRINT *,'HALFEDGE',ihedg
             curve => brep%edges(ihedg(1))%curve
             ! get current halfedge's uv polyline, oriented such that
             ! the incident face's interior is on the left
             call get_polyline_endpoints( &
                  brep, &
                  ihedg, &
                  ihead, &
                  itail, &
                  sens, &
                  np )
             PRINT *,'SENS =',sens
             allocate(uvpoly(2)%mat(2,np), uvpoly(3)%mat(2,np), xyzpoly(3,np))
             uvpoly(2)%mat(1:2,1:np) = curve%polyline%uv(1:2,sens,ihead:itail:(-1)**sens)
             uvpoly(3)%mat(1:2,1:np) = curve%polyline%uv(1:2,ihedg(2),ihead:itail:(-1)**sens)
             xyzpoly(1:3,1:np) = curve%polyline%xyz(1:3,ihead:itail:(-1)**sens)
             !
             ! compute intersection between the polyline and the uv displacement
             ninter = 0
             call intersect_2Dpolylines( &
                  uvpoly(1:2), &
                  [1,1], &
                  [2,np], &
                  ninter, &
                  ipls, &
                  lamb )
             !
             ! renormalize lambda
             lamb(1,:) = lamb(1,:) / lambmax
             !
             IF ( DEBUG ) THEN
                call debug_phf_polylines( &
                     uvpoly(1:2), &
                     ninter, &
                     ipls, &
                     lamb )
                PRINT *,'--> CHECK POLYLINES'
                PAUSE
             END IF
             !
             stat_proj = 0 !******
             return   !******
             !
          end do halfedges
       end do wires
  end do faces
  
end subroutine projection_hyperface




subroutine projection_hyperedge( &
       brep, &
       hyperedge, &
       iedge, &
       uv, &
       xyz, &
       duv, &
       dxyz, &
       iedgetmp, &
       uvtmp, &
       debug, &
       stat_proj )
    use mod_linalg
    use mod_types_brep
    use mod_brep
    use mod_halfedge
    use mod_types_intersection
    use mod_intersection
    use mod_tolerances
    use mod_hypergraph
    implicit none
    logical, intent(in) :: debug
    type(type_brep),      intent(in), target  :: brep
    type(type_hyperedge), intent(in)          :: hyperedge
    integer,              intent(in)          :: iedge
    real(kind=fp),        intent(in)          :: uv(2,2)
    real(kind=fp),        intent(in)          :: xyz(3)
    real(kind=fp),        intent(in)          :: duv(2,2)
    real(kind=fp),        intent(in)          :: dxyz(3)
    integer,              intent(out)         :: iedgetmp
    real(kind=fp),        intent(out)         :: uvtmp(2,2)
    integer,              intent(out)         :: stat_proj
    real(kind=fp)                             :: duvtmp(2,2)
    real(kind=fp), dimension(3)               :: xyzprev, xyztmp, xyztarget, dxyztmp, p
    type(type_intersection_curve), pointer    :: curve => null()
    type(type_intersection_polyline), pointer :: polyline => null()
    real(kind=fp)                             :: w, wfirst, wlast, ds, duv_ds(2,2,2), dxyz_ds(3,2)
    integer                                   :: jedge, ihedg(2), ifirst, ilast, sens, np, k, stat
    real(kind=fp), dimension(4)               :: lowerb, upperb

    IF ( DEBUG ) PRINT *,'PROJECTION HYPEREDGE'
    stat_proj = 1
    upperb(:) = 2._fp
    lowerb = -upperb
    
    do jedge = 1,hyperedge%ne
       if ( hyperedge%halfedges(1,jedge) == iedge ) then
          ihedg = hyperedge%halfedges(:,jedge)
          iedgetmp = ihedg(1)
          curve => brep%edges(iedgetmp)%curve
          polyline => curve%polyline
          call get_polyline_endpoints( &
               brep, &
               ihedg, &
               ifirst, &
               ilast, &
               sens, &
               np )
          exit
       end if
    end do

    uvtmp = uv
    duvtmp = duv
    xyztmp = xyz
    dxyztmp = dxyz
    
    do k = 1,hyperedge%ne
       IF ( DEBUG ) THEN
          PRINT *,'IHEDGE =', IHEDG
          PRINT *,' FACES =', get_face(brep, [iedgetmp,1]), get_face(brep, [iedgetmp,2])
          PRINT *,'    UV =', UVtmp
          PRINT *,'   XYZ =', XYZtmp
          PRINT *,'   DUV =', DUVtmp
          PRINT *,'  DXYZ =', DXYZtmp
       END IF

       xyztarget = xyztmp + dxyztmp
       
       p = real((-1)**sens, kind=fp) * curve%param_vector
       wfirst = dot_product(p, polyline%xyz(:,ifirst))
       wlast  = dot_product(p, polyline%xyz(:,ilast))
       w = dot_product(p, xyztarget)
       IF ( DEBUG ) THEN
          PRINT *,'WFIRST, WLAST, W =',wfirst, wlast, w
       END IF
       
       if ( is_in_closed_interval(w, wfirst, wlast, tolerance=EPSuv) ) then
          ! relax temporary point onto true intersection underlying the hyperedge
          uvtmp = uvtmp + duvtmp ! initial iterate
          xyzprev = xyztmp
          call newton_intersection_polyline( &
               curve%surf, &
               lowerb, &
               upperb, &
               xyzprev, &
               sum(dxyztmp**2), &
               stat, &
               uvtmp, &
               xyztmp )
          if ( stat > 0 ) then
             PRINT *,'projection_hyperedge : failed to reproject onto transversal intersection curve, k =',k
             PRINT *,'POLYLINE ENDPOINTS:'
             PRINT *,polyline%xyz(:,ifirst)
             PRINT *,polyline%xyz(:,ilast)
             !PAUSE
             RETURN
          else
             IF ( DEBUG ) THEN
                PRINT *,'CONVERGED :'
                PRINT *,'IHEDGE =', IHEDG
                PRINT *,' FACES =', get_face(brep, [iedgetmp,1]), get_face(brep, [iedgetmp,2])
                PRINT *,'    UV =', UVTMP
                PRINT *,'   XYZ =', XYZTMP
             END IF
             stat_proj = 0
             return
          end if
       else
          ! the target point is on another halfedge
          if ( w < wfirst ) then
             ! move on to previous halfedge on the hyperedge
             jedge = 1 + mod(jedge + hyperedge%ne - 2,hyperedge%ne) ! previous
          else
             ! move on to next halfedge on the hyperedge
             jedge = 1 + mod(jedge,hyperedge%ne) ! next
          end if

          ihedg = hyperedge%halfedges(:,jedge)
          iedgetmp = ihedg(1)
          curve => brep%edges(iedgetmp)%curve
          polyline => curve%polyline
          call get_polyline_endpoints( &
               brep, &
               ihedg, &
               ifirst, &
               ilast, &
               sens, &
               np )

          if ( w < wfirst ) then
             uvtmp = polyline%uv(:,:,ilast)
             xyztmp = polyline%xyz(:,ilast)
          else
             uvtmp = polyline%uv(:,:,ifirst)
             xyztmp = polyline%xyz(:,ifirst)
          end if

          dxyztmp = xyztarget - xyztmp
          call diffgeom_intersection( &
               curve%surf, &
               uvtmp, &
               duv_ds, &
               dxyz_ds, &
               stat )
          ds = dot_product(dxyztmp, dxyz_ds(:,1))
          duvtmp = ds * duv_ds(:,1,:)
          dxyztmp = ds * dxyz_ds(:,1)

       end if
    end do
    
    uvtmp = uv
    iedgetmp = iedge
    PRINT *,'projection_hyperedge : failed to reproject'
    RETURN
    
  end subroutine projection_hyperedge





  subroutine insert_point_in_curve( &
       curve, &
       uvinter, &
       xyzinter, &
       iprev, &
       debug, &
       stat )
    use mod_types_intersection
    use mod_intersection
    implicit none
    type(type_intersection_curve), intent(inout) :: curve
    real(kind=fp),                 intent(in)    :: uvinter(2,2)
    real(kind=fp),                 intent(in)    :: xyzinter(3)
    integer,                       intent(in)    :: iprev
    logical,                       intent(in)    :: debug
    integer,                       intent(out)   :: stat
    
    ! add polyline point
    call insert_polyline_point( &
         uvinter, &
         xyzinter, &
         stat, &
         curve%polyline, &
         iprev )
    if ( stat == 0 ) then
       IF ( DEBUG ) PRINT *,'NEW POINT INSERTED IN POLYLINE : UV =',uvinter
       ! shift downstream polyline split points
       where ( curve%isplit(2,:) > iprev ) &
            curve%isplit(2,:) = curve%isplit(2,:) + 1
    end if
    
  end subroutine insert_point_in_curve





  

  subroutine debug_phf_polylines( &
       uvpoly, &
       ninter, &
       ipls, &
       lamb )
    implicit none
    type(type_matrix), intent(in) :: uvpoly(2)
    integer,           intent(in) :: ninter
    integer,           intent(in) :: ipls(:,:)
    real(kind=fp),     intent(in) :: lamb(:,:)
    integer                       :: i

    OPEN(UNIT=13, FILE='/d/bandrieu/GitHub/FFTsurf/test/polylines/xy1.dat', ACTION='WRITE')
    WRITE (13,*) 2
    DO I = 1,2
       WRITE (13,*) UVPOLY(1)%MAT(:,I)
    END DO
    CLOSE(13)
    OPEN(UNIT=13, FILE='/d/bandrieu/GitHub/FFTsurf/test/polylines/xy2.dat', ACTION='WRITE')
    WRITE (13,*) size(UVPOLY(2)%mat,2)
    DO I = 1,size(UVPOLY(2)%mat,2)
       WRITE (13,*) UVPOLY(2)%MAT(:,I)
    END DO
    CLOSE(13)
    OPEN(UNIT=13, FILE='/d/bandrieu/GitHub/FFTsurf/test/polylines/result.dat', ACTION='WRITE')
    WRITE (13,*) NINTER
    DO I = 1,NINTER
       WRITE (13,*) IPLS(:,I)
    END DO
    DO I = 1,NINTER
       WRITE (13,*) LAMB(:,I)
    END DO
    CLOSE(13)
    
  end subroutine debug_phf_polylines


  subroutine debug_phf_diffgeom( &
       curve, &
       uvinter )
    use mod_polynomial
    use mod_types_intersection
    implicit none
    type(type_intersection_curve), intent(in) :: curve
    real(kind=fp),                 intent(in) :: uvinter(2,2)
    
    CALL WRITE_POLYNOMIAL(curve%surf(1)%ptr%x, '/d/bandrieu/GitHub/FFTsurf/test/Jouke/debug_diffgeominter_surf1.cheb')
    CALL WRITE_POLYNOMIAL(curve%surf(2)%ptr%x, '/d/bandrieu/GitHub/FFTsurf/test/Jouke/debug_diffgeominter_surf2.cheb')
    PRINT *,'UV =',UVINTER
    
  end subroutine debug_phf_diffgeom


  subroutine debug_phf_simptinv( &
       curve, &
       uvinter )
    use mod_polynomial
    use mod_types_intersection
    implicit none
    type(type_intersection_curve), intent(in) :: curve
    real(kind=fp),                 intent(in) :: uvinter(2,2)
    
    CALL WRITE_POLYNOMIAL(curve%surf(1)%ptr%x, '/d/bandrieu/GitHub/FFTsurf/test/Jouke/debug_relaxtng_surf1.cheb')
    CALL WRITE_POLYNOMIAL(curve%surf(2)%ptr%x, '/d/bandrieu/GitHub/FFTsurf/test/Jouke/debug_relaxtng_surf2.cheb')
    PRINT *,'UVinter0 =',uvinter
  end subroutine debug_phf_simptinv










  subroutine projection_surface( &
       surf, &
       xyz, &
       uv, &
       lowerb, &
       upperb, &
       stat, &
       xyzproj )
    use mod_diffgeom
    use mod_linalg
    use mod_tolerances
    implicit none
    LOGICAL, PARAMETER :: DEBUG = .false.
    integer, parameter                :: itmax = 20
    type(type_surface), intent(in)    :: surf
    real(kind=fp),      intent(in)    :: xyz(3)
    real(kind=fp),      intent(inout) :: uv(2)
    real(kind=fp),      intent(in)    :: lowerb(2)
    real(kind=fp),      intent(in)    :: upperb(2)
    integer,            intent(out)   :: stat
    real(kind=fp),      intent(out)   :: xyzproj(3)
    real(kind=fp)                     :: r(3), dxyz_duv(3,2), n(3), rn
    real(kind=fp)                     :: mat(2,2), rhs(2), cond, duv(2), erruv
    integer                           :: it, ivar
    
    stat = 1

    IF ( DEBUG ) THEN
       PRINT *,'PROJECTION SURFACE'
       PRINT *,'XYZ =',XYZ
       PRINT *,'UV0 =',UV
    END IF

    IF ( DEBUG ) PRINT *,'|R.T| , |DUV| , eps*cond'
    do it = 1,itmax
       call eval( &
            xyzproj, &
            surf, &
            uv )

       r = xyz - xyzproj
       
       do ivar = 1,2 ! <---------------+
          call evald1( &               !
               dxyz_duv(:,ivar), &     !
               surf, &                 !
               uv, &                   !
               ivar )                  !
       end do ! <----------------------+
       n = cross(dxyz_duv(:,1), dxyz_duv(:,2))
       n = n / norm2(n)
       rn = dot_product(r,n)
       
       mat = matmul(transpose(dxyz_duv), dxyz_duv)
       rhs = matmul(transpose(dxyz_duv), r)

       call linsolve_svd( &
            duv, &
            mat, &
            rhs, &
            2, &
            2, &
            1, &
            cond, &
            tol=EPSmath )

       erruv = sum(duv**2)

       IF ( DEBUG ) PRINT *,sqrt(max(0._fp,sum(r**2) - rn**2)), sqrt(erruv), EPSfp*cond

       ! modify Newton step => keep iterate inside feasible region
       call nd_box_reflexions( &
            uv, &
            lowerb, &
            upperb, &
            duv, &
            2 )

       ! update solution
       uv = uv + duv
       
       ! termination criterion
       if ( erruv < max(EPSuvsqr, EPSfpsqr*cond**2) ) then
          if ( sum(r**2) - rn**2 < EPSxyzsqr ) then
             stat = 0
             call eval( &
                  xyzproj, &
                  surf, &
                  uv )
             IF ( DEBUG ) PRINT *,'CONVERGED, UV =',UV,', XYZ =',XYZPROJ
          else
             IF ( DEBUG ) PRINT *,'STAGNATION'
          end if
          return
       end if
    end do
    
  end subroutine projection_surface
  
end module mod_projection
