module mod_projection

  use mod_math

  implicit none

contains


  

  subroutine projection_hyperface2( &
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
    use mod_linalg
    use mod_diffgeom
    use mod_types_brep
    use mod_halfedge
    use mod_types_intersection
    use mod_intersection
    use mod_tolerances
    use mod_geometry
    use mod_polynomial
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
    real(kind=fp)                             :: duvtmp(2), xyztmp(3), dxyztmp(3)
    real(kind=fp)                             :: xyztarget(3)
    logical                                   :: change_face
    type(type_matrix)                         :: uvpoly(2)
    type(type_intersection_curve), pointer    :: curve => null()
    type(type_intersection_polyline), pointer :: polyline => null()
    integer                                   :: ihead, itail, np, sens
    integer                                   :: ninter
    integer, allocatable                      :: ipls(:,:)
    real(kind=fp), allocatable                :: lamb(:,:)
    integer                                   :: iinter, i1, i2
    real(kind=fp)                             :: w, vec(2)
    logical                                   :: in_open, on_border
    integer                                   :: stat, stat_circ
    real(kind=fp)                             :: uvinter(2,2), xyzinter(3)
    real(kind=fp), dimension(4)               :: lowerb, upperb
    real(kind=fp)                             :: jac(3,2)
    logical                                   :: singular
    integer                                   :: j
    real(kind=fp)                             :: duv_ds(2,2,2), dxyz_ds(3,2)
    real(kind=fp)                             :: ao(2), h, wtmp
    logical                                   :: inside
    integer                                   :: it, iwire, ihedg(2), istart, ivar
    INTEGER :: I

    IF ( DEBUG ) PRINT *,'PROJECTION HYPERFACE'

    upperb(:) = 1._fp + EPSuv
    lowerb = -upperb
    
    ifacetmp = iface
    uvtmp = uv
    stat_proj = 1

    duvtmp = duv
    xyztmp = xyz
    dxyztmp = dxyz
    
    xyztarget = xyz + dxyz

    allocate(uvpoly(1)%mat(2,2))

    IF ( DEBUG ) THEN
       PRINT *,'DXYZ  =',dxyz
       PRINT *,'DUV   =',duv
    END IF

    it = 0
    faces : do
       IF ( DEBUG ) PRINT *,'FACE',ifacetmp
       uvpoly(1)%mat(:,1) = uvtmp
       uvpoly(1)%mat(:,2) = uvtmp + lambmax*duvtmp
       !
       change_face = .false.
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
             polyline => curve%polyline
             ! get current halfedge's uv polyline
             ihead = curve%isplit(2,brep%edges(ihedg(1))%isplit)
             itail = curve%isplit(2,brep%edges(ihedg(1))%isplit + 1)
             np = itail- ihead + 1
             sens = 1 + mod(ihedg(2),2)
             allocate(uvpoly(2)%mat(2,np))
             uvpoly(2)%mat(1:2,1:np) = polyline%uv(1:2,sens,ihead:itail)
             !
             ! compute intersection between the polyline and the uv displacement
             ninter = 0
             call intersect_2Dpolylines( &
                  uvpoly, &
                  [1,1], &
                  [2,np], &
                  ninter, &
                  ipls, &
                  lamb )
             !
             IF ( DEBUG ) THEN
                uvpoly(1)%mat(:,2) = uvtmp + duvtmp
                OPEN(UNIT=13, FILE='polylines/xy1.dat', ACTION='WRITE')
                WRITE (13,*) 2
                DO I = 1,2
                   WRITE (13,*) UVPOLY(1)%MAT(:,I)
                END DO
                CLOSE(13)

                OPEN(UNIT=13, FILE='polylines/xy2.dat', ACTION='WRITE')
                WRITE (13,*) NP
                DO I = 1,NP
                   WRITE (13,*) UVPOLY(2)%MAT(:,I)
                END DO
                CLOSE(13)

                OPEN(UNIT=13, FILE='polylines/result.dat', ACTION='WRITE')
                WRITE (13,*) NINTER
                DO I = 1,NINTER
                   WRITE (13,*) IPLS(:,I)
                END DO
                DO I = 1,NINTER
                   WRITE (13,*) LAMB(:,I)
                END DO
                CLOSE(13)
                PRINT *,'--> CHECK POLYINES'
                PAUSE
             END IF
             !
             if ( ninter > 0 ) then
                ! renormalize lambda
                lamb(1,:) = lamb(1,:) / lambmax
                ! pick point closest to target point
                iinter = maxloc(lamb(1,1:ninter),1)
                IF ( DEBUG ) PRINT *,'LAMBDA =',LAMB(:,IINTER)
                w = lamb(2,iinter)
                i1 = ihead + ipls(2,iinter) - 1
                i2 = i1 + 1
                vec = polyline%uv(1:2,sens,i2) - polyline%uv(1:2,sens,i1)
                !                
                if ( is_in_open_interval(w*norm2(vec), 0._fp, norm2(vec), EPSuv) ) then
                   IF ( DEBUG ) PRINT *,'IN OPEN INTERVAL'
                   in_open = .true.
                   ! on se ramene au sommet de la polyline le plus proche
                   if ( w > 0.5_fp ) then
                      w = 1._fp
                      uvinter = polyline%uv(:,:,i2)
                      xyztmp = polyline%xyz(:,i2)
                   else
                      w = 0._fp
                      uvinter = polyline%uv(:,:,i1)
                      xyztmp = polyline%xyz(:,i1)
                   end if
                   duvtmp = uvtmp + duvtmp - uvinter(:,sens)
                   uvtmp = uvinter(:,sens)
                   do ivar = 1,2 ! <-----------------------+
                      call evald1( &                       !
                           jac(:,ivar), &                  !
                           brep%faces(ifacetmp)%surface, & !
                           uvtmp, &                        !
                           ivar )                          !
                   end do ! <------------------------------+
                   call solve_NxN( &
                        duvtmp, &
                        matmul(transpose(jac), jac), &
                        matmul(transpose(jac), xyztarget - xyztmp), &
                        singular )
                   dxyztmp = matmul(jac, duvtmp)
                else
                   in_open = .false.
                   IF ( DEBUG ) PRINT *,'AT POLYLINE POINT'
                   uvinter = (1._fp - w)*polyline%uv(:,:,i1) + w*polyline%uv(:,:,i2)
                end if
                !
                IF ( DEBUG ) THEN
                   PRINT *,'EXACT INTERSECTION POINT (POLYLINE VERTEX):'
                   PRINT *,'UV  =',UVINTER(:,SENS)
                   PRINT *,'UVinter =',UVINTER
                   PRINT *,'XYZ =',XYZTMP
                   PRINT *,'DUV =',duvtmp
                END IF
                !
                on_border = .false.
                do j = 1,2
                   if ( abs(uvpoly(2)%mat(j,i1-ihead+1)) > 1._fp - EPSuv .and. &
                        abs(uvpoly(2)%mat(j,i2-ihead+1)) > 1._fp - EPSuv .and. &
                        uvpoly(2)%mat(j,i1-ihead+1) * uvpoly(2)%mat(j,i2-ihead+1) > 0._fp ) then
                      on_border = .true.
                      stat_circ = 2 ! infinite radius
                      ao(:) = 0._fp
                      ao(j) = -uvtmp(j)
                      ! ao is the unit, inward-pointing normal to the face boundary
                      IF ( DEBUG ) PRINT *,'ON BORDER'
                      exit 
                   end if
                end do
                !
                if ( .not.on_border ) then
                   !! use circular approximation
                   ! compute the exact tangent direction to the curve
                   call diffgeom_intersection( &
                        curve%surf, &
                        uvinter, &
                        duv_ds, &
                        dxyz_ds, &
                        stat )
                   if ( stat > 1 ) then
                      PRINT *,'STAT DIFFGEOM_INTERSECTION =',STAT
                      CALL WRITE_POLYNOMIAL(curve%surf(1)%ptr%x, 'Jouke/debug_diffgeominter_surf1.cheb')
                      CALL WRITE_POLYNOMIAL(curve%surf(2)%ptr%x, 'Jouke/debug_diffgeominter_surf2.cheb')
                      PRINT *,'UV =',UVINTER
                      STAT_PROJ = 2
                      RETURN
                   end if
                   !
                   ! pick correct orientation
                   duv_ds = sign(1._fp, dot_product(duv_ds(:,1,sens), vec)) * duv_ds
                   if ( sens == 1 ) duv_ds = -duv_ds ! the face is on the left side
                   if ( w > 0.5_fp ) vec = -vec
                   !
                   IF ( DEBUG ) THEN
                      PRINT *,'CIRCULAR APPROXIMATION'
                      PRINT *,'AB = ',VEC
                      PRINT *,'TANGENT =',DUV_DS(:,1,SENS)
                      PRINT *,'DUV =',DUVTMP
                   END IF
                   call circular_approximation( &
                        duv_ds(:,1,sens), &
                        vec, &
                        ao, &
                        stat_circ )
                end if

                if ( stat_circ < 2 ) then
                   ! finite radius
                   if ( sum(duvtmp**2) > sum(ao**2) ) then
                      PRINT *,'/!\ INVALID CIRCULAR APPROXIMATION: |DUV|/R =',norm2(duvtmp)/norm2(ao)
                      PAUSE
                      RETURN
                   end if
                   h = sum((duvtmp - ao)**2) - sum(ao**2)
                   if ( stat_circ == 1 ) h = -h ! face boundary locally concave
                   inside = ( h < -EPSmath )
                   !
                   IF ( DEBUG) THEN
                      PRINT *,'STAT_CIRC =',STAT_CIRC
                      PRINT *,'AO =',AO,', |RAD| =',NORM2(AO)
                      PRINT *,'|UV + DUV - O|/|RAD| =',SQRT(sum((duvtmp - ao)**2) / sum(ao**2))
                      PRINT *,'H/L =',h / NORM2(DUVTMP)
                      PRINT *,'L/R =',SQRT(SUM(DUVTMP**2) / SUM(AO**2))
                      PRINT *,'H/R =',h / NORM2(AO)
                   END IF
                   !
                   if ( inside .and. in_open ) then
                      ! the polyline intersection is a false-positive
                      ! refine locally the polyline by adding a new point on the current segment
                      wtmp = lamb(2,iinter)
                      IF ( DEBUG ) PRINT *,'WTMP =',WTMP
                      uvinter = (1._fp - wtmp)*polyline%uv(:,:,i1) + wtmp*polyline%uv(:,:,i2)
                      ! relax to exact intersectino curve
                      call simultaneous_point_inversions( &
                           curve%surf, &
                           lowerb, &
                           upperb, &
                           stat, &
                           uvinter, &
                           xyzinter )
                      if ( stat > 0 ) then
                         ! failed to converge :(
                         PRINT *,'FAILED TO RELAX ON TANGENTIAL INTERSECTION'
                         CALL WRITE_POLYNOMIAL(curve%surf(1)%ptr%x, 'Jouke/debug_relaxtng_surf1.cheb')
                         CALL WRITE_POLYNOMIAL(curve%surf(2)%ptr%x, 'Jouke/debug_relaxtng_surf2.cheb')
                         PRINT *,'UVinter0 =',(1._fp - wtmp)*polyline%uv(:,:,i1) + wtmp*polyline%uv(:,:,i2)
                         RETURN 
                      else
                         ! add polyline point
                         call insert_polyline_point( &
                              uvinter, &
                              xyzinter, &
                              stat, &
                              polyline, &
                              i1 )
                         if ( stat == 0 ) then
                            IF ( DEBUG ) PRINT *,'NEW POINT INSERTED IN POLYLINE : UV =',uvinter
                            ! shift downstream polyline split points
                            where ( curve%isplit(2,:) >= i2 ) curve%isplit(2,:) = curve%isplit(2,:) + 1
                         end if
                      end if
                   end if
                   !
                   if ( .not.inside .and. &
                        h < 5.d-2 * norm2(duvtmp) ) then
                      ! the target point is just outside (the displacement is nearly
                      ! tangent to the circular approximation)
                      wtmp = dot_product(duvtmp, vec) / sum(vec**2)
                      if ( w > 0.5_fp ) wtmp = 1._fp - wtmp
                      IF ( DEBUG ) PRINT *,'NEARLY TANGENTIAL DISPLACEMENT'
                      IF ( DEBUG ) PRINT *,'WTMP =',WTMP
                      uvinter = (1._fp - wtmp)*polyline%uv(:,:,i1) + wtmp*polyline%uv(:,:,i2)
                      ! relax to exact intersection curve
                      call simultaneous_point_inversions( &
                           curve%surf, &
                           lowerb, &
                           upperb, &
                           stat, &
                           uvinter, &
                           xyzinter )
                      if ( stat > 0 ) then
                         ! failed to converge :(
                         PRINT *,'FAILED TO RELAX ON TANGENTIAL INTERSECTION'
                         CALL WRITE_POLYNOMIAL(curve%surf(1)%ptr%x, 'Jouke/debug_relaxtng_surf1.cheb')
                         CALL WRITE_POLYNOMIAL(curve%surf(2)%ptr%x, 'Jouke/debug_relaxtng_surf2.cheb')
                         PRINT *,'UVinter0 =',(1._fp - wtmp)*polyline%uv(:,:,i1) + wtmp*polyline%uv(:,:,i2)
                         RETURN 
                      else
                         ! add polyline point
                         call insert_polyline_point( &
                              uvinter, &
                              xyzinter, &
                              stat, &
                              polyline, &
                              i1 )
                         if ( stat == 0 ) then
                            IF ( DEBUG ) PRINT *,'NEW POINT INSERTED IN POLYLINE : UV =',uvinter
                            ! shift downstream polyline split points
                            where ( curve%isplit(2,:) >= i2 ) curve%isplit(2,:) = curve%isplit(2,:) + 1
                         end if
                         !
                         duvtmp = uvinter(1:2,sens) - uvtmp
                         IF ( DEBUG ) PRINT *,'MODIF DUVTMP =',duvtmp
                         inside = .true.
                      end if
                   end if
                else
                   ! infinite radius
                   inside = ( dot_product(ao,duvtmp) > -EPSmath )
                end if
                IF ( DEBUG ) PRINT *,'INSIDE?',INSIDE
                change_face = .not.inside
                !
                if ( change_face ) then
                   ! change brep face
                   ifacetmp = get_face(brep, get_twin(ihedg))
                   ! uv coordinates in new face
                   uvtmp = uvinter(1:2,ihedg(2))
                   ! new point displacement
                   do ivar = 1,2 ! <-----------------------+
                      call evald1( &                       !
                           jac(:,ivar), &                  !
                           brep%faces(ifacetmp)%surface, & !
                           uvtmp, &                        !
                           ivar )                          !
                   end do ! <------------------------------+
                   call solve_NxN( &
                        duvtmp, &
                        matmul(transpose(jac), jac), &
                        matmul(transpose(jac), xyztarget - xyztmp), &
                        singular )
                   dxyztmp = matmul(jac, duvtmp)
                   IF ( DEBUG ) THEN
                      PRINT *,'CHANGE BREP FACE -->',IFACETMP
                      PRINT *,'  UVTMP =',uvtmp
                      PRINT *,' DUVTMP =',duvtmp
                      PRINT *,' XYZTMP =',xyztmp
                      PRINT *,'DXYZTMP =',dxyztmp
                   END IF
                end if
                IF ( DEBUG ) PAUSE
             end if
             !
             deallocate(uvpoly(2)%mat)
             if ( change_face ) exit wires
             !
             ! move on to the next halfedge on the wire
             ihedg = get_next(brep, ihedg)
             if ( get_orig(brep, ihedg) == istart ) exit halfedges
          end do halfedges
          !
       end do wires
       !
       if ( .not.change_face ) then
          IF ( DEBUG ) PRINT *,'CONVERGED, IT =',IT
          stat_proj = 0
          exit faces
       end if
       !
       it = it + 1
       IF ( IT > BREP%NF ) THEN
          PRINT *,'*** FAILED TO PROJECT ON HYPERFACE ***'
          DUVTMP(:) = 0._FP
          EXIT faces
       END IF
    end do faces

    uvtmp = uvtmp + duvtmp
    deallocate(uvpoly(1)%mat)
    
  end subroutine projection_hyperface2


  
  subroutine circular_approximation( &
       t, &
       ab, &
       ao, &
       stat )
    implicit none
    ! Given two points a, b and a vector t,
    ! returns the vector from a to the center of the (unique) circle
    ! passing through a and b, that is tangent to tng at a.
    ! If t and ab are colinear, the circle has inifite radius, therefore
    ! return the unit direction from a to the center (at infinite)
    ! t  : tangent to the circle at point a
    ! ab : vector from point a to point b (both on the circle)
    ! ao : vector from point a to the center o of the circle
    ! stat = 0 : positive radius (convex)
    !        1 : negative radius (concave)
    !        2 : infinite radius (flat)
    real(kind=fp), intent(in)  :: t(2)
    real(kind=fp), intent(in)  :: ab(2)
    real(kind=fp), intent(out) :: ao(2)
    integer,       intent(out) :: stat
    real(kind=fp)              :: denom, rad

    ao = [-t(2), t(1)] / norm2(t)
    denom = dot_product(ab, ao)

    if ( abs(denom) < EPSmath ) then
       stat = 2
    else
       rad = 0.5_fp * sum(ab**2) / denom
       if ( rad > 0._fp ) then
          stat = 0
       else
          stat = 1
       end if
       ao = rad*ao
    end if
    
  end subroutine circular_approximation




  subroutine is_inside_circular_approximation( &
       duv, &
       tng, &
       vec, &
       inside )
    implicit none
    real(kind=fp), intent(in)  :: duv(2)
    real(kind=fp), intent(in)  :: tng(2)
    real(kind=fp), intent(in)  :: vec(2)
    logical,       intent(out) :: inside
    real(kind=fp)              :: nor(2), denom, rad, dist

    nor = [-tng(2), tng(1)] / norm2(tng)
    denom = dot_product(vec, nor)

    if ( abs(denom) < EPSmath ) then
       inside = ( dot_product(nor,duv) > -EPSmath )
    else
       rad = 0.5_fp * sum(vec**2) / denom
       dist = sum((duv - rad*nor)**2)
       inside = ( dist < rad**2 )
       !PRINT *,'RAD =',RAD
       !PRINT *,'CTR - UV =',RAD*NOR
       PRINT *,'(DIST-R)/R =',(SQRT(DIST) - ABS(RAD))/ABS(RAD), '&
            &, |DUV|/R =',NORM2(DUV)/ABS(RAD)
       if ( rad < 0._fp ) inside = .not.inside
    end if
    
  end subroutine is_inside_circular_approximation


  

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
    use mod_linalg
    use mod_diffgeom
    use mod_types_brep
    use mod_halfedge
    use mod_types_intersection
    use mod_intersection
    use mod_tolerances
    use mod_geometry
    use mod_polynomial
    implicit none
    logical, intent(in) :: debug
    real(kind=fp), parameter                  :: lambdamax = 1._fp
    type(type_brep), intent(in)               :: brep
    integer,         intent(in)               :: iface
    real(kind=fp),   intent(in)               :: uv(2)
    real(kind=fp),   intent(in)               :: xyz(3)
    real(kind=fp),   intent(in)               :: duv(2)
    real(kind=fp),   intent(in)               :: dxyz(3)
    integer,         intent(out)              :: ifacetmp
    real(kind=fp),   intent(out)              :: uvtmp(2)
    integer,         intent(out)              :: stat_proj
    real(kind=fp)                             :: xyztarget(3)
    real(kind=fp)                             :: duvtmp(2), normduv
    real(kind=fp)                             :: xyztmp(3)
    real(kind=fp)                             :: dxyztmp(3)
    logical                                   :: change_face
    type(type_intersection_curve), pointer    :: curve => null()
    type(type_intersection_polyline), pointer :: polyline => null()
    integer                                   :: i1, i2, np, sens
    integer                                   :: ninter
    type(type_matrix)                         :: uvpoly(2)
    integer, allocatable                      :: ipls(:,:)
    real(kind=fp), allocatable                :: lambda(:,:)
    integer                                   :: iinter
    real(kind=fp)                             :: uvinter(2,2), w, sgn, det, detpoly, jac(3,2), xyzinter(3)
    real(kind=fp)                             :: tng(2), duv_ds(2,2,2), dxyz_ds(3,2)
    logical                                   :: singular
    integer                                   :: it, iwire, ihedg(2), istart, ivar, iedgeprev
    real(kind=fp), dimension(4)               :: lowerb, upperb
    integer                                   :: stat, imax, imin
    INTEGER :: I

    IF ( DEBUG ) PRINT *,'PROJECTION HYPERFACE'

    stat_proj = 1
    
    upperb(:) = 1._fp + EPSuv
    lowerb = -upperb
    
    xyztarget = xyz + dxyz

    ifacetmp = iface
    uvtmp = uv
    duvtmp = duv
    xyztmp = xyz
    dxyztmp = dxyz

    normduv = norm2(duvtmp)

    IF ( DEBUG ) THEN
       PRINT *,'DXYZ  =',dxyz
       PRINT *,'DUV   =',duv
    END IF

    allocate(uvpoly(1)%mat(2,2))


    iedgeprev = 0
    it = 0
    faces : do ! <-----------------------------------------------------------------------------------+
       IF ( DEBUG ) PRINT *,'FACE',ifacetmp
       change_face = .false.                                                                         !
       uvpoly(1)%mat(:,1) = uvtmp                                                                    !
       uvpoly(1)%mat(:,2) = uvtmp + lambdamax * duvtmp                                               !
       !                                                                                             !
       wires : do iwire = 0,brep%faces(ifacetmp)%ninner ! <---------------------------------------+  !
          ! get index of first halfedge on the wire                                               !  !
          if ( iwire == 0 ) then ! <---------------------+                                        !  !
             ihedg = brep%faces(ifacetmp)%outer          !                                        !  !
          else ! ----------------------------------------+                                        !  !
             ihedg = brep%faces(ifacetmp)%inner(:,iwire) !                                        !  !
          end if ! <-------------------------------------+                                        !  !
          ! start brep vertex                                                                     !  !
          istart = get_orig(brep, ihedg)                                                          !  !
          ! cycle the wire's halfedges                                                            !  !
          halfedges : do ! <-------------------------------------------------------------------+  !  !
             IF ( DEBUG ) PRINT *,'HALFEDGE',ihedg
             if ( ihedg(1) == iedgeprev ) GOTO 100                                             !  !  !
             curve => brep%edges(ihedg(1))%curve                                               !  !  !
             ! get current halfedge's uv polyline                                              !  !  !
             i1 = curve%isplit(2,brep%edges(ihedg(1))%isplit)                                  !  !  !
             i2 = curve%isplit(2,brep%edges(ihedg(1))%isplit + 1)                              !  !  !
             np = i2 - i1 + 1                                                                  !  !  !
             sens = 1 + mod(ihedg(2),2)                                                        !  !  !
             polyline => curve%polyline                                                        !  !  !
             allocate(uvpoly(2)%mat(2,np))                                                     !  !  !
             uvpoly(2)%mat(1:2,1:np) = polyline%uv(1:2,sens,i1:i2)                             !  !  !
             !                                                                                 !  !  !
             ! compute intersection between the polyline and the uv displacement               !  !  !
             ninter = 0                                                                        !  !  !
             call intersect_2Dpolylines( &                                                     !  !  !
                  uvpoly, &                                                                    !  !  !
                  [1,1], &                                                                     !  !  !
                  [2,np], &                                                                    !  !  !
                  ninter, &                                                                    !  !  !
                  ipls, &                                                                      !  !  !
                  lambda )                                                                     !  !  !
             IF ( DEBUG ) THEN
                OPEN(UNIT=13, FILE='polylines/xy1.dat', ACTION='WRITE')
                WRITE (13,*) 2
                DO I = 1,2
                   WRITE (13,*) UVPOLY(1)%MAT(:,I)
                END DO
                CLOSE(13)

                OPEN(UNIT=13, FILE='polylines/xy2.dat', ACTION='WRITE')
                WRITE (13,*) NP
                DO I = 1,NP
                   WRITE (13,*) UVPOLY(2)%MAT(:,I)
                END DO
                CLOSE(13)

                OPEN(UNIT=13, FILE='polylines/result.dat', ACTION='WRITE')
                WRITE (13,*) NINTER
                DO I = 1,NINTER
                   WRITE (13,*) IPLS(:,I)
                END DO
                DO I = 1,NINTER
                   WRITE (13,*) LAMBDA(:,I)
                END DO
                CLOSE(13)
                PRINT *,'--> CHECK POLYINES'
                PAUSE
             END IF
             !                                                                                 !  !  !
             if ( ninter > 0 ) then ! <-----------------------------------------------------+  !  !  !
                ! pick closest point to target uv point                                     !  !  !  !
                lambda(1,:) = lambda(1,:) / lambdamax
                iinter = maxloc(lambda(1,1:ninter), 1)                                      !  !  !  !
                IF ( DEBUG ) PRINT *,'LAMBDA =',lambda(:,iinter)
                i1 = i1 + ipls(2,iinter) - 1                                                !  !  !  !
                i2 = i1 + 1                                                                 !  !  !  !
                ! relax onto exact intersection curve                                       !  !  !  !
                w = lambda(2,iinter)                                                        !  !  !  !
                uvinter = (1._fp - w)*polyline%uv(1:2,1:2,i1) + w*polyline%uv(1:2,1:2,i2)   !  !  !  !
                call simultaneous_point_inversions( &                                       !  !  !  !
                     curve%surf, &                                                          !  !  !  !
                     lowerb, &                                                              !  !  !  !
                     upperb, &                                                              !  !  !  !
                     stat, &                                                                !  !  !  !
                     uvinter, &                                                             !  !  !  !
                     xyztmp )                                                               !  !  !  !
                if ( stat > 0 ) then ! <------------------------------------------------+   !  !  !  !
                   ! failed to converge :(                                              !   !  !  !  !
                   PRINT *,'FAILED TO RELAX ON TANGENTIAL INTERSECTION'
                   CALL WRITE_POLYNOMIAL(curve%surf(1)%ptr%x, 'Jouke/debug_relaxtng_surf1.cheb')
                   CALL WRITE_POLYNOMIAL(curve%surf(2)%ptr%x, 'Jouke/debug_relaxtng_surf2.cheb')
                   PRINT *,'UVTNG =',uvinter
                   PAUSE 
                end if ! <--------------------------------------------------------------+   !  !  !  !
                IF ( DEBUG ) THEN
                   PRINT *,'EXACT INTERSECTION POINT :'
                   PRINT *,'UV  =',UVINTER(:,SENS)
                   PRINT *,'XYZ =',XYZTMP
                END IF
                ! compute exact lambda                                                      !  !  !  !
                IF ( DEBUG ) PRINT *,'|UVTMP - UVEXACT| =',NORM2(uvinter(:,sens) - uvtmp)
                IF ( SUM((uvinter(:,sens) - uvtmp)**2) < EPSUVSQR ) THEN
                   IF ( DEBUG ) PRINT *,'CORRECTION W:',dot_product(uvinter(:,sens) - uvtmp, duvtmp) / sum(duvtmp**2),'&
                        & -->',LAMBDA(1,IINTER)
                   w = lambda(1,iinter)
                ELSE
                   w = dot_product(uvinter(:,sens) - uvtmp, duvtmp) / sum(duvtmp**2)
                END IF
                IF ( DEBUG ) PRINT *,'LAMBDA* =',W,', *|DUV| =',W*normduv
                if ( abs(w)*normduv < EPSuv ) then ! <----------------------------------+   !  !  !  !
                   ! the current point is exactly on the current face's boundary        !   !  !  !  !
                   ! check if it is entering or exiting that face                       !   !  !  !  !
                   imax = maxloc(abs(uvinter(:,sens)),1)                                !   !  !  !  !
                   imin = 1 + mod(imax,2)                                               !   !  !  !  !
                   if ( abs(uvinter(imax,sens)) > 1._fp - EPSuv ) then ! <-----------+  !   !  !  !  !
                      if ( abs(uvinter(imin,sens)) > 1._fp - EPSuv ) then ! <-----+  !  !   !  !  !  !
                         ! corner                                                 !  !  !   !  !  !  !
                         IF ( DEBUG ) PRINT *,'CORNER',duvtmp*uvinter(:,sens)
                         change_face = ( any(duvtmp*uvinter(:,sens) > EPSmath) ) 
                      else ! -----------------------------------------------------+  !  !   !  !  !  !
                         ! point along a border                                   !  !  !   !  !  !  !
                         IF ( DEBUG ) PRINT *,'BORDER',duvtmp(imax)*uvinter(imax,sens)
                         change_face = ( duvtmp(imax)*uvinter(imax,sens) > EPSmath )
                      end if ! <--------------------------------------------------+  !  !   !  !  !  !
                   else ! <----------------------------------------------------------+  !   !  !  !  !
                      call diffgeom_intersection( &                                     !   !  !  !  !
                           curve%surf, &                                                !   !  !  !  !
                           uvinter, &                                                   !   !  !  !  !
                           duv_ds, &                                                    !   !  !  !  !
                           dxyz_ds, &                                                   !   !  !  !  !
                           stat )                                                       !   !  !  !  !
                      if ( stat > 1 ) then
                         ! take action...
                         PRINT *,'STAT DIFFGEOM_INTERSECTION =',STAT
                         CALL WRITE_POLYNOMIAL(curve%surf(1)%ptr%x, 'Jouke/debug_diffgeominter_surf1.cheb')
                         CALL WRITE_POLYNOMIAL(curve%surf(2)%ptr%x, 'Jouke/debug_diffgeominter_surf2.cheb')
                         PRINT *,'UV =',UVINTER
                         STOP
                      end if
                      tng = real((-1)**sens, kind=fp) * (polyline%uv(:,sens,i2) - polyline%uv(:,sens,i1))
                      sgn = sign(1._fp, dot_product(duv_ds(:,1,sens), tng))             !   !  !  !  !
                      duv_ds = sgn * duv_ds                                             !   !  !  !  !
                      det = duv_ds(1,1,sens)*duvtmp(2) - duv_ds(2,1,sens)*duvtmp(1)     !   !  !  !  !
                      IF ( DEBUG ) THEN
                         PRINT *,'TANGENT =',DUV_DS(:,1,SENS)
                         PRINT *,'DUV =',DUVTMP
                         PRINT *,'DET =',DET
                      END IF
                      detpoly = tng(1)*duvtmp(2) - tng(2)*duvtmp(1)                     !   !  !  !  !
                      do
                         exit
                         IF ( DEBUG ) PRINT *,'DETPOLY =',DETPOLY
                         if ( detpoly*det > 0._fp ) exit
                         ! add new point to polyline
                         uvinter = 0.5_fp * ( &
                              polyline%uv(1:2,1:2,i1) + &
                              polyline%uv(1:2,1:2,i2) )
                         call simultaneous_point_inversions( &
                              curve%surf, &
                              lowerb, &
                              upperb, &
                              stat, &
                              uvinter, &
                              xyzinter )
                         if ( stat > 0 ) exit
                         call insert_polyline_point( &
                              uvinter, &
                              xyzinter, &
                              stat, &
                              polyline, &
                              i1 )
                         if ( stat > 0 ) exit
                         where ( curve%isplit(2,:) >= i2 ) curve%isplit(2,:) = curve%isplit(2,:) + 1
                         IF ( DEBUG ) PRINT *,'NEW POINT INSERTED IN POLYLINE : UV =',uvinter
                         tng = real((-1)**sens, kind=fp) * (polyline%uv(:,sens,i2) - polyline%uv(:,sens,i1))
                         detpoly = tng(1)*duvtmp(2) - tng(2)*duvtmp(1)
                      end do

                      IF ( DEBUG ) THEN
                         i1 = curve%isplit(2,brep%edges(ihedg(1))%isplit)
                         i2 = curve%isplit(2,brep%edges(ihedg(1))%isplit + 1)
                         np = i2 - i1 + 1       
                         OPEN(UNIT=13, FILE='polylines/xy3.dat', ACTION='WRITE')
                         WRITE (13,*) np
                         DO I = i1,i2
                            WRITE (13,*) polyline%uv(:,sens,i)
                         END DO
                         CLOSE(13)
                      END IF
  
                      IF ( DEBUG ) THEN
                         !PRINT *,'TANGENT =',DUV_DS(:,1,SENS)
                         !PRINT *,'DUV =',DUVTMP
                         !PRINT *,'DET =',DET
                         !PRINT *,'DETPOLY =',DETPOLY
                         PAUSE
                      END IF
                      change_face = ( det < -EPSmath )                                  !   !  !  !  !
                   end if ! <-----------------------------------------------------------+   !  !  !  !
                   !                                                                        !  !  !  !
                elseif ( is_in_closed_interval(w, 0._fp, 1._fp, EPSmath) ) then ! ------+   !  !  !  !
                   ! the target point is in another brep face                           !   !  !  !  !
                   change_face = .true.                                                 !   !  !  !  !
                end if ! <--------------------------------------------------------------+   !  !  !  !
                !                                                                           !  !  !  !
                if ( change_face ) then ! <---------------------------------------------+   !  !  !  !
                   iedgeprev = ihedg(1)                                                 !   !  !  !  !
                   ifacetmp = get_face(brep, get_twin(ihedg))                           !   !  !  !  !
                   uvtmp = uvinter(1:2,ihedg(2))                                        !   !  !  !  !
                   ! update uv displacement                                             !   !  !  !  !
                   do ivar = 1,2 ! <-----------------------+                            !   !  !  !  !
                      call evald1( &                       !                            !   !  !  !  !
                           jac(:,ivar), &                  !                            !   !  !  !  !
                           brep%faces(ifacetmp)%surface, & !                            !   !  !  !  !
                           uvtmp, &                        !                            !   !  !  !  !
                           ivar )                          !                            !   !  !  !  !
                   end do ! <------------------------------+                            !   !  !  !  !
                   call solve_NxN( &                                                    !   !  !  !  !
                        duvtmp, &                                                       !   !  !  !  !
                        matmul(transpose(jac), jac), &                                  !   !  !  !  !
                        matmul(transpose(jac), xyztarget - xyztmp), &                   !   !  !  !  !
                        singular )                                                      !   !  !  !  !
                   normduv = norm2(duvtmp)
                   ! update xyz displacement                                            !   !  !  !  !
                   dxyztmp = matmul(jac, duvtmp)                                        !   !  !  !  !
                   IF ( DEBUG ) THEN
                      PRINT *,'CHANGE BREP FACE -->',IFACETMP
                      PRINT *,'  UVTMP =',uvtmp
                      PRINT *,' DUVTMP =',duvtmp
                      PRINT *,' XYZTMP =',xyztmp
                      PRINT *,'DXYZTMP =',dxyztmp
                      PAUSE
                   END IF
                end if ! <--------------------------------------------------------------+   !  !  !  !
                !                                                                           !  !  !  !
             end if ! <---------------------------------------------------------------------+  !  !  !
             deallocate(uvpoly(2)%mat)                                                         !  !  !
             if ( change_face ) exit wires                                                     !  !  !
             ! move on to the next halfedge on the wire                                        !  !  !
             100 ihedg = get_next(brep, ihedg)                                                 !  !  !
             if ( get_orig(brep, ihedg) == istart ) exit halfedges                             !  !  !
             !                                                                                 !  !  !
          end do halfedges ! <-----------------------------------------------------------------+  !  !
          !                                                                                       !  !
       end do wires ! <---------------------------------------------------------------------------+  !
       !
       if ( .not.change_face ) then
          IF ( DEBUG ) PRINT *,'CONVERGED, IT =',IT
          stat_proj = 0
          exit faces
       end if
       it = it + 1
       IF ( IT > BREP%NF ) THEN
          PRINT *,'*** FAILED TO PROJECT ON HYPERFACE ***'
          DUVTMP(:) = 0._FP
          EXIT faces
       END IF
    end do faces ! <----------------------------------------------------------------------------+

    uvtmp = uvtmp + duvtmp
    deallocate(uvpoly(1)%mat)




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
       debug )
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
    real(kind=fp)                             :: duvtmp(2,2)
    real(kind=fp), dimension(3)               :: xyzprev, xyztmp, xyztarget, dxyztmp, p
    type(type_intersection_curve), pointer    :: curve => null()
    type(type_intersection_polyline), pointer :: polyline => null()
    real(kind=fp)                             :: w, wfirst, wlast, ds, duv_ds(2,2,2), dxyz_ds(3,2)
    integer                                   :: jedge, ihedg(2), ifirst, ilast, sens, np, k, stat
    real(kind=fp), dimension(4)               :: lowerb, upperb

    IF ( DEBUG ) PRINT *,'PROJECTION HYPEREDGE'
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

    IF ( DEBUG ) THEN
       PRINT *,'IHEDGE =', IHEDG
       PRINT *,' FACES =', get_face(brep, [iedgetmp,1]), get_face(brep, [iedgetmp,2])
       PRINT *,'    UV =', UV
       PRINT *,'   XYZ =', XYZ
       PRINT *,'  DXYZ =', DXYZ
    END IF
    
    do k = 1,hyperedge%ne
       xyztarget = xyztmp + dxyztmp
       
       p = real((-1)**sens, kind=fp) * curve%param_vector
       wfirst = dot_product(p, polyline%xyz(:,ifirst))
       wlast  = dot_product(p, polyline%xyz(:,ilast))
       w = dot_product(p, xyztarget)
       
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
             PAUSE
          else
             IF ( DEBUG ) THEN
                PRINT *,'CONVERGED :'
                PRINT *,'IHEDGE =', IHEDG
                PRINT *,' FACES =', get_face(brep, [iedgetmp,1]), get_face(brep, [iedgetmp,2])
                PRINT *,'    UV =', UVTMP
                PRINT *,'   XYZ =', XYZTMP
             END IF
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
    STOP 'projection_hyperedge : failed to reproject'
    
  end subroutine projection_hyperedge








end module mod_projection
