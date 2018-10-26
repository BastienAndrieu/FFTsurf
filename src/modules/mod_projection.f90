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
    real(kind=fp)                             :: w, vec(2), normvec
    logical                                   :: in_open, on_border
    integer                                   :: stat, stat_circ
    real(kind=fp)                             :: uvinter(2,2), xyzinter(3)
    real(kind=fp), dimension(4)               :: lowerb, upperb
    real(kind=fp)                             :: jac(3,2)
    logical                                   :: singular
    integer                                   :: j
    real(kind=fp)                             :: duv_ds(2,2,2), dxyz_ds(3,2)
    real(kind=fp)                             :: ao(2), h, wtmp, wq, aux(2)
    logical                                   :: inside
    real(kind=fp), dimension(2)               :: p, q, pu, pq
    integer                                   :: it, iwire, ihedg(2), istart, ivar

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
       ! the provisional point uvtmp is inside the parametric domain of
       ! the face ifacetmp (possibly on its boundary)
       ! check whether the point uvtmp + duvtmp is inside/outside this domain
       uvpoly(1)%mat(:,1) = uvtmp
       uvpoly(1)%mat(:,2) = uvtmp + lambmax*duvtmp
       !
       change_face = .false.
       ! to do so, intersect the segment [uvtmp, uvtmp+duvtmp] with all the curve segments
       ! that form the boundary of the domain (outer an inner wires of the face)
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
                call debug_phf_polylines( &
                     uvpoly, &
                     ninter, &
                     ipls, &
                     lamb )
                PRINT *,'--> CHECK POLYLINES'
                PAUSE
             END IF
             !
             if ( ninter > 0 ) then
                ! renormalize lambda
                lamb(1,:) = lamb(1,:) / lambmax
                ! pick point closest to target point
                iinter = maxloc(lamb(1,1:ninter),1)
                w = lamb(2,iinter)
                i1 = ihead + ipls(2,iinter) - 1
                i2 = i1 + 1
                IF ( DEBUG ) THEN
                   PRINT *,'LAMBDA =',LAMB(:,IINTER)
                   PRINT *,'POLYLINE SEGMENT:'
                   PRINT *,'XYZ ='
                   PRINT *,polyline%XYZ(1:3,i1)
                   PRINT *,polyline%XYZ(1:3,i2)
                   PRINT *,'UV ='
                   PRINT *,polyline%UV(1:2,sens,i1)
                   PRINT *,polyline%UV(1:2,sens,i2)
                END IF
                vec = polyline%uv(1:2,sens,i2) - polyline%uv(1:2,sens,i1)
                normvec = norm2(vec)
                !
                ! check if the intersection point is in the interior of a polyline segment
                ! or at one of its vertices
                if ( is_in_open_interval(w*normvec, 0._fp, normvec, EPSuv) ) then ! <---+
                   ! the intersection point in in the interior of a segment             !
                   IF ( DEBUG ) PRINT *,'IN OPEN INTERVAL'                              !
                   in_open = .true.                                                     !
                   ! snap uvtmp to the closest polyline vertex (i.e. to the face's      !
                   ! boundary)                                                          ! 
                   if ( w > 0.5_fp ) then ! <--------------+                            !
                      w = 1._fp                            !                            !
                      uvinter = polyline%uv(:,:,i2)        !                            !
                      xyztmp = polyline%xyz(:,i2)          !                            !
                   else ! ---------------------------------+                            !
                      w = 0._fp                            !                            !
                      uvinter = polyline%uv(:,:,i1)        !                            !
                      xyztmp = polyline%xyz(:,i1)          !                            !
                   end if ! <------------------------------+                            !
                   ! update uv displacement                                             !
                   duvtmp = uvtmp + duvtmp - uvinter(:,sens)                            !
                   ! update uv coordinates                                              !
                   uvtmp = uvinter(:,sens)                                              !
                   ! update xyz displacement                                            !
                   do ivar = 1,2 ! <-----------------------+                            !
                      call evald1( &                       !                            !
                           jac(:,ivar), &                  !                            !
                           brep%faces(ifacetmp)%surface, & !                            !
                           uvtmp, &                        !                            !
                           ivar )                          !                            !
                   end do ! <------------------------------+                            !
                   call solve_NxN( &                                                    !
                        duvtmp, &                                                       !
                        matmul(transpose(jac), jac), &                                  !
                        matmul(transpose(jac), xyztarget - xyztmp), &                   !
                        singular )                                                      !
                   dxyztmp = matmul(jac, duvtmp)                                        !
                else ! -----------------------------------------------------------------+
                   ! the intersection point is a polyline vertex, so it lies exactly on !
                   ! the face's boundary                                                !
                   in_open = .false.                                                    !
                   IF ( DEBUG ) PRINT *,'AT POLYLINE POINT'                             !
                   uvinter = (1._fp - w)*polyline%uv(:,:,i1) + w*polyline%uv(:,:,i2)    !
                   duvtmp = uvtmp + duvtmp - uvinter(:,sens)                            !
                   uvtmp = uvinter(:,sens)                                              !
                end if ! <--------------------------------------------------------------+
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
                do j = 1,2 ! <---------------------------------------------------------+
                   aux = uvpoly(2)%mat(j,[i1,i2]-ihead+1)                              !
                   if ( all(abs(aux) > 1._fp - EPSuv) .and. &                          !
                        product(aux) > 0._fp ) then ! <-----------------------------+  !
                      on_border = .true.                                            !  !
                      stat_circ = 2 ! infinite radius                               !  !
                      ao(:) = 0._fp                                                 !  !
                      ao(j) = -uvtmp(j)                                             !  !
                      ! ao is the unit, inward-pointing normal to the face boundary !  !
                      IF ( DEBUG ) THEN                                             !  !
                         PRINT *,'ON BORDER'                                        !  !
                         PRINT *,'INWARD NORMAL =',AO                               !  !
                      END IF                                                        !  !
                      exit                                                          !  !
                   end if ! <-------------------------------------------------------+  !
                end do ! <-------------------------------------------------------------+
                !
                if ( .not.on_border ) then ! <-----------------------------------------+
                   !! use circular approximation                                       !
                   ! compute the exact tangent direction to the curve                  !
                   call diffgeom_intersection( &                                       !
                        curve%surf, &                                                  !
                        uvinter, &                                                     !
                        duv_ds, &                                                      !
                        dxyz_ds, &                                                     !
                        stat )                                                         !
                   if ( stat > 1 ) then
                      PRINT *,'STAT DIFFGEOM_INTERSECTION =',STAT
                      call debug_phf_diffgeom( &
                           curve, &
                           uvinter )
                      STAT_PROJ = 2
                      RETURN
                   end if
                   !                                                                   !
                   ! pick correct orientation                                          !
                   duv_ds = sign(1._fp, dot_product(duv_ds(:,1,sens), vec)) * duv_ds   !
                   if ( sens == 1 ) duv_ds = -duv_ds ! the face is on the left side    !
                   if ( w > 0.5_fp ) vec = -vec                                        !
                   !                                                                   !
                   IF ( DEBUG ) THEN
                      PRINT *,'CIRCULAR APPROXIMATION'
                      PRINT *,'AB = ',VEC
                      PRINT *,'TANGENT =',DUV_DS(:,1,SENS)
                      PRINT *,'DUV =',DUVTMP
                   END IF
                   call circular_approximation( &                                      !
                        duv_ds(:,1,sens), &                                            !
                        vec, &                                                         !
                        ao, &                                                          !
                        stat_circ )                                                    !
                end if ! <-------------------------------------------------------------+
                !
                if ( stat_circ < 2 ) then ! <------------------------------------------+
                   ! finite radius                                                     !
                   if ( sum(duvtmp**2) > sum(ao**2) ) then
                      PRINT *,'/!\ INVALID CIRCULAR APPROXIMATION: |DUV|/R =',norm2(duvtmp)/norm2(ao)
                      PAUSE
                      RETURN
                   end if
                   h = sum((duvtmp - ao)**2) - sum(ao**2)                              !
                   if ( stat_circ == 1 ) h = -h ! face boundary locally concave        !
                   inside = ( h < -EPSmath )                                           !
                   !                                                                   !
                   IF ( DEBUG) THEN
                      PRINT *,'STAT_CIRC =',STAT_CIRC
                      PRINT *,'AO =',AO,', |RAD| =',NORM2(AO)
                      PRINT *,'|UV + DUV - O|/|RAD| =',SQRT(sum((duvtmp - ao)**2) / sum(ao**2))
                      PRINT *,'H/L =',h / NORM2(DUVTMP)
                      PRINT *,'L/R =',SQRT(SUM(DUVTMP**2) / SUM(AO**2))
                      PRINT *,'H/R =',h / NORM2(AO)
                   END IF
                   !                                                                   !
                   if ( inside ) then ! <-------------------------------------------+
                      do ! <------------------------------------------------------+
                         aux = uvtmp + duvtmp - polyline%uv(:,sens,i1)            !
                         vec = polyline%uv(1:2,sens,i2) - polyline%uv(1:2,sens,i1)!
                         if ( sens == 1 ) vec = -vec                              !
                         IF ( DEBUG ) THEN
                            PRINT *,'POLYLINE SEGMENT :'
                            PRINT *,polyline%uv(1:2,sens,i1)
                            PRINT *,polyline%uv(1:2,sens,i2)
                            PRINT *,'SENS =',SENS
                            PRINT *,'DET =',vec(1)*aux(2) - vec(2)*aux(1)
                         END IF
                         if ( vec(1)*aux(2) - vec(2)*aux(1) > -EPSmath ) exit     !
                         if ( sens == 1 ) vec = -vec                              !
                         ! target point inside circular approximation but outside !
                         ! linear approximation (happens only if the domain       !
                         ! is locally convex)                                     !
                         wtmp = dot_product(aux, vec) / sum(vec**2)               !
                         IF ( DEBUG ) PRINT *,'WTMP =',WTMP
                         uvinter = (1._fp - wtmp)*polyline%uv(:,:,i1) + &         !
                              wtmp*polyline%uv(:,:,i2)                            !
                         p = uvinter(:,sens)                                      !
                         ! relax to exact intersection curve                      !
                         call simultaneous_point_inversions( &                    !
                              curve%surf, &                                       !
                              lowerb, &                                           !
                              upperb, &                                           !
                              stat, &                                             !
                              uvinter, &                                          !
                              xyzinter )                                          !
                         if ( stat > 0 ) then
                            ! failed to converge :(
                            PRINT *,'FAILED TO RELAX ON TANGENTIAL INTERSECTION'
                            PRINT *,'UVTMP  =',UVTMP
                            PRINT *,'DUVTMP =',DUVTMP
                            PRINT *,'UVPOLY ='
                            PRINT *,POLYLINE%UV(:,SENS,I1)
                            PRINT *,POLYLINE%UV(:,SENS,I2)
                            call debug_phf_simptinv( &
                                 curve, &
                                 (1._fp - wtmp)*polyline%uv(:,:,i1) + &
                                 wtmp*polyline%uv(:,:,i2) )
                            RETURN 
                         else
                            q = uvinter(:,sens)
                            wq = dot_product(q - polyline%uv(:,sens,i1), vec) &
                                 / sum(vec**2)
                            ! refine polyline
                            call insert_point_in_curve( &
                                 curve, &
                                 uvinter, &
                                 xyzinter, &
                                 i1, &
                                 debug, &
                                 stat )
                            !
                            IF ( DEBUG ) PRINT *,'WQ =',WQ
                            if ( wq < wtmp ) then ! <---+
                               i1 = i1 + 1              !
                               i2 = i2 + 1              !
                            end if ! <------------------+
                            pq = q - p
                            pu = uvtmp + duvtmp - p
                            IF ( DEBUG ) THEN
                               PRINT *,'U =',UVTMP
                               PRINT *,'D =',DUVTMP
                               PRINT *,'P =',P
                               PRINT *,'Q =',Q
                               PRINT *,'(PU.PQ)/(PQ.PQ) =',dot_product(pu,pq) / sum(pq**2)
                               PAUSE
                            END IF
                            if ( dot_product(pu,pq) > sum(pq**2) ) then
                               PRINT *,'U =',UVTMP
                               PRINT *,'D =',DUVTMP
                               PRINT *,'P =',P
                               PRINT *,'Q =',Q
                               PRINT *,'(PU.PQ)/(PQ.PQ) =',dot_product(pu,pq) / sum(pq**2)
                               PAUSE
                               inside = .false.
                               RETURN
                            end if
                         end if
                      end do
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
                      uvinter = (1._fp - wtmp)*polyline%uv(:,:,i1) + &
                           wtmp*polyline%uv(:,:,i2)
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
                         call debug_phf_simptinv( &
                                 curve, &
                                 (1._fp - wtmp)*polyline%uv(:,:,i1) + &
                                 wtmp*polyline%uv(:,:,i2) )
                         RETURN 
                      else
                         duvtmp = uvinter(1:2,sens) - uvtmp
                         IF ( DEBUG ) PRINT *,'MODIF DUVTMP =',duvtmp
                         wtmp = dot_product(duvtmp, vec) / sum(vec**2)
                         if ( w > 0.5_fp ) wtmp = 1._fp - wtmp
                         if ( wtmp > 1._fp ) then
                            i1 = i1 + 1
                         elseif ( wtmp < 0._fp ) then
                            i1 = i1 - 1
                         end if
                         ! refine polyline
                         call insert_point_in_curve( &
                              curve, &
                              uvinter, &
                              xyzinter, &
                              i1, &
                              debug, &
                              stat )
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
                if ( change_face ) then ! <----------------------------+
                   ! change brep face                                  !
                   ifacetmp = get_face(brep, get_twin(ihedg))          !
                   ! uv coordinates in new face                        !
                   uvtmp = uvinter(1:2,ihedg(2))                       !
                   ! new point displacement                            !
                   do ivar = 1,2 ! <-----------------------+           !
                      call evald1( &                       !           !
                           jac(:,ivar), &                  !           !
                           brep%faces(ifacetmp)%surface, & !           !
                           uvtmp, &                        !           !
                           ivar )                          !           !
                   end do ! <------------------------------+           !
                   call solve_NxN( &                                   !
                        duvtmp, &                                      !
                        matmul(transpose(jac), jac), &                 !
                        matmul(transpose(jac), xyztarget - xyztmp), &  !
                        singular )                                     !
                   dxyztmp = matmul(jac, duvtmp)                       !
                   IF ( DEBUG ) THEN ! -----------------------+        !
                      PRINT *,'CHANGE BREP FACE -->',IFACETMP !        !
                      PRINT *,'  UVTMP =',uvtmp               !        !
                      PRINT *,' DUVTMP =',duvtmp              !        !
                      PRINT *,' XYZTMP =',xyztmp              !        !
                      PRINT *,'DXYZTMP =',dxyztmp             !        !
                   END IF ! <---------------------------------+        !
                end if ! <---------------------------------------------+
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
    
  end subroutine projection_hyperface



  
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
             PAUSE
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

       IF ( DEBUG ) PRINT *,sqrt(sum(r**2) - rn**2), sqrt(erruv), EPSfp*cond

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
