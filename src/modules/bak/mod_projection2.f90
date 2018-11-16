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
    integer                                   :: iinter, jinter, i1, i2
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
             polyline => curve%polyline
             ! get current halfedge's uv polyline
             ! (oriented such that the interior of the face is on the left)
             ihead = curve%isplit(2,brep%edges(ihedg(1))%isplit)
             itail = curve%isplit(2,brep%edges(ihedg(1))%isplit + 1)
             np = itail- ihead + 1
             sens = 1 + mod(ihedg(2),2)
             allocate(uvpoly(2)%mat(2,np), uvpoly(3)%mat(2,np))
             uvpoly(2)%mat(1:2,1:np) = polyline%uv(1:2,sens,ihead:itail)
             uvpoly(3)%mat(1:2,1:np) = polyline%uv(1:2,ihedg(2),ihead:itail)
             allocate(xyzpoly(3,np))
             xyzpoly(1:3,1:np) = polyline%xyz(1:3,ihead:itail)
             surfpair(1)%ptr => curve%surf(sens)%ptr     ! check ordre
             surfpair(2)%ptr => curve%surf(ihedg(2))%ptr ! check ordre
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
             IF ( DEBUG ) THEN
                uvpoly(1)%mat(:,2) = uvtmp + duvtmp
                call debug_phf_polylines( &
                     uvpoly(1:2), &
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
                !iinter = 1
                !do jinter = 2,ninter
                !   if ( lamb(1,jinter) > lamb(1,iinter) - EPSfp .and. &
                !        lamb(2,jinter) < lamb(2,iinter) ) iinter = jinter
                !end do
                w = lamb(2,iinter)  ! local abscissa along intersected polyline segment
                i1 = ipls(2,iinter) ! origin vertex of the intersected polyline segment
                i2 = i1 + 1         ! destination vertex of the intersected polyline segment
                IF ( DEBUG ) THEN
                   PRINT *,'LAMBDA =',LAMB(:,IINTER)
                   PRINT *,'POLYLINE SEGMENT:'
                   PRINT *,'XYZ ='
                   PRINT *,xyzpoly(1:3,i1)
                   PRINT *,xyzpoly(1:3,i2)
                   PRINT *,'UV ='
                   PRINT *,uvpoly(2)%mat(1:2,i1)
                   PRINT *,uvpoly(2)%mat(1:2,i2)
                END IF
                ! vector from origin to destination
                vec = uvpoly(2)%mat(1:2,i2) - uvpoly(2)%mat(1:2,i1)
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
                   else ! ---------------------------------+                            !
                      w = 0._fp                            !                            !
                   end if ! <------------------------------+                            !
                   uvinter(1:2,sens) = &                                                !
                        (1._fp - w)*uvpoly(2)%mat(1:2,i1) + w*uvpoly(2)%mat(1:2,i2)     !
                   uvinter(1:2,ihedg(2)) = &                                            !
                        (1._fp - w)*uvpoly(3)%mat(1:2,i1) + w*uvpoly(3)%mat(1:2,i2)     !
                   xyztmp = (1._fp - w)*xyzpoly(1:3,i1) + w*xyzpoly(1:3,i2)             !
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
                   uvinter(1:2,sens) = &                                                !
                        (1._fp - w)*uvpoly(2)%mat(1:2,i1) + w*uvpoly(2)%mat(1:2,i2)     !
                   uvinter(1:2,ihedg(2)) = &                                            !
                        (1._fp - w)*uvpoly(3)%mat(1:2,i1) + w*uvpoly(3)%mat(1:2,i2)     !
                   xyztmp = (1._fp - w)*xyzpoly(1:3,i1) + w*xyzpoly(1:3,i2)             !
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
                   aux = uvpoly(2)%mat(j,[i1,i2])                                      !
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
                        surfpair, &                                                    !
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
                   duv_ds = sign(1._fp, dot_product(duv_ds(:,1,1), vec)) * duv_ds      !
                   
             end if
          end do halfedges
          
       end do wires
    end do faces
  

end module mod_projection
