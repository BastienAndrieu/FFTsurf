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
    use mod_tolerances
    use mod_intersection
    implicit none
    logical, intent(inout) :: debug
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
    integer                                   :: it
    logical                                   :: change_face
    integer                                   :: iwire, ihedg(2), istart
    type(type_intersection_curve), pointer    :: curve => null()
    type(type_matrix)                         :: uvpoly(3)
    real(kind=fp), allocatable                :: xyzpoly(:,:)
    integer                                   :: ihead, itail, np, sens
    integer                                   :: ninter
    integer, allocatable                      :: ipls(:,:)
    real(kind=fp), allocatable                :: lamb(:,:)
    integer                                   :: iinter, jpls, kpls, lpls
    real(kind=fp)                             :: apls
    real(kind=fp)                             :: uvinter(2,2)
    logical                                   :: on_border
    real(kind=fp)                             :: aux(2)
    integer                                   :: ivar
    integer                                   :: stat_tangent
    real(kind=fp)                             :: duv_ds(2,2,2), dxyz_ds(3,2)
    real(kind=fp), dimension(2)               :: vecjk, vecjo
    integer                                   :: stat_circ
    integer                                   :: subit
    real(kind=fp)                             :: aplstmp
    integer                                   :: stat_ptinv, stat_insert
    real(kind=fp)                             :: deltar, det
    logical                                   :: inside
    real(kind=fp)                             :: jac(3,2)
    logical                                   :: singular
    logical                                   :: skipnext
    
    IF ( DEBUG ) PRINT *,'PROJECTION HYPERFACE'

    upperb(1:4) = 1._fp + EPSuv
    lowerb = -upperb

    ifacetmp = iface
    uvtmp = uv
    stat_proj = 1

    if ( norm2(duv) < EPSuv ) then
       IF ( DEBUG ) THEN
          PRINT *,'|DUV| << 1'
          PAUSE
       END IF
       uvtmp = uv
       stat_proj = 0
       return
    end if

    duvtmp = duv
    xyztmp = xyz
    dxyztmp = dxyz

    xyztarget = xyz + dxyz

    allocate(uvpoly(1)%mat(2,2))

    it = 0
    faces : do
       IF ( DEBUG ) PRINT *,'FACE',ifacetmp
       ! the provisional point uvtmp is inside the parametric domain of
       ! the face ifacetmp (possibly on its boundary)
       ! check whether the point uvtmp + duvtmp is inside/outside this domain
       !uvpoly(1)%mat(:,1) = uvtmp
       !uvpoly(1)%mat(:,2) = uvtmp + lambmax*duvtmp

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
          skipnext = .false.
          halfedges : do
             uvpoly(1)%mat(:,1) = uvtmp
             uvpoly(1)%mat(:,2) = uvtmp + lambmax*duvtmp
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
             IF ( DEBUG ) PRINT *,'SENS =',sens
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
             if ( ninter > 0 ) lamb(1,:) = lamb(1,:) / lambmax ! renormalize lambda
             !
             IF ( DEBUG ) THEN
                call debug_phf_polylines( &
                     uvpoly(1:2), &
                     ninter, &
                     ipls, &
                     lamb )
                PRINT *,'UVTMP  =',UVTMP
                PRINT *,'DUVTMP =',DUVTMP
                PRINT *,'--> CHECK POLYLINES'
                PAUSE
             END IF
             !                                                                                            
             if ( ninter > 0 ) then ! <----------------------------------------------------------------------+
                ! pick intersection point closest to target point                                            !
                iinter = maxloc(lamb(1,1:ninter),1)                                                          !
                apls = lamb(2,iinter) ! local abscissa along intersected polyline segment                    !
                ! if this point is in the interior of a polyline segment, snap it to the                     !
                ! nearest polyline vertex (with index jpls)                                                  !
                jpls = ipls(2,iinter)                                                                        !
                if ( apls > 0.5_fp ) jpls = jpls + 1                                                         !
                uvtmp = (1._fp - apls)*uvpoly(2)%mat(1:2,ipls(2,iinter)) + &                                 !
                     apls*uvpoly(2)%mat(1:2,ipls(2,iinter)+1)                                                !
                IF ( DEBUG ) PRINT *,'UVPOLY =',uvtmp
                xyztmp = xyzpoly(1:3,jpls)                                                                   !
                uvinter(1:2,sens) = uvpoly(2)%mat(1:2,jpls)                                                  !
                uvinter(1:2,ihedg(2)) = uvpoly(3)%mat(1:2,jpls)                                              !
                ! correct displacement                                                                       !
                duvtmp = uvtmp + duvtmp - uvinter(1:2,sens)                                                  !
                uvtmp = uvinter(1:2,sens)                                                                    !
                do ivar = 1,2 ! <-----------------------+                                                    !
                   call evald1( &                       !                                                    !
                        jac(:,ivar), &                  !                                                    !
                        brep%faces(ifacetmp)%surface, & !                                                    !
                        uvtmp, &                        !                                                    !
                        ivar )                          !                                                    !
                end do ! <------------------------------+                                                    !
                call solve_NxN( &                                                                            !
                     duvtmp, &                                                                               !
                     matmul(transpose(jac), jac), &                                                          !
                     matmul(transpose(jac), xyztarget - xyztmp), &                                           !
                     singular )                                                                              !
                dxyztmp = matmul(jac, duvtmp)                                                                !
                !dxyztmp = xyztarget - xyztmp                                                                 !
                !                                                                                            !
                IF ( DEBUG ) THEN
                   PRINT *,'UVTMP* =', UVTMP
                   PRINT *,'DUVTMP  =',duvtmp
                END IF
                !                                                                                            !
                ! special treatment if the point is on a natural border of the parametric square             !
                ! (prevents robustness issues in tangential surface intersections)                           !
                on_border = .false.
                stat_circ = 0                                                                                !
                do ivar = 1,2 ! <---------------------------------------------------------------------+      !
                   aux(1) = uvtmp(ivar)                                                               !      !
                   kpls = jpls + 1                                                                    !      !
                   if ( jpls == np ) kpls = np - 1                                                    !      !
                   aux(2) = uvpoly(2)%mat(ivar,kpls)                                                  !      !
                   if ( all(abs(aux) > 1._fp - EPSuv) .and. product(aux) > 0._fp ) then ! <---------+ !      !
                      vecjk = uvpoly(2)%mat(1:2,kpls) - uvpoly(2)%mat(1:2,jpls)                     ! !      !
                      if ( kpls < jpls ) vecjk = -vecjk                                             ! !      !
                      duv_ds(1:2,1,sens) = vecjk                                                    ! !      !
                      ! pick correct polyline segment                                               ! !      !
                      if ( dot_product(duv_ds(1:2,1,sens), duvtmp) > 0._fp ) then ! <--+            ! !      !
                         kpls = jpls + 1                                               !            ! !      !
                      else ! ----------------------------------------------------------+            ! !      !
                         kpls = jpls - 1                                               !            ! !      !
                      end if ! <-------------------------------------------------------+            ! !      !
                      if ( jpls == np ) kpls = np - 1                                               ! !      !
                      if ( jpls == 1 ) kpls = 2                                                     ! !      !
                      vecjk = uvpoly(2)%mat(1:2,kpls) - uvpoly(2)%mat(1:2,jpls)                     ! !      !
                      !                                                                             ! !      !
                      on_border = .true.
                      stat_circ = 2                                                                 ! !      !
                      vecjo = [-duv_ds(2,1,sens), duv_ds(1,1,sens)]                                 ! !      !
                      exit                                                                          ! !      !
                   end if ! <-----------------------------------------------------------------------+ !      !
                end do ! <----------------------------------------------------------------------------+      !
                !                                                                                            !
                if ( stat_circ == 0 ) then ! <--------------------------------------------------------+      !
                   ! compute the exact tangent direction to the curve                                 !      !
                   call diffgeom_intersection( &                                                      !      !
                        curve%surf, &                                                                 !      !
                        uvinter, &                                                                    !      !
                        duv_ds, &                                                                     !      !
                        dxyz_ds, &                                                                    !      !
                        stat_tangent )                                                                !      !
                   ! pick correct orientation (incident face's interior on the left)                  !      !
                   kpls = jpls + 1                                                                    !      !
                   if ( jpls == np ) kpls = np - 1                                                    !      !
                   vecjk = uvpoly(2)%mat(1:2,kpls) - uvpoly(2)%mat(1:2,jpls)                          !      !
                   if ( kpls < jpls ) vecjk = -vecjk                                                  !      !
                   if ( dot_product(duv_ds(1:2,1,sens), vecjk) < 0._fp ) duv_ds = -duv_ds             !      !
                   !                                                                                  !      !
                   ! pick correct polyline segment                                                    !      !
                   if ( dot_product(duv_ds(1:2,1,sens), duvtmp) > 0._fp ) then ! <------------------+ !      !
                      kpls = jpls + 1                                                               ! !      !
                      !if ( jpls == np ) goto 100
                      if ( jpls == np ) then
                         deallocate(uvpoly(2)%mat, uvpoly(3)%mat, xyzpoly)
                         ! move on to the next halfedge on the wire
                         ihedg = get_next(brep, ihedg)
                         IF ( DEBUG ) PRINT *,'MOVE ON TO NEXT HALFEDGE...'
                         if ( get_orig(brep, ihedg) == istart ) exit halfedges
                         cycle halfedges
                      end if
                   else ! --------------------------------------------------------------------------+ !      !
                      kpls = jpls - 1                                                               ! !      !
                      !if ( jpls == 1 ) goto 100
                      if ( jpls == 1 ) then
                         deallocate(uvpoly(2)%mat, uvpoly(3)%mat, xyzpoly)
                         skipnext = .true.
                         ! return to the previous halfedge on the wire
                         ihedg = get_prev(brep, ihedg)
                         IF ( DEBUG ) PRINT *,'BACK TO PREVIOUS HALFEDGE...'
                         !DEBUG = .TRUE.
                         !PAUSE
                         cycle halfedges
                      end if
                   end if ! <-----------------------------------------------------------------------+ !      !
                   if ( jpls == np ) kpls = np - 1                                                    !      !
                   if ( jpls == 1 ) kpls = 2                                                          !      !
                   vecjk = uvpoly(2)%mat(1:2,kpls) - uvpoly(2)%mat(1:2,jpls)                          !      !
                   !                                                                                  !      !
                   ! use circular approximation                                                       !      !
                   call circular_approximation( &                                                     !      !
                        duv_ds(1:2,1,sens), &                                                         !      !
                        vecjk, &                                                                      !      !
                        vecjo, &                                                                      !      !
                        stat_circ )                                                                   !      !
                end if ! <----------------------------------------------------------------------------+      !
                IF ( DEBUG ) THEN
                   PRINT *,'ON BORDER?',ON_BORDER
                   PRINT *,'TANGENT =',DUV_DS(1:2,1,SENS)
                   PRINT *,'TNG.DUV =',dot_product(duv_ds(1:2,1,sens), duvtmp)
                   PRINT *,'J, K    =',JPLS, KPLS,'/',NP
                   PRINT *,'VECJK   =',VECJK
                   PRINT *,'VECJO   =',VECJO
                   PAUSE
                END IF
                !                                                                                            !
                if ( stat_circ < 2 ) then ! <--------------------------------------------------------------+ !
                   ! finite radius                                                                         ! !
                   deltar = sum((duvtmp - vecjo)**2) - sum(vecjo**2)                                       ! !
                   if ( stat_circ == 1 ) deltar = -deltar ! domain boundary locally concave                ! !
                   inside = ( deltar < -EPSmath )                                                          ! !
                   !                                                                                       ! !
                   if ( inside ) then ! <----------------------------------------------------------------+ ! !
                      subit = 0                                                                          ! ! !
                      do ! <---------------------------------------------------------------------------+ ! ! !
                         aux = uvtmp + duvtmp - uvpoly(2)%mat(1:2,jpls)                                ! ! ! ! 
                         vecjk = uvpoly(2)%mat(1:2,kpls) - uvpoly(2)%mat(1:2,jpls)                     ! ! ! !
                         det = vecjk(1)*aux(2) - vecjk(2)*aux(1)                                       ! ! ! !
                         if ( kpls < jpls ) det = -det                                                 ! ! ! !
                         IF ( DEBUG ) THEN
                            PRINT *,'J, K =',JPLS, KPLS
                            PRINT *,'AUX  =',AUX
                            PRINT *,'VJK  =', VECJK
                            PRINT *,'DET  =',DET
                            PAUSE
                         END IF
                         if ( det > EPSmath ) exit                                                     ! ! ! !
                         ! target point inside circular approximation but outside linear approx.       ! ! ! !
                         ! (happens only if the domain is locally convex)                              ! ! ! !
                         aplstmp = dot_product(aux, vecjk) / sum(vecjk**2)                             ! ! ! !
                         uvinter(1:2,sens) = (1._fp - aplstmp)*uvpoly(2)%mat(1:2,jpls) + &             ! ! ! !
                              aplstmp*uvpoly(2)%mat(1:2,kpls)                                          ! ! ! !
                         uvinter(1:2,ihedg(2)) = (1._fp - aplstmp)*uvpoly(3)%mat(1:2,jpls) + &         ! ! ! !
                              aplstmp*uvpoly(3)%mat(1:2,kpls)                                          ! ! ! !
                         aux = uvinter(1:2,sens)                                                       ! ! ! !
                         ! relax to exact intersection curve                                           ! ! ! !
                         call simultaneous_point_inversions( &                                         ! ! ! !
                              curve%surf, &                                                            ! ! ! !
                              lowerb, &                                                                ! ! ! !
                              upperb, &                                                                ! ! ! !
                              stat_ptinv, &                                                            ! ! ! !
                              uvinter, &                                                               ! ! ! !
                              xyztmp )                                                                 ! ! ! !
                         !                                                                             ! ! ! !
                         if ( stat_ptinv > 0 ) then ! <----------------------------------------------+ ! ! ! !
                            PRINT *,'projection_hyperface: FAILED TO RELAX TANGENTIAL INTERSECTION'  ! ! ! ! !
                            uvinter(1:2,sens) = (1._fp - aplstmp)*uvpoly(2)%mat(1:2,jpls) + &        ! ! ! ! !
                                 aplstmp*uvpoly(2)%mat(1:2,kpls)                                     ! ! ! ! !
                            uvinter(1:2,ihedg(2)) = (1._fp - aplstmp)*uvpoly(3)%mat(1:2,jpls) + &    ! ! ! ! !
                                 aplstmp*uvpoly(3)%mat(1:2,kpls)                                     ! ! ! ! !
                            PRINT *,'UVPOLY =',uvinter                                               ! ! ! ! ! 
                            call debug_phf_simptinv( &                                               ! ! ! ! !
                                 curve, &                                                            ! ! ! ! !
                                 uvinter )                                                           ! ! ! ! !
                            PAUSE
                            RETURN                                                                   ! ! ! ! !
                         else ! ---------------------------------------------------------------------+ ! ! ! !
                            apls = dot_product(uvinter(1:2,sens) - uvpoly(2)%mat(1:2,jpls),&
                                 vecjk) / sum(vecjk**2)     ! ! ! ! !
                            IF ( DEBUG ) PRINT *,'ATMP, A =',APLSTMP, APLS
                            !apls = dot_product(uvinter(1:2,sens) - uvtmp, vecjk) / sum(vecjk**2)     ! ! ! ! !
                            ! refine polyline                                                        ! ! ! ! !
                            if ( sens == 1 ) then
                               lpls = np - max(jpls, kpls) + 1
                            else
                               lpls = min(jpls, kpls)
                            end if
                            call insert_point_in_curve( &                                            ! ! ! ! !
                                 curve, &                                                            ! ! ! ! !
                                 uvinter, &                                                          ! ! ! ! !
                                 xyztmp, &                                                           ! ! ! ! !
                                 lpls, &                                                             ! ! ! ! !
                                 debug, &                                                            ! ! ! ! !
                                 stat_insert )                                                       ! ! ! ! !
                            deallocate(uvpoly(2)%mat, uvpoly(3)%mat, xyzpoly)
                            call get_polyline_endpoints( &
                                 brep, &
                                 ihedg, &
                                 ihead, &
                                 itail, &
                                 sens, &
                                 np )
                            allocate(uvpoly(2)%mat(2,np), uvpoly(3)%mat(2,np), xyzpoly(3,np))
                            uvpoly(2)%mat(1:2,1:np) = curve%polyline%uv(1:2,sens,ihead:itail:(-1)**sens)
                            uvpoly(3)%mat(1:2,1:np) = curve%polyline%uv(1:2,ihedg(2),ihead:itail:(-1)**sens)
                            xyzpoly(1:3,1:np) = curve%polyline%xyz(1:3,ihead:itail:(-1)**sens)
                            !                                                                        ! ! ! ! !
                            if ( apls > aplstmp ) then ! <----+                                      ! ! ! ! !
                               if ( jpls > kpls ) then ! <--+ !                                      ! ! ! ! !
                                  jpls = jpls + 1           ! !                                      ! ! ! ! !
                                  kpls = jpls - 1           ! !                                      ! ! ! ! !
                               end if ! <-------------------+ !                                      ! ! ! ! !
                            else ! ---------------------------+                                      ! ! ! ! !
                               if ( jpls < kpls ) then ! <--+ !                                      ! ! ! ! !
                                  jpls = jpls + 1           ! !                                      ! ! ! ! !
                                  kpls = jpls - 1           ! !                                      ! ! ! ! !
                               end if ! <-------------------+ !                                      ! ! ! ! !
                            end if ! <------------------------+                                      ! ! ! ! !
                            !                                                                        ! ! ! ! !
                            if ( dot_product(uvtmp + duvtmp - aux, uvinter(1:2,sens) - aux) > &      ! ! ! ! !
                                 sum((uvinter(1:2,sens) - aux)**2) ) then ! <----------------------+ ! ! ! ! !
                               PRINT *, '(PU.PQ)/(PQ.PQ) =', &                                     ! ! ! ! ! !
                                    dot_product(uvtmp + duvtmp - aux, uvinter(1:2,sens) - aux)     ! ! ! ! ! !
                               PAUSE                                                               ! ! ! ! ! !
                               inside = .false.                                                    ! ! ! ! ! !
                               return                                                              ! ! ! ! ! !
                            end if ! <-------------------------------------------------------------+ ! ! ! ! !
                         end if ! <------------------------------------------------------------------+ ! ! ! !
                         subit = subit + 1                                                             ! ! ! !
                         if ( subit > 10 ) RETURN                                                      ! ! ! !
                      end do ! <-----------------------------------------------------------------------+ ! ! !
                      !                                                                                  ! ! !
                   end if ! <----------------------------------------------------------------------------+ ! !
                   !                                                                                       ! !
                   if ( .not.inside .and. deltar < 5.d-2 * norm2(duvtmp) ) then ! <----------------------+ ! !
                      IF ( DEBUG ) PRINT *,'NEARLY TANGENTIAL DISPLACEMENT'                              ! ! !
                      IF ( DEBUG ) PRINT *,'DR/R =',deltar/norm2(vecjo),', DR/|DUV| =',deltar/norm2(duvtmp)
                      ! the target point is just outside the local circular approximation of the domain  ! ! !
                      ! (the displacement is nearly tangent to the circular boundary)                    ! ! !
                      aplstmp = dot_product(duvtmp, vecjk) / sum(vecjk**2)                               ! ! !
                      uvinter(1:2,sens) = (1._fp - aplstmp)*uvpoly(2)%mat(1:2,jpls) + &                  ! ! !
                           aplstmp*uvpoly(2)%mat(1:2,kpls)                                               ! ! !
                      uvinter(1:2,ihedg(2)) = (1._fp - aplstmp)*uvpoly(3)%mat(1:2,jpls) + &              ! ! !
                           aplstmp*uvpoly(3)%mat(1:2,kpls)                                               ! ! !
                      aux = uvinter(1:2,sens)                                                            ! ! !
                      ! relax to exact intersection curve                                                ! ! !
                      call simultaneous_point_inversions( &                                              ! ! !
                           curve%surf, &                                                                 ! ! !
                           lowerb, &                                                                     ! ! !
                           upperb, &                                                                     ! ! !
                           stat_ptinv, &                                                                 ! ! !
                           uvinter, &                                                                    ! ! !
                           xyztmp )                                                                      ! ! !
                      if ( stat_ptinv > 0 ) then ! <---------------------------------------------------+ ! ! !
                         PRINT *,'projection_hyperface: FAILED TO RELAX ON TANGENTIAL INTERSECTION'    ! ! ! !
                         uvinter(1:2,sens) = (1._fp - aplstmp)*uvpoly(2)%mat(1:2,jpls) + &             ! ! ! !
                              aplstmp*uvpoly(2)%mat(1:2,kpls)                                          ! ! ! !
                         uvinter(1:2,ihedg(2)) = (1._fp - aplstmp)*uvpoly(3)%mat(1:2,jpls) + &         ! ! ! !
                              aplstmp*uvpoly(3)%mat(1:2,kpls)                                          ! ! ! !
                         PRINT *,'UVPOLY =',uvinter                                                    ! ! ! ! 
                         call debug_phf_simptinv( &                                                    ! ! ! !
                              curve, &                                                                 ! ! ! !
                              uvinter )                                                                ! ! ! !
                         PAUSE
                         RETURN                                                                        ! ! ! !
                      else ! --------------------------------------------------------------------------+ ! ! !
                         duvtmp = uvinter(1:2,sens) - uvtmp                                            ! ! ! !
                         !duvtmp(:) = 0._fp!uvtmp + duvtmp - uvinter(1:2,sens)
                         !uvtmp = uvinter(1:2,sens)
                         apls = dot_product(duvtmp, vecjk) / sum(vecjk**2)                             ! ! ! !
                         if ( apls < 0._fp ) then ! <-----+                                            ! ! ! !
                            jpls = jpls - 1               !                                            ! ! ! !
                            kpls = kpls - 1               !                                            ! ! ! !
                         elseif ( apls > 1._fp ) then ! --+                                            ! ! ! !
                            jpls = jpls + 1               !                                            ! ! ! !
                            kpls = kpls + 1               !                                            ! ! ! !
                         end if ! <-----------------------+                                            ! ! ! !
                         if ( sens == 1 ) then ! <--------+                                            ! ! ! !
                            jpls = np - jpls + 1          !                                            ! ! ! !
                            kpls = np - kpls + 1          !                                            ! ! ! !
                         end if ! <-----------------------+                                            ! ! ! !
                         ! refine polyline                                                             ! ! ! !
                         call insert_point_in_curve( &                                                 ! ! ! !
                              curve, &                                                                 ! ! ! !
                              uvinter, &                                                               ! ! ! !
                              xyztmp, &                                                                ! ! ! !
                              min(jpls, kpls), &                                                       ! ! ! !
                              debug, &                                                                 ! ! ! !
                              stat_insert )                                                            ! ! ! !
                         inside = .true.                                                               ! ! ! !
                         !RETURN!****
                      end if ! <-----------------------------------------------------------------------+ ! ! !
                   end if ! <----------------------------------------------------------------------------+ ! !
                   !                                                                                       ! !
                else ! ------------------------------------------------------------------------------------+ !
                   ! infinite radius                                                                       ! !
                   !inside = ( dot_product(vecjo, duvtmp) > EPSmath )                                       ! !
                   vecjo = vecjo / hypot(vecjo(1), vecjo(2))
                   det = dot_product(vecjo, duvtmp)
                   IF ( DEBUG ) PRINT *,'DOT =',det
                   if ( abs(det) < EPSuv ) then
                      if ( det < 0._fp ) then
                         duvtmp = duvtmp - 2._fp*det*vecjo
                         det = dot_product(vecjo, duvtmp) ! -det
                         dxyztmp = matmul(jac, duvtmp)
                         xyztarget = xyztmp + dxyztmp
                         IF ( DEBUG ) PRINT *,'CORRECTION DUVTMP, DOT =',dot_product(vecjo, duvtmp)
                      end if
                      inside = .true.
                   end if
                   inside = ( det > 0._fp )                      
                end if ! <---------------------------------------------------------------------------------+ !
                !                                                                                            !
                IF ( DEBUG ) PRINT *,'INSIDE?',INSIDE                                                        !
                change_face = .not.inside                                                                    !
                !                                                                                            !
                if ( change_face ) then ! <----------------------------------------------------------------+ !
                   ! change brep face                                                                      ! !
                   ifacetmp = get_face(brep, get_twin(ihedg))                                              ! !
                   ! uv coordinates in new face                                                            ! !
                   uvtmp = uvinter(1:2,ihedg(2))                                                           ! !
                   ! new point displacement                                                                ! !
                   do ivar = 1,2 ! <-----------------------+                                               ! !
                      call evald1( &                       !                                               ! !
                           jac(:,ivar), &                  !                                               ! !
                           brep%faces(ifacetmp)%surface, & !                                               ! !
                           uvtmp, &                        !                                               ! !
                           ivar )                          !                                               ! !
                   end do ! <------------------------------+                                               ! !
                   call solve_NxN( &                                                                       ! !
                        duvtmp, &                                                                          ! !
                        matmul(transpose(jac), jac), &                                                     ! !
                        matmul(transpose(jac), xyztarget - xyztmp), &                                      ! !
                        singular )                                                                         ! !
                   dxyztmp = matmul(jac, duvtmp)                                                           ! !
                   IF ( DEBUG ) THEN ! -----------------------+                                            ! !
                      PRINT *,'CHANGE BREP FACE -->',IFACETMP !                                            ! !
                      PRINT *,'  UVTMP =',uvtmp               !                                            ! !
                      PRINT *,' DUVTMP =',duvtmp              !                                            ! !
                      PRINT *,' XYZTMP =',xyztmp              !                                            ! !
                      PRINT *,'DXYZTMP =',dxyztmp             !                                            ! !
                   END IF ! <---------------------------------+                                            ! !
                end if ! <---------------------------------------------------------------------------------+ !
                IF ( DEBUG ) PAUSE                                                                           !
                !                                                                                            !
             end if ! <--------------------------------------------------------------------------------------+
             !
!100          deallocate(uvpoly(2)%mat, uvpoly(3)%mat, xyzpoly)
             deallocate(uvpoly(2)%mat, uvpoly(3)%mat, xyzpoly)
             if ( change_face ) exit wires
             !
             ! move on to the next halfedge on the wire
             ihedg = get_next(brep, ihedg)
             if ( get_orig(brep, ihedg) == istart ) exit halfedges
             if ( skipnext ) then
                ihedg = get_next(brep, ihedg)
                skipnext = .false.
                IF ( DEBUG ) PRINT *,'SKIP HALFEDGE...'
                if ( get_orig(brep, ihedg) == istart ) exit halfedges
             end if
          end do halfedges

       end do wires

       if ( .not.change_face ) then
          IF ( DEBUG ) PRINT *,'CONVERGED, IT =',it
          IF ( DEBUG ) PRINT *,'UVNEW =',uvtmp + duvtmp
          stat_proj = 0
          exit faces
       end if

       it = it + 1
       if ( it > 2*brep%nf ) then
          PRINT *,'*** FAILED TO PROJECT ON HYPERFACE ***'
          duvtmp(1:2) = 0._fp
          ifacetmp = iface
          uvtmp = uv
          exit faces
       end if

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
       xyztmp, &
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
    real(kind=fp),        intent(out)         :: xyztmp(3)
    integer,              intent(out)         :: stat_proj
    real(kind=fp)                             :: duvtmp(2,2)
    real(kind=fp), dimension(3)               :: xyzprev, xyztarget, dxyztmp, p
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
       IF ( ifirst < 0 .OR. ifirst > size(polyline%xyz,2) .or. &
         ilast < 0 .OR. ilast > size(polyline%xyz,2) ) THEN
         PRINT *,'BREP EDGE#', iedgetmp, ', POLYLINE HEAD, TAIL:', ifirst, ilast
       END IF
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
    xyztmp = xyz
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
       IF ( DEBUG ) PRINT *,'NEW POINT INSERTED IN POLYLINE : UV =',uvinter, ', IVERT =',iprev+1
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
    integer                           :: stat_refl

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
            2, &
            stat_refl )
       if ( stat_refl > 0 ) return

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
