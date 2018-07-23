recursive subroutine intersect_surface_pair( &
     surfroot, &
     region, &
     interdata, &
     uvxyz, &
     nuvxyz, &
     stat_degeneracy )
  use mod_math
  use mod_bernstein
  use mod_diffgeom
  use mod_regiontree
  use mod_tolerances
  use mod_types_intersection
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  type(ptr_surface),            intent(in)    :: surfroot(2)
  type(ptr_region),             intent(inout) :: region(2)
  type(type_intersection_data), intent(inout) :: interdata
  real(kind=fp), allocatable,   intent(inout) :: uvxyz(:,:)
  integer,                      intent(inout) :: nuvxyz
  integer,                      intent(inout) :: stat_degeneracy
  logical                                     :: overlap
  integer                                     :: stat_collineal
  real(kind=fp)                               :: uv_collineal(2,2)
  real(kind=fp), dimension(3)                 :: xyz_collineal, n_collineal
  integer                                     :: stat_loopdetection
  real(kind=fp)                               :: tmp(7)
  real(kind=fp)                               :: uv_subdiv(2,2)
  type(ptr_polynomial)                        :: poly(2)
  real(kind=fp)                               :: param_vector(3)
  integer                                     :: stat_singularpoint
  real(kind=fp)                               :: dxyz_duv(3,2,2), n(3,2), dxyz_ds(3,2), DUV_DS(2,2,2), CURVATURE(2)
  type(ptr_region)                            :: newregion(2)
  integer, dimension(2)                       :: stat_subdiv, nchild
  integer                                     :: isurf, ipt, ivar, ichild, jchild, ipoly

  if ( stat_degeneracy /= 0 ) return ! a degeneracy has been encountered

  IF ( DEBUG ) THEN
      PRINT *,''; PRINT *,'';
      PRINT *,'INTERSECT_SURFACE_PAIR'
      PRINT *,'UVBOXES ='
      DO ISURF = 1,2 ! <-----------------+
         PRINT *,REGION(ISURF)%PTR%UVBOX !
      END DO ! <-------------------------+
      PRINT *,'IPTS ='
      DO ISURF = 1,2
         IF ( REGION(ISURF)%PTR%NPTS < 1 ) THEN
            PRINT *,'N/A'
         ELSE
            PRINT *,REGION(ISURF)%PTR%IPTS(1:REGION(ISURF)%PTR%NPTS)
         END IF
      END DO
  END IF

  ! inherit from parent regions all intersection points contained in current regions
  if ( nuvxyz > 0 ) then ! <----------------------------+
     do isurf = 1,2 ! <---------------------------+     !
        call inherit_points( &                    !     !
             region(isurf)%ptr, &                 !     !
             uvxyz(2*isurf-1:2*isurf,1:nuvxyz), & !     !
             nuvxyz )                             !     !
     end do ! <-----------------------------------+     !
  end if ! <--------------------------------------------+


  ! compute bounding boxes for each region...
  do isurf = 1,2 ! <------------------------------------------------------+
     if ( .not.associated(region(isurf)%ptr%xyzbox) ) then ! <----+       !
        allocate( region(isurf)%ptr%xyzbox )                      !       !
        call bernOBB2( &                                          !       !
             region(isurf)%ptr%poly(1)%ptr%coef, &                !       !
             region(isurf)%ptr%poly(1)%ptr%degr, &                !       !
             region(isurf)%ptr%xyzbox )                           !       !
     end if ! <---------------------------------------------------+       !
  end do ! <--------------------------------------------------------------+

  ! ... and check if these bounding boxes overlap
  call overlap_OBBs( &
       region(1)%ptr%xyzbox, &
       region(2)%ptr%xyzbox, &
       overlap )

  if ( .not.overlap ) return ! disjoint bounding boxes => empty intersection

  ! check if there is a pair of collineal corner points
  call find_collineal_corners( &
       region, &
       stat_collineal, &
       uv_collineal, &
       n_collineal, &
       xyz_collineal )
  IF ( DEBUG ) PRINT *,'STAT_COLLINEAL(CORNERS) =',stat_collineal

  if ( stat_collineal <= 0 ) then ! <----------------------------+
     ! there is a pair of collineal corner points, we will       !
     ! subdivide both region at their parametric center point    !
     do isurf = 1,2 ! <----------------------------+             !
        uv_subdiv(:,isurf) = 0.5_fp * ( &          !             !
             region(isurf)%ptr%uvbox([1,3]) + &    !             !
             region(isurf)%ptr%uvbox([2,4]) )      !             !
     end do ! <------------------------------------+             !
  elseif ( stat_collineal < 2 ) then ! --------------------------+
     n_collineal(:) = 0._fp                                      !
  else ! --------------------------------------------------------+
     ! one surface has a singular corner                         !
     stat_degeneracy = 20                                        !
     return                                                      !
  end if ! <-----------------------------------------------------+


  ! check Hohmeyer's loop detection criterion
  do isurf = 1,2 ! <----------------------------------+
     poly(isurf)%ptr => region(isurf)%ptr%poly(2)%ptr !
  end do ! <------------------------------------------+
  call loop_detection_criterion( &
       poly, &
       stat_loopdetection, &
       param_vector, &
       ( stat_collineal <= 0 ), &
       n_collineal, &
       randomize=.false. )
  nullify(poly(1)%ptr, poly(2)%ptr)
  IF ( DEBUG ) PRINT *,'STAT_LOOPDETECTION =',stat_loopdetection

  if ( stat_loopdetection < 0 ) then ! <--------------+
     ! a degeneracy has been encountered              !
     stat_degeneracy = stat_loopdetection             !
     return                                           !
  end if ! <------------------------------------------+

  ! If the loop detection criterion IS satisfied, then the intersection of the 
  ! current surface regions is a (possibly empty) set of isolated points and/or 
  ! open curve branches. All the isolated points as well as the curves' endpoints
  ! are located on the boundary of a surface region.
  if ( stat_loopdetection == 0 ) then ! <---------------------------------------+
     ! first we start new temporary region trees, rooted at copies of current   !
     ! regions (all the data contained in the regions is copied, except the     !
     ! pointers to children)                                                    !
     do isurf = 1,2 ! <-------------------+                                     !
        allocate( newregion(isurf)%ptr )  !                                     !
        call copy_region( &               !                                     !
             region(isurf)%ptr, &         !                                     !
             newregion(isurf)%ptr )       !                                     !
     end do ! <---------------------------+                                     !
     !                                                                          !
     !IF ( DEBUG ) PRINT *,'BEFORE SIMPLE_SURFACES, NC =',INTERDATA%NC, &
     !     ', ALLOCATED?',ALLOCATED(INTERDATA%CURVES)
     call intersect_simple_surfaces( &                                          !
          surfroot, &                                                           !
          newregion, &                                                          !
          param_vector, &                                                       !
          interdata, &                                                          !
          uvxyz, &                                                              !
          nuvxyz, &                                                             !
          stat_degeneracy )                                                     !
     !IF ( DEBUG ) PRINT *,'BACK TO INTERSECT_SURFACE_PAIR'
     !                                                                          !
     do isurf = 1,2 ! <------------------------------------------------------+  !
        ! copy the data in temporay regions back to the current regions      !  !
        ! (essentially indices of newly discovered intersection points)...   !  !
        if ( newregion(isurf)%ptr%npts > 0 ) then ! <---------------------+  !  !
           !IF ( DEBUG ) PRINT *,'COPY NEW POINTS, SURF',ISURF
           call append_n( &                                               !  !  !
                region(isurf)%ptr%ipts, &                                 !  !  !
                region(isurf)%ptr%npts, &                                 !  !  !
                newregion(isurf)%ptr%ipts(1:newregion(isurf)%ptr%npts), & !  !  !
                newregion(isurf)%ptr%npts, &                              !  !  !
                unique=.true. )                                           !  !  !
           !IF ( DEBUG ) PRINT *,'OK'
        end if ! <--------------------------------------------------------+  !  !
        ! ...then free the temporay region trees...                          !  !
        !IF ( DEBUG ) PRINT *,'FREE MEMORY, SURF',ISURF
        nullify( &                                                           !  !
             newregion(isurf)%ptr%xyzbox,      &                             !  !
             newregion(isurf)%ptr%poly(1)%ptr, &                             !  !
             newregion(isurf)%ptr%poly(2)%ptr )                              !  !
        deallocate( newregion(isurf)%ptr%poly )                              !  !
        call free_region_tree( newregion(isurf)%ptr )                        !  !
        deallocate( newregion(isurf)%ptr )                                   !  !
        !IF ( DEBUG ) PRINT *,'OK'
     end do ! <--------------------------------------------------------------+  !
     !                                                                          !
     ! ... and finally return (there cannot be other intersection points/curves ! 
     ! between the current surface regions)                                     !
     return                                                                     !
     !                                                                          !
  end if ! <--------------------------------------------------------------------+


  if ( stat_collineal > 0 ) then ! <--------------------------------------------+
     ! There are no pair of collineal corners have, so we now                   !
     ! search for a pair of interior collineal points.                          !
     ! First, set initial iterate to the parametric center of both regions ...  !
     do isurf = 1,2 ! <------------------------------------------+              !
        uv_collineal(:,isurf) = 0.5_fp * ( &                     !              !
             region(isurf)%ptr%uvbox([1,3]) + &                  !              !
             region(isurf)%ptr%uvbox([2,4]) )                    !              !
     end do ! <--------------------------------------------------+              !
     ! ... then run box-constrained Newton-Raphson algorithm                    !
     call find_collineal_points( &                                              !
          surfroot, &                                                           !
          [region(1)%ptr%uvbox([1,3]), region(2)%ptr%uvbox([1,3])] - EPSuv, &   !
          [region(1)%ptr%uvbox([2,4]), region(2)%ptr%uvbox([2,4])] + EPSuv, &   !
          stat_collineal, &                                                     !
          uv_collineal, &                                                       !
          n_collineal, &                                                        !
          xyz_collineal )                                                       !
     IF ( DEBUG ) PRINT *,'STAT_COLLINEAL =',stat_collineal
     !                                                                          !
     ! if a pair of collineal points have been discovered, we will subdivide    !
     ! the surface regions at that point                                        !
     if ( stat_collineal <= 0 ) uv_subdiv = uv_collineal                        !
  end if ! <--------------------------------------------------------------------+


  if ( stat_collineal < 0 ) then  ! <-------------------------------------------+
     ! The collineal points are coincident: this is a tangential contact point. !
     ! first, add the point to the collection...                                !
     if ( nuvxyz > 0) then ! <------------+                                     !
        call check_unicity( &             !                                     !
             xyz_collineal, &             !                                     !
             3, &                         !                                     !
             uvxyz(5:7,1:nuvxyz), &       !                                     !
             nuvxyz, &                    !                                     !
             EPSxyz, &                    !                                     !
             ipt )                        !                                     !
     else ! ------------------------------!                                     !
        ipt = 1                           !                                     !
     end if ! <---------------------------+                                     !
     if ( ipt > nuvxyz ) then ! <--------------------------+                    !
        tmp(1:2) = uv_collineal(:,1)                       !                    !
        tmp(3:4) = uv_collineal(:,2)                       !                    !
        tmp(5:7) = xyz_collineal                           !                    !
        call append_vector( &                              !                    !
             tmp(1:7), &                                   !                    !
             7, &                                          !                    !
             uvxyz, &                                      !                    !
             nuvxyz )                                      !                    !
     end if ! <--------------------------------------------+                    !
     do isurf = 1,2 ! <------------------------+                                !
        call add_points_bottom_up( &           !                                !
             region(isurf)%ptr, &              !                                !
             [nuvxyz], &                       !                                !
             1 )                               !                                !
     end do ! <--------------------------------+                                !
     !                                                                          !
     ! ... then determine the nature of that singular point                     !
     ! compute tangent and normal vectors                                       !
     do isurf = 1,2 ! <----------------------------------+                      !
        do ivar = 1,2 ! <------------------+             !                      !
           call evald1( &                  !             !                      !
                dxyz_duv(:,ivar,isurf), &  !             !                      !
                surfroot(isurf)%ptr, &     !             !                      !
                uv_collineal(:,isurf), &   !             !                      !
                ivar )                     !             !                      !
        end do ! <-------------------------+             !                      !
        n(:,isurf) = cross( &                            !                      !
             dxyz_duv(:,1,isurf), dxyz_duv(:,2,isurf) )  !                      !
     end do ! <------------------------------------------+                      !
     call characterize_tangential_intersection_point( &                         !
          surfroot, &                                                           !
          uv_collineal, &                                                       !
          dxyz_duv, &                                                           !
          n, &                                                                  !
          stat_singularpoint, &                                                 !
          dxyz_ds )                                                             !
     !                                                                          !
     select case (stat_singularpoint) ! <-------------------------------+       !
        case (1) ! -----------------------------------------------------+       !
           ! point on tangential intersection curve                     !       !
           PRINT *,'>>> POINT ON TAGENTIAL INTSERSECTION CURVE'         !       !
           call diffgeom_intersection_curve( &
                surfroot, &
                uv_collineal, &
                duv_ds, &
                dxyz_ds, &
                stat_singularpoint, &
                curvature )
           PRINT *,'CURVATURE =',CURVATURE(1)
           PRINT *,'DXYZ_DS =',dxyz_ds(:,1)
           PRINT *,'DUV_DS='
           CALL PRINT_MAT(DUV_DS(:,1,:))
        case (2) ! -----------------------------------------------------+       !
           ! branch point (two curves meet at that point)               !       !
           stat_loopdetection = 3 ! force subdivision of both surfaces  !       !
           PRINT *,'>>> BRANCH POINT'                                   !       !
           call diffgeom_intersection_curve( &
                surfroot, &
                uv_collineal, &
                duv_ds, &
                dxyz_ds, &
                stat_singularpoint, &
                curvature )
           PRINT *,'CURVATURE =',CURVATURE
           PRINT *,'DXYZ_DS ='
           CALL PRINT_MAT(dxyz_ds)
           PRINT *,'DUV_DS,1='
           CALL PRINT_MAT(DUV_DS(:,1,:))
           PRINT *,'DUV_DS,2='
           CALL PRINT_MAT(DUV_DS(:,2,:))
        case (3) ! -----------------------------------------------------+       !
           ! isolated tangential contact point                          !       !
           stat_loopdetection = 3 ! force subdivision of both surfaces  !       !
           PRINT *,'>>> ISOLATED TANGENTIAL CONTACT POINT'              !       !
        case (4) ! -----------------------------------------------------+       !
           ! high-order contact point                                   !       !
           PRINT *,'>>> HIGH-ORDER CONTACT POINT'                       !       !
     end select ! <-----------------------------------------------------+       !
     !                                                                          !
  elseif ( stat_collineal > 0 ) then ! -----------------------------------------+
     ! no pair of collineal points has been found, we subdivide the surface     !
     ! with largest Gauss map at its center point                               !
     do isurf = 1,2 ! <------------------------------------------+              !
        uv_subdiv(:,isurf) = 0.5_fp * ( &                        !              !
             region(isurf)%ptr%uvbox([1,3]) + &                  !              !
             region(isurf)%ptr%uvbox([2,4]) )                    !              !
     end do ! <--------------------------------------------------+              !
     !                                                                          !
  elseif ( .FALSE. ) THEN !stat_collineal == 0 ) then ! ----------------------------------------+
     if ( stat_loopdetection /= 3 ) then
        isurf = stat_loopdetection
        do ivar = 1,2
           if ( is_in_open_interval( &
                uv_collineal(ivar,isurf), &
                region(isurf)%ptr%uvbox(2*ivar-1), &
                region(isurf)%ptr%uvbox(2*ivar), &
                tolerance=EPSregion ) ) exit
        end do
        if ( ivar > 2 ) then! <---------------------+
           uv_subdiv(:,isurf) = 0.5_fp * ( &        !
                region(isurf)%ptr%uvbox([1,3]) + &  !
                region(isurf)%ptr%uvbox([2,4]) )    !
        end if ! <----------------------------------+
     end if
  end if ! <--------------------------------------------------------------------+

  ! subdivide the surface regions
  nchild(:) = 1
  stat_subdiv(:) = 2
  do isurf = 1,2 ! <------------------------------------------------------------+
     if ( stat_loopdetection == isurf .or. stat_loopdetection == 3 ) then ! <-+ !
        call subdiv_region( &                                                 ! !
             region(isurf)%ptr, &                                             ! !
             uv_subdiv(:,isurf), &                                            ! !
             stat_subdiv(isurf) )                                             ! !
        !                                                                     ! !
        if ( stat_subdiv(isurf) == 2 ) then ! <------------------------+      ! !
           ! the subdivision point is at one of the surface's corners  !      ! !
           nchild(isurf) = 1                                           !      ! !
        else ! --------------------------------------------------------+      ! !
           ! the subdivision point is interior to the surface          !      ! !
           nchild(isurf) = size(region(isurf)%ptr%child)               !      ! !
        end if ! <-----------------------------------------------------+      ! !
     end if ! <---------------------------------------------------------------+ !
  end do ! <--------------------------------------------------------------------+

  if ( all(nchild < 2) ) then  ! <---------------------------------+
     IF ( .FALSE. ) THEN
        PRINT *,'STAT_COLLINEAL =',STAT_COLLINEAL
        PRINT *,'UVBOXES ='
        DO ISURF = 1,2 ! <-----------------+
           PRINT *,REGION(ISURF)%PTR%UVBOX !
        END DO ! <-------------------------+
        PRINT *,'UV_COLLINEAL =',UV_COLLINEAL
        STOP '%%%%%%%%%%%%%%%%%%%'
     END IF
     if ( stat_loopdetection == 3 ) then ! <------------------+    !
        ! recursion terminates prematurely because both       !    !
        ! regions have ranges smaller than 2*EPSregion        !    !
        stat_degeneracy = 21                                  !    !
        return                                                !    !
     else ! --------------------------------------------------+    !
        isurf = 3 - stat_loopdetection                        !    !
        call subdiv_region( &                                 !    !
             region(isurf)%ptr, &                             !    !
             uv_subdiv(:,isurf), &                            !    !
             stat_subdiv(isurf) )                             !    !
        if ( stat_subdiv(isurf) == 2 ) then ! <-----------+   !    !
           ! recursion terminates prematurely because     !   !    !
           ! both regions have ranges smaller than        !   !    !
           ! 2*EPSregion                                  !   !    !
           stat_degeneracy = 22                           !   !    !
           return                                         !   !    !
        else ! -------------------------------------------+   !    !
           nchild(isurf) = size(region(isurf)%ptr%child)  !   !    !
        end if ! <----------------------------------------+   !    !
     end if ! <-----------------------------------------------+    !
  end if ! <-------------------------------------------------------+


  ! compute Bezier control points for children regions
  do isurf = 1,2 ! <------------------------------------------------------------+
     if ( nchild(isurf) < 2 ) cycle ! no children                               !
     if ( stat_subdiv(isurf) < 0 ) cycle ! the children are already allocated   !
     !                                                                          !
     ! map uv_subdiv to the local frame ([0,1]^2) of the region                 !
     do ivar = 1,2 ! <--------------------------------+                         !
        uv_subdiv(ivar,isurf) = 0.5_fp * ( 1._fp + &  !                         !
             ab2n1p1( &                               !                         !
             uv_subdiv(ivar,isurf), &                 !                         !
             region(isurf)%ptr%uvbox(2*ivar-1), &     !                         !
             region(isurf)%ptr%uvbox(2*ivar) ) )      !                         !
     end do ! <---------------------------------------+                         !
     !                                                                          !
     ! allocate polynomials                                                     !
     do ichild = 1,nchild(isurf) ! <---------------------------+                !
        allocate(region(isurf)%ptr%child(ichild)%poly(2)    )  !                !
        allocate(region(isurf)%ptr%child(ichild)%poly(1)%ptr)  !                !
        allocate(region(isurf)%ptr%child(ichild)%poly(2)%ptr)  !                !
     end do ! <------------------------------------------------+                !
     !                                                                          !
     if ( stat_subdiv(isurf) == 0 ) then ! <------------------------------+     !
        ! 4 children                                                      !     !
        do ipoly = 1,2 ! <------------------------------------------+     !     !
           call subdiv_bezier2( &                                   !     !     !
                region(isurf)%ptr%poly(ipoly)%ptr, &                !     !     !
                uv_subdiv(:,isurf), &                               !     !     !
                bsw=region(isurf)%ptr%child(1)%poly(ipoly)%ptr, &   !     !     !
                bse=region(isurf)%ptr%child(2)%poly(ipoly)%ptr, &   !     !     !
                bnw=region(isurf)%ptr%child(3)%poly(ipoly)%ptr, &   !     !     !
                bne=region(isurf)%ptr%child(4)%poly(ipoly)%ptr )    !     !     !
        end do ! <--------------------------------------------------+     !     !
        !                                                                 !     !
     elseif ( stat_subdiv(isurf) == 1 ) then ! ---------------------------+     !
        ! 2 children                                                      !     !
        if ( region(isurf)%ptr%child(2)%uvbox(1) <= &                     !     !
             region(isurf)%ptr%uvbox(1) + EPSregion ) then ! <---------+  !     !
           ! the subdivision point is on a iso-u boundary              !  !     !
           do ipoly = 1,2 ! <---------------------------------------+  !  !     !
              call subdiv_bezier2_only_v( &                         !  !  !     !
                   region(isurf)%ptr%poly(ipoly)%ptr, &             !  !  !     !
                   v=uv_subdiv(2,isurf), &                          !  !  !     !
                   bs=region(isurf)%ptr%child(1)%poly(ipoly)%ptr, & !  !  !     !
                   bn=region(isurf)%ptr%child(2)%poly(ipoly)%ptr )  !  !  !     !
           end do ! <-----------------------------------------------+  !  !     !
           !                                                           !  !     !
        else ! --------------------------------------------------------+  !     !
           ! the subdivision point is on a iso-v boundary              !  !     !
           do ipoly = 1,2 ! <---------------------------------------+  !  !     !
              call subdiv_bezier2_only_u( &                         !  !  !     !
                   region(isurf)%ptr%poly(ipoly)%ptr, &             !  !  !     !
                   u=uv_subdiv(1,isurf), &                          !  !  !     !
                   bw=region(isurf)%ptr%child(1)%poly(ipoly)%ptr, & !  !  !     !
                   be=region(isurf)%ptr%child(2)%poly(ipoly)%ptr )  !  !  !     !
           end do ! <-----------------------------------------------+  !  !     !
           !                                                           !  !     !
        end if ! <-----------------------------------------------------+  !     !
        !                                                                 !     !
     end if ! <-----------------------------------------------------------+     !
     !                                                                          !
  end do ! <--------------------------------------------------------------------+


  ! carry on the recursion with pairs of children regions
  do jchild = 1,nchild(2) ! <--------------------------------------+
     if ( nchild(2) == 1 ) then ! <-------------------------+      !
        newregion(2)%ptr => region(2)%ptr                   !      !
     else ! ------------------------------------------------+      !
        newregion(2)%ptr => region(2)%ptr%child(jchild)     !      !
     end if ! <---------------------------------------------+      !
     !                                                             !
     do ichild = 1,nchild(1) ! <-------------------------------+   !
        if ( nchild(1) == 1 ) then ! <----------------------+  !   !
           newregion(1)%ptr => region(1)%ptr                !  !   !
        else ! ---------------------------------------------+  !   !
           newregion(1)%ptr => region(1)%ptr%child(ichild)  !  !   !
        end if ! <------------------------------------------+  !   !
        !                                                      !   !
        call intersect_surface_pair( &                         !   !
             surfroot, &                                       !   !
             newregion, &                                      !   !
             interdata, &                                      !   !
             uvxyz, &                                          !   !
             nuvxyz, &                                         !   !
             stat_degeneracy )                                 !   !
        !                                                      !   !
     end do ! <------------------------------------------------+   !
  end do ! <-------------------------------------------------------+

end subroutine intersect_surface_pair
