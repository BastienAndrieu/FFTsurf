recursive subroutine intersect_surface_surface( &
     surfroot, &
     region, &
     interdat, &
     stat_degeneracy ) 
  use mod_math
  use mod_bernstein2
  use mod_diffgeom2
  use mod_regiontree
  use mod_tolerances
  implicit none
  LOGICAL :: DEBUG = .true.
  type(ptr_surface),            intent(in)    :: surfroot(2)
  type(ptr_region),             intent(inout) :: region(2)
  type(type_intersection_data), intent(inout) :: interdat
  integer,                      intent(inout) :: stat_degeneracy
  logical                                     :: overlap
  real(kind=fp)                               :: param_vector(3)
  real(kind=fp), allocatable                  :: uvxyz(:,:)
  integer                                     :: nuvxyz
  integer                                     :: stat_loopdetection
  real(kind=fp)                               :: uv_collineal(2,2), xyz_collineal(3), normal_collineal(3)!, nor(3)
  integer                                     :: stat_collineal
  logical                                     :: separable
  real(kind=fp), dimension(4)                 :: lowerb, upperb, ranges
  integer                                     :: stat_singularpoint
  real(kind=fp)                               :: dxyz_duv(3,2,2), n(3,2), dxyz_ds(3,2)
  real(kind=fp)                               :: uv_subdiv(2,2)
  integer                                     :: stat_subdiv, nchild(2)
  type(ptr_region)                            :: newregion(2)
  integer                                     :: isurf, ivar, ichild, jchild, idpt
  CHARACTER                                   :: STR1, STR2
  logical                                     :: isonboundary(2)

  if ( stat_degeneracy > 1 ) return

  IF ( DEBUG ) THEN
     PRINT *,''; PRINT *,''; PRINT *,''
     PRINT *,'UVBOXES ='
     DO ISURF = 1,2
        PRINT *,REGION(ISURF)%PTR%UVBOX
     END DO
  END IF

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
       uv_collineal, &
       normal_collineal, &
       stat_collineal )
  IF ( DEBUG ) THEN
     PRINT *,'COLLINEAL CORNERS?', (STAT_COLLINEAL <= 0)
     IF ( STAT_COLLINEAL <= 0 ) THEN
        PRINT *,' UV =',uv_collineal
     END IF
  END IF

  if ( stat_collineal < 0 ) then
     call eval( &
          xyz_collineal, &
          surfroot(1)%ptr, &
          uv_collineal(:,1) )
     call intersect_surface_surface_elsewhere( &
          region(1)%ptr%poly(1)%ptr, &
          region(2)%ptr%poly(1)%ptr, &
          xyz_collineal, &
          separable, &
          randomize=.true. )
     PRINT *,'SEPARABLE ?',SEPARABLE
     !IF ( .NOT.SEPARABLE ) STOP 'debug'
     if ( separable ) return
  end if
  !IF ( STAT_COLLINEAL > 0 ) THEN
  !   PRINT *,'NO PAIR OF COLLINEAL CORNERS'
  !END IF

  ! use Hohmeyer's loop detection technique
  call hohmeyer_loop_detection( &
       region(1)%ptr%poly(2)%ptr, &
       region(2)%ptr%poly(2)%ptr, &
       param_vector, &
       stat_loopdetection, &
       ( stat_collineal <= 0 ),  &
       normal_collineal )

  IF ( DEBUG ) PRINT *,'STAT_LOOPDETECTION =',STAT_LOOPDETECTION
  IF ( STAT_LOOPDETECTION < 0 ) THEN
     DO ISURF = 1,2
        PRINT *,'SURF #',ISURF
        PRINT *,'UVBOX =',REGION(ISURF)%PTR%UVBOX
        WRITE (STR1,'(I1)') ISURF
        CALL WRITE_POLYNOMIAL( &
             REGION(ISURF)%PTR%POLY(1)%PTR, &
             'dev_intersection_surface_surface/region_' // str1 // '_x.bern' )
        CALL WRITE_POLYNOMIAL( &
             REGION(ISURF)%PTR%POLY(2)%PTR, &
             'dev_intersection_surface_surface/region_' // str1 // '_pn.bern' )
        IF ( ASSOCIATED(REGION(ISURF)%PTR%CHILD) ) THEN
           PRINT *,SIZE(REGION(ISURF)%PTR%CHILD),' CHILDREN'
           DO ICHILD = 1,SIZE(REGION(ISURF)%PTR%CHILD)
              WRITE (STR2,'(I1)') ICHILD
              CALL WRITE_POLYNOMIAL( &
                   REGION(ISURF)%PTR%CHILD(ICHILD)%POLY(1)%PTR, &
                   'dev_intersection_surface_surface/child_' // str1 // '_' // str2 // '_x.bern' )
              CALL WRITE_POLYNOMIAL( &
                   REGION(ISURF)%PTR%CHILD(ICHILD)%POLY(2)%PTR, &
                   'dev_intersection_surface_surface/child_' // str1 // '_' // str2 // '_pn.bern' )
           END DO
        END IF
     END DO
     STOP
  END IF

  ! if the loop detection criterion IS satisfied, there can be no singular
  ! intersection point, or closed intersection curve. 
  ! Run the algorithm for finding the potential simple intersection branches.
  if ( stat_loopdetection == 0 ) then ! <---------------------------------+
     IF ( .false. ) THEN
        PRINT *,'INTERSECT SIMPLE SURFACES, P =',PARAM_VECTOR
     ELSE
        !IF ( INTERDAT%NC > 0 ) STOP 'DEBUG PARAM_VECTOR'
        IF (.FALSE.) THEN
           PRINT *,'INTERSECT SIMPLE SURFACES, P =',PARAM_VECTOR
           PRINT *,'UVBOXES ='
           DO ISURF = 1,2
              PRINT *,REGION(ISURF)%PTR%UVBOX
           END DO
           CALL WRITE_POLYNOMIAL( REGION(1)%PTR%POLY(1)%PTR, 'trace_intersection_polyline/x1.bern' )
           CALL WRITE_POLYNOMIAL( REGION(2)%PTR%POLY(1)%PTR, 'trace_intersection_polyline/x2.bern' )
           CALL WRITE_POLYNOMIAL( REGION(1)%PTR%POLY(2)%PTR, 'trace_intersection_polyline/pn1.bern' )
           CALL WRITE_POLYNOMIAL( REGION(2)%PTR%POLY(2)%PTR, 'trace_intersection_polyline/pn2.bern' )
        END IF

        do isurf = 1,2
           allocate( newregion(isurf)%ptr )
           call copy_region( &
                region(isurf)%ptr, &
                newregion(isurf)%ptr )
        end do

        nuvxyz = 0
        call intersect_simple_surfaces( &                                    !
             surfroot, &                                                     !
             newregion, &                                                    !
             param_vector, &                                                 !
             interdat, &                                                     !
             uvxyz, &                                                        !
             nuvxyz, &                                                       !
             stat_degeneracy )                                               !
        if ( allocated(uvxyz) ) deallocate( uvxyz )                          !

        do isurf = 1,2
           if ( newregion(isurf)%ptr%npts > 0 ) then
              call append_n( &
                   region(isurf)%ptr%ipts, &
                   region(isurf)%ptr%npts, &
                   newregion(isurf)%ptr%ipts(1:newregion(isurf)%ptr%npts), &
                   2, &
                   unique=.true. )
           end if
           nullify( &
                newregion(isurf)%ptr%xyzbox,      &
                newregion(isurf)%ptr%poly(1)%ptr, &
                newregion(isurf)%ptr%poly(2)%ptr )
           deallocate( newregion(isurf)%ptr%poly )
           call free_region_tree( newregion(isurf)%ptr )
           deallocate( newregion(isurf)%ptr )
        end do
        !PRINT *,'';PRINT *,'';PRINT *,''
     END IF
     return                                                               !
  end if ! <--------------------------------------------------------------+


  if ( stat_collineal <= 0 ) then ! <-------------------------------------+
     ! the pair of collineal points discovered is a pair of corners, so   !
     ! we subividide each region at its center point                      !
     do isurf = 1,2 ! <------------------------------------------+        !
        uv_subdiv(:,isurf) = 0.5_fp * ( &                        !        !
             region(isurf)%ptr%uvbox([1,3]) + &                  !        !
             region(isurf)%ptr%uvbox([2,4]) )                    !        !
     end do ! <--------------------------------------------------+        !
  else ! -----------------------------------------------------------------+
     ! ... if no pair of collineal corners, run Newton-Raphson algorithm. !
     ! Define feasible region for the Newton algorithm ( uvboxes of the   !
     ! surface regions )                                                  !
     do isurf = 1,2 ! <------------------------------------------+        !
        uv_collineal(:,isurf) = 0.5_fp * ( &                     !        !
             region(isurf)%ptr%uvbox([1,3]) + &                  !        !
             region(isurf)%ptr%uvbox([2,4]) )                    !        !
        lowerb(2*isurf-[1,0]) = region(isurf)%ptr%uvbox([1,3])   !        !
        upperb(2*isurf-[1,0]) = region(isurf)%ptr%uvbox([2,4])   !        !
     end do ! <--------------------------------------------------+        !
     ! expand the feasible region by a factor ( 1 + 2*EPSuv )             !
     ranges = 0.5_fp * (upperb - lowerb)                                  !
     lowerb = lowerb - EPSuv * ranges                                     !
     upperb = upperb + EPSuv * ranges                                     !
     !                                                                    !
     ! run box-constrained Newton-Raphson algorithm                       !
     call find_collineal_points( &                                        !
          surfroot, &                                                     !
          uv_collineal, &                                                 !
          normal_collineal, &                                             !
          lowerb, &                                                       !
          upperb, &                                                       !
          stat_collineal )                                                !

     if ( stat_collineal <= 0 ) uv_subdiv = uv_collineal
     IF ( DEBUG ) THEN
        PRINT *,'COLLINEAL INTERIOR POINTS?', (STAT_COLLINEAL == 0)
        IF ( STAT_COLLINEAL == 0 ) PRINT *,UV_COLLINEAL
     END IF
  end if ! <--------------------------------------------------------------+

  !PRINT *, 'STAT_COLLINEAL = ',STAT_COLLINEAL
  !IF ( STAT_COLLINEAL <= 0 ) PRINT *,'COLLINEAL POINT FOUND AT UV=',REAL(UV_COLLINEAL)

  if ( stat_collineal > 0 ) then ! <-------------------------------------------------------------------------+
     ! if no pair of collineal points has been found, subdivide the surface                                  !
     ! with largest Gauss map (define "largest" : approx. solid angle, area                                  !
     ! of spherical polygon...) , or both surfaces if each one's Gauss map                                   !
     ! does not fit in one single hemisphere) and carry on the recursion on                                  !
     ! the new pairs of surface regions                                                                      !
     ! ( with no additions to the 'intersection_data' structure )                                            !
     do isurf = 1,2 ! <------------------------------------------+                                           !
        uv_subdiv(:,isurf) = 0.5_fp * ( &                     !                                           !
             region(isurf)%ptr%uvbox([1,3]) + &                  !                                           !
             region(isurf)%ptr%uvbox([2,4]) )                    !                                           !
     end do ! <--------------------------------------------------+                                           !
     !                                                                                                       !
  elseif ( stat_collineal < 0 ) then ! ----------------------------------------------------------------------+
     ! else, if the pair of collineal points are coincident, we have a singular                              !
     ! point (tangential contact). We need to characterize the nature of that point                          !
     ! - isolated contact point ( indicates a potential non-manifoldness )                                   !
     ! - point on a tangential intersection curve                                                            !
     ! - branch point ( two transversal intersection branches meet at that point )                           !
     ! - higher-order contact point ( need to investigate higher-order derivatives )                         !
     ! Once characterized, the intersection point in appended to the 'intersection_data' structure.          !
     ! Proper action needs to be taken, depending on the nature of the point we just found:                  !
     !  - BRANCH POINT => subdivide both surfaces at that point and use Hohmeyer's method for determining    !
     !                    wether the children's Gauss maps can intersect at other points (p. 75)             !
     !  - ISOLATED CONTACT POINT: - if interior => error (non-manifoldness)                                  !
     !                            - if boundary => ignore (safe?)                                            !
     !  - TANGENTIAL INTERSECTION CURVE => error (should have been detected in the intial set of surfaces)   !
     !  - HIGHER ORDER CONTACT POINT: ???                                                                    !
     !                                                                                                       !
     !                                                                                                       !
     ! compute tangent and normal vectors, required for singular point characterization                      !
     do isurf = 1,2 ! <-------------------------------------+                                                !
        do ivar = 1,2 ! <------------------+                !                                                !
           call evald1( &                  !                !                                                !
                dxyz_duv(:,ivar,isurf), &  !                !                                                !
                surfroot(isurf)%ptr, &     !                !                                                !
                uv_collineal(:,isurf), &   !                !                                                !
                ivar )                     !                !                                                !
        end do ! <-------------------------+                !                                                !
        n(:,isurf) = cross( &                               !                                                !
             dxyz_duv(:,1,isurf), dxyz_duv(:,2,isurf) )     !                                                !
     end do ! <---------------------------------------------+                                                !
     call characterize_tangential_intersection_point( &                                                      !
          surfroot, &                                                                                        !
          uv_collineal, &                                                                                    !
          dxyz_duv, &                                                                                        !
          n, &                                                                                               !
          dxyz_ds, &                                                                                         !
          stat_singularpoint )                                                                               !
     !                                                                                                       !
     call eval( &
          xyz_collineal, &
          surfroot(1)%ptr, &
          uv_collineal(:,1) )
     call add_intersection_point( &
          uv_collineal, &
          xyz_collineal, &
          surfroot, &
          interdat, &
          idpt )
     SELECT CASE (stat_singularpoint)                                                                        !
     CASE (0)                                                                                                !
        PRINT *,'>>> tangential intersection curve <<<'                                                      !
     CASE (1)                                                                                                !
        PRINT *,'>>> branch point <<<'                                                                       !
        stat_loopdetection = 3
     CASE (2)                                                                                                !
        PRINT *,'>>> isolated tangential contact point <<<'                                                  !
     CASE (3)                                                                                                !
        PRINT *,'>>> higher-order contact point <<<'                                                         !
     END SELECT                                                                                              !
     !                                                                                                       !
  else ! ----------------------------------------------------------------------------------------------------+
     ! (stat_collineal = 0)                                                                                  !
     !stat_loopdetection = 0
     !do isurf = 1,2
     !   if ( maxval(region(isurf)%ptr%poly(2)%ptr%degr) > 0 ) then
     !      stat_loopdetection = stat_loopdetection + isurf
     !   end if
     !end do
     uv_subdiv = uv_collineal
     if ( stat_loopdetection /= 3 ) then
        isurf = stat_loopdetection
        do ivar = 1,2
           isonboundary(ivar) = ( &
                abs( uv_collineal(ivar,isurf) - region(isurf)%ptr%uvbox(2*ivar-1) ) < EPSregion .or. &
                abs( uv_collineal(ivar,isurf) - region(isurf)%ptr%uvbox(2*ivar) )   < EPSregion )
        end do
        if ( all(isonboundary) ) then
           uv_subdiv(:,isurf) = 0.5_fp * ( region(isurf)%ptr%uvbox([1,3]) + region(isurf)%ptr%uvbox([2,4]) )
        end if
     end if
  end if ! <-------------------------------------------------------------------------------------------------+


  ! subdivide both surfaces/gauss maps ...
  !PRINT *,'SUBDIVIDE AT UV =',UV_COLLINEAL
  nchild(:) = 1
  do isurf = 1,2 ! <-------------------------------------------------------------------------------------------------+
     if ( stat_loopdetection == isurf .or. stat_loopdetection == 3 ) then ! <-------------------------------------+  !
        call subdiv_region( &                                                                                     !  !
             region(isurf)%ptr, &                                                                                 !  !
             uv_subdiv(:,isurf), &                                                                             !  !
             stat_subdiv )                                                                                        !  !
        !PRINT *,'SURF #',ISURF,', STAT_SUBDIV =',STAT_SUBDIV
        !                                                                                                         !  !
        if ( stat_subdiv == 2 ) then ! <-----------------------------------------------------------------------+  !  !
           ! the subdivision point is at one of the surface's corners                                          !  !  !
           nchild(isurf) = 1                                                                                   !  !  !
        else ! ------------------------------------------------------------------------------------------------+  !  !
           ! the subdivision point is interior to the surface                                                  !  !  !
           nchild(isurf) = size(region(isurf)%ptr%child)                                                       !  !  !
           do ivar = 1,2 ! <-------------------------------------+                                             !  !  !
              uv_subdiv(ivar,isurf) = 0.5_fp * ( 1._fp + &    !                                             !  !  !
                   ab2n1p1( &                                    !                                             !  !  !
                   uv_subdiv(ivar,isurf), &                   !                                             !  !  !
                   region(isurf)%ptr%uvbox(2*ivar-1), &          !                                             !  !  !
                   region(isurf)%ptr%uvbox(2*ivar) ) )           !                                             !  !  !
           end do ! <--------------------------------------------+                                             !  !  !
           if ( stat_subdiv == 0 ) then ! <-----------------------------------------------------------------+  !  !  !
              ! 4 children                                                                                  !  !  !  !
              allocate( &                                                                                   !  !  !  !
                   region(isurf)%ptr%child(1)%poly(2), region(isurf)%ptr%child(2)%poly(2), &                !  !  !  !
                   region(isurf)%ptr%child(3)%poly(2), region(isurf)%ptr%child(4)%poly(2) )                 !  !  !  !
              ! subdivide surface                                                                           !  !  !  !
              allocate( &                                                                                   !  !  !  !
                   region(isurf)%ptr%child(1)%poly(1)%ptr, region(isurf)%ptr%child(2)%poly(1)%ptr, &        !  !  !  !
                   region(isurf)%ptr%child(3)%poly(1)%ptr, region(isurf)%ptr%child(4)%poly(1)%ptr )         !  !  !  !
              call subdiv_bezier2( &                                                                        !  !  !  !
                   region(isurf)%ptr%poly(1)%ptr, &                                                         !  !  !  !
                   uv_subdiv(:,isurf), &                                                                 !  !  !  !
                   bsw=region(isurf)%ptr%child(1)%poly(1)%ptr, &                                            !  !  !  !
                   bse=region(isurf)%ptr%child(2)%poly(1)%ptr, &                                            !  !  !  !
                   bnw=region(isurf)%ptr%child(3)%poly(1)%ptr, &                                            !  !  !  !
                   bne=region(isurf)%ptr%child(4)%poly(1)%ptr )                                             !  !  !  !
              ! subdivide pseudo-normal patch                                                               !  !  !  !
              allocate( &                                                                                   !  !  !  !
                   region(isurf)%ptr%child(1)%poly(2)%ptr, region(isurf)%ptr%child(2)%poly(2)%ptr, &        !  !  !  !
                   region(isurf)%ptr%child(3)%poly(2)%ptr, region(isurf)%ptr%child(4)%poly(2)%ptr )         !  !  !  !
              call subdiv_bezier2( &                                                                        !  !  !  !
                   region(isurf)%ptr%poly(2)%ptr, &                                                         !  !  !  !
                   uv_subdiv(:,isurf), &                                                                 !  !  !  !
                   bsw=region(isurf)%ptr%child(1)%poly(2)%ptr, &                                            !  !  !  !
                   bse=region(isurf)%ptr%child(2)%poly(2)%ptr, &                                            !  !  !  !
                   bnw=region(isurf)%ptr%child(3)%poly(2)%ptr, &                                            !  !  !  !
                   bne=region(isurf)%ptr%child(4)%poly(2)%ptr )                                             !  !  !  !
           elseif ( stat_subdiv == 1 ) then ! --------------------------------------------------------------+  !  !  !
              ! 2 children                                                                                  !  !  !  !
              allocate( region(isurf)%ptr%child(1)%poly(2), region(isurf)%ptr%child(2)%poly(2) )            !  !  !  !
              allocate( region(isurf)%ptr%child(1)%poly(1)%ptr, region(isurf)%ptr%child(2)%poly(1)%ptr )    !  !  !  !
              allocate( region(isurf)%ptr%child(1)%poly(2)%ptr, region(isurf)%ptr%child(2)%poly(2)%ptr )    !  !  !  !
              !                                                                                             !  !  !  !
              if ( region(isurf)%ptr%child(2)%uvbox(1) <= &                                                 !  !  !  !
                   region(isurf)%ptr%uvbox(1) + EPSregion ) then ! <-------------------------+              !  !  !  !
                 call subdiv_bezier2_only_v( &                                               !              !  !  !  !
                      region(isurf)%ptr%poly(1)%ptr, &                                       !              !  !  !  !
                      v=uv_subdiv(2,isurf), &                                             !              !  !  !  !
                      bs=region(isurf)%ptr%child(1)%poly(1)%ptr, &                           !              !  !  !  !
                      bn=region(isurf)%ptr%child(2)%poly(1)%ptr )                            !              !  !  !  !
                 call subdiv_bezier2_only_v( &                                               !              !  !  !  !
                      region(isurf)%ptr%poly(2)%ptr, &                                       !              !  !  !  !
                      v=uv_subdiv(2,isurf), &                                             !              !  !  !  !
                      bs=region(isurf)%ptr%child(1)%poly(2)%ptr, &                           !              !  !  !  !
                      bn=region(isurf)%ptr%child(2)%poly(2)%ptr )                            !              !  !  !  !
              else ! ------------------------------------------------------------------------+              !  !  !  !
                 call subdiv_bezier2_only_u( &                                               !              !  !  !  !
                      region(isurf)%ptr%poly(1)%ptr, &                                       !              !  !  !  !
                      u=uv_subdiv(1,isurf), &                                             !              !  !  !  !
                      bw=region(isurf)%ptr%child(1)%poly(1)%ptr, &                           !              !  !  !  !
                      be=region(isurf)%ptr%child(2)%poly(1)%ptr )                            !              !  !  !  !
                 call subdiv_bezier2_only_u( &                                               !              !  !  !  !
                      region(isurf)%ptr%poly(2)%ptr, &                                       !              !  !  !  !
                      u=uv_subdiv(1,isurf), &                                             !              !  !  !  !
                      bw=region(isurf)%ptr%child(1)%poly(2)%ptr, &                           !              !  !  !  !
                      be=region(isurf)%ptr%child(2)%poly(2)%ptr )                            !              !  !  !  !
              end if ! <---------------------------------------------------------------------+              !  !  !  !
              !                                                                                             !  !  !  !
           end if ! <---------------------------------------------------------------------------------------+  !  !  !
           !                                                                                                   !  !  !
        end if ! <---------------------------------------------------------------------------------------------+  !  !
     end if ! <---------------------------------------------------------------------------------------------------+  !
  end do ! <---------------------------------------------------------------------------------------------------------+

  if ( all( nchild < 2 ) ) then
     PRINT *,'intersect_surface_surface : NO MORE SUBDIVISIONS...'
     stat_degeneracy = 111

     DO ISURF = 1,2
        PRINT *,'SURF #',ISURF
        PRINT *,'UVBOX =',REGION(ISURF)%PTR%UVBOX
        WRITE (STR1,'(I1)') ISURF
        CALL WRITE_POLYNOMIAL( &
             REGION(ISURF)%PTR%POLY(1)%PTR, &
             'dev_intersection_surface_surface/region_' // str1 // '_x.bern' )
        CALL WRITE_POLYNOMIAL( &
             REGION(ISURF)%PTR%POLY(2)%PTR, &
             'dev_intersection_surface_surface/region_' // str1 // '_pn.bern' )
        IF ( ASSOCIATED(REGION(ISURF)%PTR%CHILD) ) THEN
           PRINT *,SIZE(REGION(ISURF)%PTR%CHILD),' CHILDREN'
           DO ICHILD = 1,SIZE(REGION(ISURF)%PTR%CHILD)
              WRITE (STR2,'(I1)') ICHILD
              CALL WRITE_POLYNOMIAL( &
                   REGION(ISURF)%PTR%CHILD(ICHILD)%POLY(1)%PTR, &
                   'dev_intersection_surface_surface/child_' // str1 // '_' // str2 // '_x.bern' )
              CALL WRITE_POLYNOMIAL( &
                   REGION(ISURF)%PTR%CHILD(ICHILD)%POLY(2)%PTR, &
                   'dev_intersection_surface_surface/child_' // str1 // '_' // str2 // '_pn.bern' )
           END DO
        END IF
     END DO

     RETURN
  end if






  ! *** DEBUG *******************************
  IF ( .false. ) THEN
     !IF ( ASSOCIATED(REGION(1)%PTR%PARENT) .AND. ASSOCIATED(REGION(2)%PTR%PARENT) ) THEN
     DO ISURF = 1,2
        PRINT *,'SURF #',ISURF
        PRINT *,'UVBOX =',REGION(ISURF)%PTR%UVBOX
        WRITE (STR1,'(I1)') ISURF
        CALL WRITE_POLYNOMIAL( &
             REGION(ISURF)%PTR%POLY(1)%PTR, &
             'dev_intersection_surface_surface/region_' // str1 // '_x.bern' )
        CALL WRITE_POLYNOMIAL( &
             REGION(ISURF)%PTR%POLY(2)%PTR, &
             'dev_intersection_surface_surface/region_' // str1 // '_pn.bern' )
        IF ( ASSOCIATED(REGION(ISURF)%PTR%CHILD) ) THEN
           PRINT *,SIZE(REGION(ISURF)%PTR%CHILD),' CHILDREN'
           DO ICHILD = 1,SIZE(REGION(ISURF)%PTR%CHILD)
              WRITE (STR2,'(I1)') ICHILD
              CALL WRITE_POLYNOMIAL( &
                   REGION(ISURF)%PTR%CHILD(ICHILD)%POLY(1)%PTR, &
                   'dev_intersection_surface_surface/child_' // str1 // '_' // str2 // '_x.bern' )
              CALL WRITE_POLYNOMIAL( &
                   REGION(ISURF)%PTR%CHILD(ICHILD)%POLY(2)%PTR, &
                   'dev_intersection_surface_surface/child_' // str1 // '_' // str2 // '_pn.bern' )
           END DO
        END IF
     END DO
     STOP '---> DEBUG SUBDIV'
  END IF
  ! *****************************************





  ! ...and carry on the recursion on the new pairs of regions
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
        call intersect_surface_surface( &                      !   !
             surfroot, &                                       !   !
             newregion, &                                      !   !
             interdat, &                                       !   !
             stat_degeneracy )                                 !   !
        !                                                      !   !
     end do ! <------------------------------------------------+   !
  end do ! <-------------------------------------------------------+

end subroutine intersect_surface_surface
