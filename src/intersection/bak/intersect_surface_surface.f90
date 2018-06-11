recursive subroutine intersect_surface_surface( &
     surfroot, &
     surfc, &
     surfpn, &
     region, &
     !stat_collineal_parents, &
     !uv_collineal_parents, &
     interdat ) 
  use mod_math
  use mod_chebyshev
  use mod_diffgeom
  use mod_types_intersection
  use mod_tolerances
  ! ...
  implicit none
  type(ptr_parametric_surface),  intent(in)    :: surfroot(2)
  type(ptr_chebyshev_series2),   intent(in)    :: surfc(2)
  type(ptr_chebyshev_series2),   intent(in)    :: surfpn(2)
  !integer,                       intent(in)    :: stat_collineal_parents
  !real(kind=MATHpr),             intent(in)    :: uv_collineal_parents(2,2)
  type(ptr_surface_region),      intent(inout) :: region(2)
  type(type_intersection_data),  intent(inout) :: interdat
  logical                                      :: overlap
  real(kind=MATHpr)                            :: uv_collineal(2,2)
  real(kind=MATHpr), dimension(4)              :: lowerb, upperb, ranges
  real(kind=MATHpr)                            :: param_vector(3)
  integer                                      :: stat_loopdetection
  integer                                      :: stat_collineal
  integer                                      :: stat_singularpoint
  real(kind=MATHpr)                            :: s1(3,2,2), n(3,2), xyz_s(3,2)
  integer                                      :: stat_subdiv
  integer                                      :: nchild(2), nchildtotal
  type(ptr_surface_region)                     :: newregion(2)
  type(type_chebyshev_series2), allocatable    :: childsurfc(:), childsurfpn(:)
  type(ptr_chebyshev_series2)                  :: newsurfc(2), newsurfpn(2)
  !type(type_chebyshev_series2), allocatable    :: newsurfc(:), newsurfpn(:)
  !type(ptr_chebyshev_series2), allocatable     :: newsurfc(:), newsurfpn(:)
  integer                                      :: isurf, ichild, jchild, ivar

  CHARACTER                                    :: STR1, STR2

  ! ...

  ! compute OBBs
  do isurf = 1,2
     if ( .not.associated(region(isurf)%ptr%xyzbox) ) then
        allocate( region(isurf)%ptr%xyzbox )
        !PRINT *,'ASSOCIATED?',ASSOCIATED(SURFC(ISURF)%PTR)
        !IF ( ASSOCIATED(SURFC(ISURF)%PTR) ) THEN
        !   PRINT *,SURFC(ISURF)%PTR
        !   PRINT *,'   COEF ALLOCATED?', ALLOCATED( SURFC(ISURF)%PTR%COEF )
        !   PRINT *,'SIZE=',SIZE( SURFC(ISURF)%PTR%COEF )
        !END IF
        call chebOBB2( &
             surfc(isurf)%ptr%coef(1:surfc(isurf)%ptr%degr(1)+1, 1:surfc(isurf)%ptr%degr(2)+1,:), &
             surfc(isurf)%ptr%degr, &
             region(isurf)%ptr%xyzbox )
     end if
  end do

  ! test OBBs for overlap
  call overlap_OBBs( &
       [ region(1)%ptr%xyzbox , region(2)%ptr%xyzbox ], &
       overlap )

  if ( .not.overlap ) return ! disjoint OBBs means empty intersection


  ! use Hohmeyer's loop detection technique
  !if ( stat_collineal_parents > 0 ) then
  do isurf = 1,2
     if ( .not.associated(region(isurf)%ptr%pnbox) ) then
        allocate( region(isurf)%ptr%pnbox )
        call chebOBB2( &
             surfpn(isurf)%ptr%coef(1:surfpn(isurf)%ptr%degr(1)+1, 1:surfpn(isurf)%ptr%degr(2)+1,:), &
             surfpn(isurf)%ptr%degr, &
             region(isurf)%ptr%pnbox )
     end if
  end do
  call hohmeyer_loop_detection( &
       [ region(1)%ptr%pnbox , region(2)%ptr%pnbox ], &
       param_vector, &
       stat_loopdetection ) 
  !else
  !   ! the Gauss maps of the two regions are known to intersect at the point 'uv_collineal_parents',
  !   ! check if they can intersect at other points
     

  !end if ! ( stat_collineal_parents > 0 )

  !PRINT *,'stat_loopdetection =', stat_loopdetection

  ! if the loop detection criterion IS satisfied, there can be no singular
  ! intersection point, or closed intersection curve. 
  ! Run the algorithm for finding the potential simple intersection branches.
  if ( stat_loopdetection == 0 ) then
     PRINT *,'PARAM_VECTOR =',REAL(PARAM_VECTOR)
     call intersect_simple_surfaces( &
          surfroot, &
          surfc, &
          region, &
          param_vector, &
          interdat )
     ! (...)
     return
  end if



  ! if the loop detection criterion IS NOT satisfied, look for a potential 
  ! pair of collineal points ...
  !  ... check surface region corners first ...
  call find_collineal_corners( &
       surfroot, &
       uv_collineal, &
       stat_collineal )

  if ( stat_collineal > 0 ) then
     ! ... if no pair of collineal corners, run Newton-Raphson algorithm.
     ! Define feasible region for the Newton algorithm ( uvboxes of the 
     ! surface regions )
     do isurf = 1,2
        uv_collineal(:,isurf) = 0.5_MATHpr * sum( region(isurf)%ptr%uvbox, 1 )
        lowerb(2*isurf-[1,0]) = region(isurf)%ptr%uvbox(1,:)
        upperb(2*isurf-[1,0]) = region(isurf)%ptr%uvbox(2,:)
     end do
     ! expand the feasible region by a factor ( 1 + 2*EPSuv )
     ranges = 0.5_MATHpr * (upperb - lowerb)
     lowerb = lowerb - EPSuv * ranges
     upperb = upperb + EPSuv * ranges

     ! run box-constrained Newton-Raphson algorithm
     call find_collineal_points( &
          surfroot, &
          uv_collineal, &
          lowerb, &
          upperb, &
          stat_collineal )

  end if



  if ( stat_collineal > 0 ) then
     ! if no pair of collineal points has been found, subdivide the surface
     ! with largest Gauss map (define "largest" : approx. solid angle, area 
     ! of spherical polygon...) , or both surfaces if each one's Gauss map 
     ! does not fit in one single hemisphere) and carry on the recursion on 
     ! the new pairs of surface regions
     ! ( with no additions to the 'intersection_data' structure )
     do isurf = 1,2
        uv_collineal(:,isurf) = 0.5_MATHpr * sum( region(isurf)%ptr%uvbox, 1 )
     end do

  elseif ( stat_collineal < 0 ) then
     
     ! else, if the pair of collineal points are coincident, we have a singular
     ! point (tangential contact). We need to characterize the nature of that point
     ! - isolated contact point ( indicates a potential non-manifoldness )
     ! - point on a tangential intersection curve
     ! - branch point ( two transversal intersection branches meet at that point )
     ! - higher-order contact point ( need to investigate higher-order derivatives )
     ! Once characterized, the intersection point in appended to the 'intersection_data' structure.
     ! Proper action needs to be taken, depending on the nature of the point we just found:
     !  - BRANCH POINT => subdivide both surfaces at that point and use Hohmeyer's method for determining
     !                    wether the children's Gauss maps can intersect at other points (p. 75)
     !  - ISOLATED CONTACT POINT: - if interior => error (non-manifoldness)
     !                            - if boundary => ignore (safe?)
     !  - TANGENTIAL INTERSECTION CURVE => error (should have been detected in the intial set of surfaces)
     !  - HIGHER ORDER CONTACT POINT: ???
     stat_loopdetection = 0
     do isurf = 1,2
        if ( any( surfpn(isurf)%ptr%degr > 0 ) ) stat_loopdetection = stat_loopdetection + isurf
     end do

     ! compute tangent and normal vectors, required for singular point characterization
     do isurf = 1,2
        do ivar = 1,2
           call evald( s1(:,ivar,isurf), surfroot(isurf)%ptr, uv_collineal(:,isurf), ivar )
        end do
        n(:,isurf) = cross( s1(:,1,isurf), s1(:,2,isurf) )
     end do
     call characterize_tangential_intersection_point( &
          surfroot, &
          uv_collineal, &
          xyz_s, &
          stat_singularpoint, &
          s1, &
          n )

     SELECT CASE (stat_singularpoint)
     CASE (0)
        PRINT *,'tangential intersection curve'
     CASE (1)
        PRINT *,'branch point'
     CASE (2)
        PRINT *,'isolated tangential contact point'
     CASE (3)
        PRINT *,'higher-order contact point'
     END SELECT

     !call add_intersection_point( &
     !     surfroot, &
     !     region, &
     !     uv_collineal, &
     !     interdat, &
     !     id )

  end if

  ! else (if the two collineal points are not coincident), subdivide both surfaces at that point
  ! and carry on the recursion on the new pairs of surface regions
  ! use Hohmeyer's method for determining wether the children's Gauss maps can intersect at other points
  ! ( with no additions to the 'intersection_data' structure )








  ! >>>>> block repeated three times, put in a separate subroutine or use GOTO -----------
  ! .... special cases where one or obth surfaces is not subdivided NEEDS to be implemented
  nchild(:) = 0
  do isurf = 1,2
     if ( stat_loopdetection == isurf .or. stat_loopdetection == 3 ) then
        call subdiv_surface_region( &
             region(isurf)%ptr, &
             uv_collineal(:,isurf), &
             stat_subdiv )
        if ( stat_subdiv < 3 ) then
           nchild(isurf) = size(region(isurf)%ptr%child)
        else
           nchild(isurf) = 1
        end if
     else
        nchild(isurf) = 1
     end if
  end do
  nchildtotal = sum( nchild )
  if ( all( nchild < 2 ) ) then
     PRINT *,'NO MORE SUBDIVISIONS...'
     return
  end if

  ! compute changes of variables for position vector and pseudo-normal patch
  allocate( newsurfc(nchildtotal), newsurfpn(nchildtotal) )
  do isurf = 1,2
     if ( nchild(isurf) == 1 ) then
        !newsurfc( (isurf-1)*nchild(1) + 1 )%ptr => surfc(isurf)%ptr
        !newsurfpn( (isurf-1)*nchild(1) + 1 )%ptr => surfpn(isurf)%ptr
        !newsurfc( (isurf-1)*nchild(1) + 1 ) = surfc(isurf)
        !newsurfpn( (isurf-1)*nchild(1) + 1 ) = surfpn(isurf)
        childsurfc( (isurf-1)*nchild(1) + 1 ) = surfc(isurf)
        childsurfpn( (isurf-1)*nchild(1) + 1 ) = surfpn(isurf)
     else
        do ichild = 1,nchild(isurf)
           !allocate( newsurfc( (isurf-1)*nchild(1) + ichild )%ptr, &
           !     newsurfpn( (isurf-1)*nchild(1) + ichild )%ptr )

           !PRINT *,' AVANT, ASSOCIATED?', ASSOCIATED(newsurfc( (isurf-1)*nchild(1) + ichild )%ptr)
           call chgvar2( &
                surfroot(isurf)%ptr%s, &                      ! surfc(isurf), & (WRONG)
                childsurfc( (isurf-1)*nchild(1) + ichild ), &
                !newsurfc( (isurf-1)*nchild(1) + ichild )%ptr, &
                region(isurf)%ptr%child(ichild)%uvbox(1,:), &
                region(isurf)%ptr%child(ichild)%uvbox(2,:) )
           !PRINT *,' APRES, ASSOCIATED?', ASSOCIATED(newsurfc( (isurf-1)*nchild(1) + ichild )%ptr)

           call chgvar2( &
                surfroot(isurf)%ptr%pn, &                     ! surfpn(isurf), & (WRONG) 
                childsurfpn( (isurf-1)*nchild(1) + ichild ), &
                !newsurfpn( (isurf-1)*nchild(1) + ichild )%ptr, &
                region(isurf)%ptr%child(ichild)%uvbox(1,:), &
                region(isurf)%ptr%child(ichild)%uvbox(2,:) )
        end do
     end if
  end do


  ! *** DEBUG *******************************
  IF ( .FALSE. ) THEN
     IF ( ASSOCIATED(REGION(1)%PTR%PARENT) .OR. ASSOCIATED(REGION(2)%PTR%PARENT) ) THEN
        DO ISURF = 1,2
           WRITE (STR1,'(I1)') ISURF
           CALL WRITE_CHEBYSHEV_SERIES2( &
                SURFROOT(ISURF)%PTR%S, &
                'pre_intersection/debug_subdivision/surf_' // STR1 // '_root.cheb' )

           CALL WRITE_CHEBYSHEV_SERIES2( &
                SURFC( ISURF )%ptr, &
                'pre_intersection/debug_subdivision/surf_' // STR1 // '_parent.cheb' )

           DO ICHILD = 1,NCHILD(ISURF)
              WRITE (STR2,'(I1)') ICHILD
              CALL WRITE_CHEBYSHEV_SERIES2( &
                   !NEWSURFC( (ISURF-1)*NCHILD(1) + ICHILD )%ptr, &
                   childSURFC( (ISURF-1)*NCHILD(1) + ICHILD ), &
                   'pre_intersection/debug_subdivision/surf_' // STR1 // '_child_' // STR2 // '.cheb' )
           END DO

           CALL EXPORT_SURFACE_REGION_TREE( &
                REGION(ISURF)%PTR, &
                'pre_intersection/debug_subdivision/surf_' // STR1 // '_tree.dat' )
        END DO
        STOP
     END IF
  END IF
  ! *****************************************


  ! call intersect_surface_surface on pairs of child surface region
  do jchild = 1,nchild(2)

     if ( nchild(2) == 1 ) then
        newregion(2)%ptr => region(2)%ptr
        newsurfc(2)%ptr => surfc(2)%ptr
     else
        newregion(2)%ptr => region(2)%ptr%child(jchild)
        newsurfc(2)%ptr => childsurfc(jchild,2)%ptr
     end if

     do ichild = 1,nchild(1)

        if ( nchild(1) == 1 ) then
           newregion(1)%ptr => region(1)%ptr
        else
           newregion(1)%ptr => region(1)%ptr%child(ichild)
        end if

        call intersect_surface_surface( &
             surfroot, &
             ![ newsurfc( ichild ), newsurfc( nchild(1) + jchild ) ], &
             newsurfc( [ichild, nchild(1) + jchild] ), &
             ![ newsurfpn( ichild ), newsurfpn( nchild(1) + jchild ) ], &
             newsurfpn( [ichild, nchild(1) + jchild] ), &
             newregion, &
             !stat_collineal, &
             !uv_collineal, &
             interdat )

     end do
  end do

  deallocate( newsurfc, newsurfpn )

  ! -------------------------------------------------------------------------------- <<<<<





end subroutine intersect_surface_surface
