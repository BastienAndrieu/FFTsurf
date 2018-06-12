recursive subroutine intersect_surface_surface_bezier( &
     surfroot, &
     bs, &
     bpn, &
     region, &
     interdat ) 
  use mod_math
  use mod_bernstein
  use mod_diffgeom
  use mod_types_intersection
  use mod_tolerances
  ! ...
  type(ptr_parametric_surface),  intent(in)    :: surfroot(2)
  type(ptr_bernstein_series2),   intent(in)    :: bs(2)
  type(ptr_bernstein_series2),   intent(in)    :: bpn(2)
  type(ptr_surface_region),      intent(inout) :: region(2)
  type(type_intersection_data),  intent(inout) :: interdat
  logical                                      :: overlap
  integer                                      :: nbpn(2)
  real(kind=MATHpr)                            :: uv_collineal(2,2)
  real(kind=MATHpr), dimension(4)              :: lowerb, upperb, ranges
  real(kind=MATHpr)                            :: param_vector(3)
  integer                                      :: stat_loopdetection
  integer                                      :: stat_collineal
  integer                                      :: stat_singularpoint
  real(kind=MATHpr)                            :: s1(3,2,2), n(3,2), xyz_s(3,2)
  integer                                      :: stat_subdiv
  integer                                      :: nchild(2)
  type(ptr_surface_region)                     :: newregion(2)
  real(kind=MATHpr)                            :: st_subdiv(2)
  type(type_bernstein_series2), dimension(4,2), target :: childbs, childbpn
  type(ptr_bernstein_series2), dimension(2)    :: newbs, newbpn
  integer                                      :: isurf, ichild
  CHARACTER                                    :: STR1, STR2
  INTEGER                                      :: COUNTER
  SAVE COUNTER
  DATA COUNTER /0/
  
  COUNTER = COUNTER + 1

  ! compute OBBs
  do isurf = 1,2
     if ( .not.associated(region(isurf)%ptr%xyzbox) ) then
        allocate( region(isurf)%ptr%xyzbox )
        call bernOBB2( &
             bs(isurf)%ptr%coef(1:bs(isurf)%ptr%degr(1)+1, 1:bs(isurf)%ptr%degr(2)+1,:), &
             bs(isurf)%ptr%degr, &
             region(isurf)%ptr%xyzbox )
     end if
  end do

  ! test OBBs for overlap
  call overlap_OBBs( &
       region(1)%ptr%xyzbox , &![ region(1)%ptr%xyzbox , region(2)%ptr%xyzbox ], &
       region(2)%ptr%xyzbox, &
       overlap )
  IF (.false.) THEN
     PRINT *,'OVERLAP ?',overlap
     DO ISURF = 1,2
        WRITE (STR1,'(I1)') ISURF
        CALL WRITE_OBB( REGION(ISURF)%PTR%XYZBOX, 'pre_intersection/xyzbox_debugobb_' // STR1 // '.dat' )
        CALL WRITE_BERNSTEIN_SERIES2( BS(ISURF)%PTR, 'pre_intersection/bs_debugobb_' // STR1 // '.bern' )
     END DO
     STOP
  END IF

  if ( .not.overlap ) return ! disjoint OBBs means empty intersection

  ! use Hohmeyer's loop detection technique
  do isurf = 1,2
     nbpn(isurf) = ( bpn(isurf)%ptr%degr(1) + 1 )*( bpn(isurf)%ptr%degr(2) + 1 )
     !PRINT *,"NBPN =",NBPN(ISURF)
  end do
  call hohmeyer_loop_detection_bezier( &
       reshape( bpn(1)%ptr%coef(1:bpn(1)%ptr%degr(1)+1,1:bpn(1)%ptr%degr(2)+1,:), [nbpn(1),3] ) , &
       nbpn(1), &
       reshape( bpn(2)%ptr%coef(1:bpn(2)%ptr%degr(1)+1,1:bpn(2)%ptr%degr(2)+1,:), [nbpn(2),3] ) , &
       nbpn(2), &
       param_vector, &
       stat_loopdetection ) 

! if the loop detection criterion IS satisfied, there can be no singular
  ! intersection point, or closed intersection curve. 
  ! Run the algorithm for finding the potential simple intersection branches.
  if ( stat_loopdetection == 0 ) then
     PRINT *,'PARAM_VECTOR =',REAL(PARAM_VECTOR)
     !call intersect_simple_surfaces( &
     !     surfroot, &
     !     surfc, &
     !     region, &
     !     param_vector, &
     !     interdat )
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

  IF ( STAT_COLLINEAL <= 0 ) PRINT *,'COLLINEAL POINT FOUND AT UV=',REAL(UV_COLLINEAL)

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
     
     !!stat_loopdetection = 0
     !!do isurf = 1,2
     !!   if ( any( surfpn(isurf)%ptr%degr > 0 ) ) stat_loopdetection = stat_loopdetection + isurf
     !!end do
     
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
        PRINT *,'>>> tangential intersection curve <<<'
     CASE (1)
        PRINT *,'>>> branch point <<<'
     CASE (2)
        PRINT *,'>>> isolated tangential contact point <<<'
     CASE (3)
        PRINT *,'>>> higher-order contact point <<<'
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



  ! SUBDIVISIONS
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




  ! compute Bezier patches for subdivisions
  do isurf = 1,2
     st_subdiv(1) = ab2n1p1( uv_collineal(1,isurf), region(isurf)%ptr%uvbox(1,1), region(isurf)%ptr%uvbox(2,1) )
     st_subdiv(2) = ab2n1p1( uv_collineal(2,isurf), region(isurf)%ptr%uvbox(1,2), region(isurf)%ptr%uvbox(2,2) )
     st_subdiv = 0.5_MATHpr * ( 1_MATHpr + st_subdiv )
     !st_subdiv = 0.5_MATHpr * ( 1_MATHpr + &
     !     ( uv_collineal(:,isurf) - region(isurf)%ptr%uvbox(1,:) ) / &
     !     ( region(isurf)%ptr%uvbox(2,:) - region(isurf)%ptr%uvbox(1,:) ) )
     if ( nchild(isurf) > 1 ) then
        !IF ( NCHILD(ISURF) < 4 ) STOP 'DEGENERATE SUBDIVISION (NOT YET IMPLEMENTED FOR BEZIER)'
        if ( nchild(isurf) == 2 ) then
           ! 2 children
           if ( abs( region(isurf)%ptr%child(1)%uvbox(2,1) - &
                region(isurf)%ptr%child(2)%uvbox(2,1) ) < EPSregion ) then
              ! subdivide along v
              call subdiv2_along_v( &
                   bs(isurf)%ptr, &
                   st_subdiv(2), &
                   childbs(1,isurf), &
                   childbs(2,isurf) )
              call subdiv2_along_v( &
                   bpn(isurf)%ptr, &
                   st_subdiv(2), &
                   childbpn(1,isurf), &
                   childbpn(2,isurf) )
           else
              ! subdivide along u
              call subdiv2_along_u( &
                   bs(isurf)%ptr, &
                   st_subdiv(1), &
                   childbs(1,isurf), &
                   childbs(2,isurf) )
              call subdiv2_along_u( &
                   bpn(isurf)%ptr, &
                   st_subdiv(1), &
                   childbpn(1,isurf), &
                   childbpn(2,isurf) )
           end if
        else
           ! 4 children
           call subdiv2( &
                bs(isurf)%ptr, &
                st_subdiv, &
                bsw=childbs(1,isurf), &
                bse=childbs(2,isurf), &
                bnw=childbs(3,isurf), &
                bne=childbs(4,isurf) )
           call subdiv2( &
                bpn(isurf)%ptr, &
                st_subdiv, &
                bsw=childbpn(1,isurf), &
                bse=childbpn(2,isurf), &
                bnw=childbpn(3,isurf), &
                bne=childbpn(4,isurf) )
        end if
     end if
  end do


  ! *** DEBUG *******************************
  IF ( .false. .AND. COUNTER > 6 ) THEN
     !IF ( ASSOCIATED(REGION(1)%PTR%PARENT) .OR. ASSOCIATED(REGION(2)%PTR%PARENT) ) THEN
        DO ISURF = 1,2
           WRITE (STR1,'(I1)') ISURF
           CALL WRITE_CHEBYSHEV_SERIES2( &
                SURFROOT(ISURF)%PTR%S, &!SURFROOT(ISURF)%PTR%PN, &
                'pre_intersection/debug_subdivision/surf_' // STR1 // '_root.cheb' )


           CALL WRITE_BERNSTEIN_SERIES2( &
                BS( ISURF )%PTR, &!BPN( ISURF )%PTR, &
                'pre_intersection/debug_subdivision/surf_' // STR1 // '_parent.bern' )

           IF ( NCHILD(ISURF) > 1 ) THEN
              DO ICHILD = 1,NCHILD(ISURF)
                 WRITE (STR2,'(I1)') ICHILD
                 CALL WRITE_BERNSTEIN_SERIES2( &
                      CHILDBS(ICHILD,ISURF), &!CHILDBPN(ICHILD,ISURF), &
                      'pre_intersection/debug_subdivision/surf_' // STR1 // '_child_' // STR2 // '.bern' )
              END DO
           END IF

           CALL EXPORT_SURFACE_REGION_TREE( &
                REGION(ISURF)%PTR, &
                'pre_intersection/debug_subdivision/surf_' // STR1 // '_tree.dat' )
        END DO
        PRINT *,'THIS MUST STOP HERE'
        STOP
     !END IF
  END IF
  ! *****************************************


  ! call intersect_surface_surface on pairs of child surface region
  do jchild = 1,nchild(2)

     if ( nchild(2) == 1 ) then
        newregion(2)%ptr => region(2)%ptr
        newbs(2)%ptr  => bs(2)%ptr
        newbpn(2)%ptr => bpn(2)%ptr
     else
        newregion(2)%ptr => region(2)%ptr%child(jchild)
        newbs(2)%ptr  => childbs(jchild,2)
        newbpn(2)%ptr => childbpn(jchild,2)
     end if

     do ichild = 1,nchild(1)

        if ( nchild(1) == 1 ) then
           newregion(1)%ptr => region(1)%ptr
           newbs(1)%ptr  => bs(1)%ptr
           newbpn(1)%ptr => bpn(1)%ptr
        else
           newregion(1)%ptr => region(1)%ptr%child(ichild)
           newbs(1)%ptr  => childbs(ichild,1)
           newbpn(1)%ptr => childbpn(ichild,1)
        end if

        call intersect_surface_surface_bezier( &
             surfroot, &
             newbs, &
             newbpn, &
             newregion, &
             interdat )

     end do
  end do



end subroutine intersect_surface_surface_bezier
