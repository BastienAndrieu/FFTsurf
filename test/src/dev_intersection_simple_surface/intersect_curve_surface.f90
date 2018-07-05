recursive subroutine intersect_curve_surface( &
     root_c, &
     root_s, &
     region_c, &
     region_s, &
     coords, &
     npts, &
     stat_newpoint )
  use mod_util    
  use mod_math
  use mod_diffgeom2
  use mod_obb
  use mod_regiontree
  use mod_tolerances
  implicit none
  type(type_curve),           intent(in)            :: root_c
  type(type_surface),         intent(in)            :: root_s
  type(type_region),          intent(inout), target :: region_c
  type(type_region),          intent(inout), target :: region_s
  real(kind=fp), allocatable, intent(inout)         :: coords(:,:) ! t,u,v,x,y,z
  integer,                    intent(inout)         :: npts
  integer,                    intent(inout)         :: stat_newpoint
  integer, allocatable                              :: sharedpts(:)
  integer                                           :: n_sharedpts
  logical                                           :: separable, overlap
  real(kind=fp)                                     :: tuv_subdiv(3)
  logical                                           :: interior(2)
  real(kind=fp)                                     :: tuv(3), xyz(3)
  integer                                           :: stat_subdiv, nchild(2)
  type(type_region), pointer                        :: newregion_c
  type(type_region), pointer                        :: newregion_s
  integer                                           :: ipt, jpt, i, j, k, ichild, jchild

  if ( stat_degeneracy > 1 ) return

  !PRINT *,''
  !PRINT *,''
  !PRINT *,' TBOX =',REGION_C%UVBOX
  !PRINT *,'UVBOX =',REGION_S%UVBOX

  ! inherit from parents all already discovered points contained inside the current regions
  if ( npts > 0 ) then
     call inherit_points( &
          region_c, &
          coords(1,1:npts), &
          npts )

     call inherit_points( &
          region_s, &
          coords(2:3,1:npts), &
          npts )
  end if

  ! get list of already discovered points common to both the current curve and surface regions
  n_sharedpts = 0
  if ( region_c%npts > 0 .and. region_s%npts > 0 ) then
     call intersection_arrays( &
          region_c%ipts(1:region_c%npts), &
          region_s%ipts(1:region_s%npts), &
          sharedpts )
     if ( allocated(sharedpts) ) n_sharedpts = size(sharedpts)
  end if
  !PRINT *,'N_SHAREDPTS =',n_sharedpts
  !IF ( N_SHAREDPTS > 0 ) PRINT *,sharedpts


  IF (.FALSE.) THEN
     if ( n_sharedpts == 0 ) then
        ! check endpoint/corner pairs for possible intersection point
        do k = 1,region_s%poly(1)%ptr%degr(2)+1,region_s%poly(1)%ptr%degr(2)
           do j = 1,region_s%poly(1)%ptr%degr(1)+1,region_s%poly(1)%ptr%degr(1)
              do i = 1,region_c%poly(1)%ptr%degr(1)+1,region_c%poly(1)%ptr%degr(1)
                 if ( sum( (region_c%poly(1)%ptr%coef(i,:,1) - region_s%poly(1)%ptr%coef(j,k,:))**2 ) < EPSxyzsqr ) then
                    xyz = 0.5_fp * ( region_c%poly(1)%ptr%coef(i,:,1) + region_s%poly(1)%ptr%coef(j,k,:) )
                    tuv(1) = real( i-1, kind=fp ) / real( region_c%poly(1)%ptr%degr(1), kind=fp )
                    tuv(2) = real( j-1, kind=fp ) / real( region_s%poly(1)%ptr%degr(1), kind=fp )
                    tuv(3) = real( k-1, kind=fp ) / real( region_s%poly(1)%ptr%degr(2), kind=fp )
                    tuv = -1._fp + 2._fp * tuv

                    call append_vector( &
                         [tuv,xyz], &
                         6, &
                         coords, &
                         npts ) 

                    call append( &
                         sharedpts, &
                         npts, &
                         newlength=n_sharedpts )

                    call add_point_bottom_up( region_c, npts )
                    call add_point_bottom_up( region_s, npts )
                 end if
              end do
           end do
        end do
     end if
  END IF


  if ( n_sharedpts == 1 ) then ! <--------------------------------------------------------------------------------------------+
     ipt = sharedpts(1)                                                                                                       !
     interior(1) = is_in_interval_strict( coords(1,ipt), region_c%uvbox(1), region_c%uvbox(2) )                               !
     interior(2) = ( &                                                                                                        !
          is_in_interval_strict( coords(2,ipt), region_s%uvbox(1), region_s%uvbox(2) ) .and. &                                !
          is_in_interval_strict( coords(3,ipt), region_s%uvbox(3), region_s%uvbox(4) ) )                                      !
     !                                                                                                                        !
     if ( any(interior) ) then ! <--------------------------------------------------------------------+                       !
        separable = .false.                                                                           !                       !
     else ! ------------------------------------------------------------------------------------------+                       !
        ! check if the curve and surface regions can intersect at other (not yet discovered) points   !                       !
        call intersect_curve_surface_elsewhere( &                                                     !                       !
             region_c%poly(1)%ptr, &                                                                         !                       !
             region_s%poly(1)%ptr, &                                                                         !                       !
             coords(4:6,sharedpts(1)), &                                                              !                       !
             separable, &                                                                             !                       !
             randomize=.true. )                                                                       !                       !
     end if ! <---------------------------------------------------------------------------------------+                       !
     !                                                                                                                        !
     if ( separable ) return ! the pair of regions of curve-surface cannot intersect except at the point already discovered   !
  end if ! <------------------------------------------------------------------------------------------------------------------+

  if ( n_sharedpts > 0 ) then ! <---------------------------------------------------------------------------------------------+
     ! are any of the already discovered points interior to the curve/surface region                                          !
     do jpt = 1,n_sharedpts ! <---------------------------------------------------------------------+                         !
        ipt = sharedpts(jpt)                                                                        !                         !
        interior(1) = is_in_interval_strict( coords(1,ipt), region_c%uvbox(1), region_c%uvbox(2) )  !                         !
        interior(2) = ( &                                                                           !                         !
             is_in_interval_strict( coords(2,ipt), region_s%uvbox(1), region_s%uvbox(2) ) .and. &   !                         !
             is_in_interval_strict( coords(3,ipt), region_s%uvbox(3), region_s%uvbox(4) ) )         !                         !
        if ( any(interior) ) then ! <-----------+                                                   !                         !
           tuv_subdiv(1) = coords(1,ipt)        !                                                   !                         !
           tuv_subdiv(2:3) = coords(2:3,ipt)    !                                                   !                         !
           exit                                 !                                                   !                         !
        end if ! <------------------------------+                                                   !                         !
     end do ! <-------------------------------------------------------------------------------------+                         !
     !                                                                                                                        !
  else ! ---------------------------------------------------------------------------------------------------------------------+
     ! n_sharedpts ==  0                                                                                                      !
     interior(:) = .false.                                                                                                    !
     !                                                                                                                        !
     ! compute bounding boxes for each region...                                                                              !
     if ( .not.associated(region_c%xyzbox) ) then ! <----------------------------------------+                                !
        allocate( region_c%xyzbox )                                                          !                                !
        call bernOBB1( &                                                                     !                                !
             region_c%poly(1)%ptr%coef(1:region_c%poly(1)%ptr%degr(1)+1,1:3,1), &                          !                                !
             region_c%poly(1)%ptr%degr(1), &                                                        !                                !
             region_c%xyzbox )                                                               !                                !
     end if ! <------------------------------------------------------------------------------+                                !
     !                                                                                                                        !
     if ( .not.associated(region_s%xyzbox) ) then ! <----------------------------------------+                                !
        allocate( region_s%xyzbox )                                                          !                                !
        call bernOBB2( &                                                                     !                                !
             region_s%poly(1)%ptr%coef(1:region_s%poly(1)%ptr%degr(1)+1,1:region_s%poly(1)%ptr%degr(2)+1,1:3), &  !                                !
             region_s%poly(1)%ptr%degr, &                                                           !                                !
             region_s%xyzbox )                                                               !                                !
     end if ! <------------------------------------------------------------------------------+                                !
     !                                                                                                                        !
     ! ... and check if these bounding boxes overlap                                                                          !
     call overlap_OBBs( &                                                                                                     !
          region_c%xyzbox, &                                                                                                  !
          region_s%xyzbox, &                                                                                                  !
          overlap )                                                                                                           !
     !                                                                                                                        !
     if ( .not.overlap ) return ! disjoint bounding boxes => empty intersection                                               !
     !                                                                                                                        !
  end if ! <------------------------------------------------------------------------------------------------------------------+


  if ( all(.not.interior) ) then ! <------------------------------------------------------------------------------------------+
     !! Search for a new intersection point using a Newton-Raphson algorithm                                                  !
     tuv = 0.5_fp * [ &                                                                                                       !
          region_c%uvbox(1) + region_c%uvbox(2), &                                                                            !
          region_s%uvbox([1,3]) + region_s%uvbox([2,4]) ]                                                                     !
     !PRINT *,'NEWTON <--- TUV0 =',TUV                                                                                        !
     !PRINT *,'------- NEWTON -------'                                                                                        !
     call newton_curve_surface( &                                                                                             !
          root_c, &                                                                                                           !
          root_s, &                                                                                                           !
          region_c%uvbox, &                                                                                                   !
          region_s%uvbox, &                                                                                                   !
          tuv, &                                                                                                              !
          stat_newpoint, &                                                                                                    !
          xyz )                                                                                                               !
     !PRINT *,'----------------------'                                                                                        !
     !IF ( STAT_NEWPOINT == 0 ) THEN
     !   PRINT *,'CONVERGED TO ', tuv, xyz
     !END IF
     !PRINT *,'STAT_NEWPOINT =',STAT_NEWPOINT                                                                                 !
     !                                                                                                                        !
     ! if we just found a degenerate point, return and report degeneracy                                                      !
     if ( stat_newpoint > 1 ) return                                                                                          !
     !                                                                                                                        !
     if ( stat_newpoint == 0 ) then ! <--------------------------------------------------------------------------------+      !
        ! check if the "new" point has not already been discovered                                                     !      !
        do ipt = 1,npts ! <-------------------------------------------------------------------+                        !      !
           if ( sum( ( xyz - coords(4:6,ipt) )**2 ) < EPSxyzsqr ) then ! <---------+          !                        !      !
              stat_newpoint = 1                                                    !          !                        !      !
              exit                                                                 !          !                        !      !
           end if ! <--------------------------------------------------------------+          !                        !      !
        end do ! <----------------------------------------------------------------------------+                        !      !
     end if ! <--------------------------------------------------------------------------------------------------------+      !
     !                                                                                                                        !
     if ( stat_newpoint == 0 ) then ! <-----------------------------------------------------------------+                     !
        ! if we just found a new intersection point, add it to the lists, and subdivide at that point   !                     !
        !PRINT *,'NEW POINT =',TUV, XYZ                                                                 !                     !
        !PRINT *,'IN   TBOX =',REGION_C%TBOX                                                            !                     !
        !PRINT *,'IN  UVBOX =',REGION_S%UVBOX                                                           !                     !
        !PRINT *,'';PRINT *,'';PRINT *,'';                                                              !                     !
        call append_vector( &                                                                           !                     !
             [tuv,xyz], &                                                                               !                     !
             6, &                                                                                       !                     !
             coords, &                                                                                  !                     !
             npts )                                                                                     !                     !
        !PRINT *,'COORDS ='                                                                             !                     !
        !CALL PRINT_MAT( TRANSPOSE(COORDS(:,1:NPTS)) )                                                  !                     !
        !                                                                                               !                     !
        ! add the new point to the current region and all their ascendants                              !                     !
        call add_point_bottom_up( region_c, npts )                                                      !                     !
        call add_point_bottom_up( region_s, npts )                                                      !                     !
        !                                                                                               !                     !
        tuv_subdiv = tuv                                                                                !                     !
     else ! --------------------------------------------------------------------------------------------+                     !
        ! else, subdivide at the parametric midpoint                                                    !                     !
        tuv_subdiv = 0.5_fp * [ &                                                                       !                     !
             region_c%uvbox(1) + region_c%uvbox(2), &                                                   !                     !
             region_s%uvbox([1,3]) + region_s%uvbox([2,4]) ]                                            !                     !
     end if ! <-----------------------------------------------------------------------------------------+                     !
  else ! ---------------------------------------------------------------------------------------------------------------------+
     IF ( N_SHAREDPTS < 1 ) STOP 'INTERIOR BUT N_SHAREDPTS < 1'                                                               !
     !tuv_subdiv = coords(1:3,sharedpts(1))                                                                                   !
     !PRINT *,'TUV_SUBDIV =',REAL( TUV_SUBDIV )                                                                               !
  end if ! <------------------------------------------------------------------------------------------------------------------+



  !PRINT *,'TUV_SUBDIV =',TUV_SUBDIV
  !! Subdivide the curve region
  call subdiv_region( &
       region_c, &
       tuv_subdiv(1), &
       stat_subdiv )

  !PRINT *,'STAT_SUBDIV C =',STAT_SUBDIV
  tuv_subdiv(1) = 0.5_fp * ( ab2n1p1( tuv_subdiv(1), region_c%uvbox(1), region_c%uvbox(2) ) + 1._fp )
  !PRINT *,'SIZE(CHILD) =',SIZE(REGION_C%CHILD)

  if ( stat_subdiv == 1 ) then ! <----------------------------------------------+
     ! the subdivision point is at one of the curve's endpoints                 !
     nchild(1) = 1                                                              !
  else ! -----------------------------------------------------------------------+
     ! the subdivision point is interior to the curve                           !
     nchild(1) = size(region_c%child)                                           !
     IF ( NCHILD(1) /= 2 ) THEN                                                 !
        PRINT *,'ERROR : NCHILD(1) =',NCHILD(1)                                 !
        STOP                                                                    !
     END IF                                                                     !
     !                                                                          !
     if ( stat_subdiv == 0 ) then ! <-----------------------------------+       !
        ! the curve region has no children yet                          !       !
        allocate( region_c%child(1)%poly(1), region_c%child(2)%poly(1) )      !       !
        allocate( region_c%child(1)%poly(1)%ptr, region_c%child(2)%poly(1)%ptr )
        call subdiv_bezier1( &                                          !       !
             region_c%poly(1)%ptr, &                                           !       !
             tuv_subdiv(1), &                                           !       !
             bl=region_c%child(1)%poly(1)%ptr, &                               !       !
             br=region_c%child(2)%poly(1)%ptr )                                !       !
     end if ! <---------------------------------------------------------+       !
     !                                                                          !
  end if ! <--------------------------------------------------------------------+


  !! Subdivide the surface region
  call subdiv_region( &
       region_s, &
       tuv_subdiv(2:3), &
       stat_subdiv )

  !PRINT *,'STAT_SUBDIV S =',STAT_SUBDIV
  tuv_subdiv(2) = 0.5_fp * ( ab2n1p1( tuv_subdiv(2), region_s%uvbox(1), region_s%uvbox(2) ) + 1._fp )
  tuv_subdiv(3) = 0.5_fp * ( ab2n1p1( tuv_subdiv(3), region_s%uvbox(3), region_s%uvbox(4) ) + 1._fp )

  if ( stat_subdiv == 2 ) then ! <----------------------------------------------+
     ! the subdivision point is at one of the surface's corners                 !
     nchild(2) = 1                                                              !
  else ! -----------------------------------------------------------------------+
     ! the subdivision point is interior to the surface                         !
     nchild(2) = size(region_s%child)                                           !
  end if ! <--------------------------------------------------------------------+

  if ( stat_subdiv == 0 ) then ! <------------------------------------------------------+
     ! 4 children                                                                       !
     allocate( &                                                                        !
          region_s%child(1)%poly(1), region_s%child(2)%poly(1), &                             !
          region_s%child(3)%poly(1), region_s%child(4)%poly(1) )                              !
     allocate( &                                                                        !
          region_s%child(1)%poly(1)%ptr, region_s%child(2)%poly(1)%ptr, &                             !
          region_s%child(3)%poly(1)%ptr, region_s%child(4)%poly(1)%ptr )                              !
     call subdiv_bezier2( &                                                             !
          region_s%poly(1)%ptr, &                                                              !
          tuv_subdiv(2:3), &                                                            !
          bsw=region_s%child(1)%poly(1)%ptr, &                                                 !
          bse=region_s%child(2)%poly(1)%ptr, &                                                 !
          bnw=region_s%child(3)%poly(1)%ptr, &                                                 !
          bne=region_s%child(4)%poly(1)%ptr )                                                  !
  elseif ( stat_subdiv == 1 ) then ! ---------------------------------------------------+
     ! 2 children                                                                       !
     allocate( region_s%child(1)%poly(1), region_s%child(2)%poly(1) )                         !
     allocate( region_s%child(1)%poly(1)%ptr, region_s%child(2)%poly(1)%ptr )                         !
     if ( region_s%child(2)%uvbox(1) <= region_s%uvbox(1) + EPSregion ) then ! <---+    !
        call subdiv_bezier2_only_v( &                                              !    !
             region_s%poly(1)%ptr, &                                                      !    !
             v=tuv_subdiv(3), &                                                    !    !
             bs=region_s%child(1)%poly(1)%ptr, &                                          !    !
             bn=region_s%child(2)%poly(1)%ptr )                                           !    !
     else ! -----------------------------------------------------------------------+    !
        call subdiv_bezier2_only_u( &                                              !    !
             region_s%poly(1)%ptr, &                                                      !    !
             u=tuv_subdiv(2), &                                                    !    !
             bw=region_s%child(1)%poly(1)%ptr, &                                          !    !
             be=region_s%child(2)%poly(1)%ptr )                                           !    !
     end if ! <--------------------------------------------------------------------+    !
     !                                                                                  !
  end if ! <----------------------------------------------------------------------------+



  !! Carry on the recursion with the children
  if ( all(nchild < 2) ) then  
     PRINT *,' TBOX =',REGION_C%UVBOX
     PRINT *,'UVBOX =',REGION_S%UVBOX
     PRINT *,N_SHAREDPTS,' SHARED POINTS'
     IF ( N_SHAREDPTS > 0) CALL PRINT_MAT( TRANSPOSE( COORDS(:,SHAREDPTS(1:N_SHAREDPTS) ) ) )
     PRINT *,'INTERIOR?',INTERIOR
     PRINT *,'separable?',separable
     PRINT *,'STAT_NEWPOINT =',STAT_NEWPOINT
     STAT_NEWPOINT = 99
     PRINT *,'intersect_curve_surface : NO MORE SUBDIVISION !!!!!'
     CALL WRITE_POLYNOMIAL( REGION_C%POLY(1)%ptr, 'dev_intersection_simple_surface/region_c_bezier.bern' )
     CALL WRITE_POLYNOMIAL( REGION_S%POLY(1)%ptr, 'dev_intersection_simple_surface/region_s_bezier.bern' )
     IF ( ASSOCIATED( REGION_C%XYZBOX ) ) THEN
        CALL WRITE_OBB( REGION_C%XYZBOX, 'dev_intersection_simple_surface/xyzbox_c.dat' )
     ELSE
        PRINT *,'OBB_C N/A'
     END IF
     IF ( ASSOCIATED( REGION_S%XYZBOX ) ) THEN
        CALL WRITE_OBB( REGION_S%XYZBOX, 'dev_intersection_simple_surface/xyzbox_s.dat' )
     ELSE
        PRINT *,'OBB_S N/A'
     END IF
     RETURN
  end if


  do jchild = 1,nchild(2) ! <--------------------------------+
     if ( nchild(2) == 1 ) then ! <-----------------+        !
        newregion_s => region_s                     !        !
     else ! ----------------------------------------+        !
        newregion_s => region_s%child(jchild)       !        !
     end if ! <-------------------------------------+        !
     !                                                       !
     do ichild = 1,nchild(1) ! <-------------------------+   !
        if ( nchild(1) == 1 ) then ! <--------------+    !   !
           newregion_c => region_c                  !    !   !
        else ! -------------------------------------+    !   ! 
           newregion_c => region_c%child(ichild)    !    !   !
        end if ! <----------------------------------+    !   !
        !                                                !   !
        call intersect_curve_surface( &                  !   !
             root_c, &                                   !   !
             root_s, &                                   !   !
             newregion_c, &                              !   !
             newregion_s, &                              !   !
             coords, &                                   !   !
             npts, &                                     !   !
             stat_degeneracy )                           !   !
        !                                                !   !
     end do ! <------------------------------------------+   !
     !                                                       !
  end do ! <-------------------------------------------------+

end subroutine intersect_curve_surface
