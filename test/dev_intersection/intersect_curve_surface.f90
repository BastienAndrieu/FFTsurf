recursive subroutine intersect_curve_surface( &
     root_c, &
     root_s, &
     region_c, &
     region_s, &
     tuvxyz, &
     ntuvxyz )

  implicit none
  type(type_curve),           intent(in)            :: root_c
  type(type_surface),         intent(in)            :: root_s
  type(region),               intent(inout), target :: region_c
  type(region),               intent(inout), target :: region_s
  real(kind=fp), allocatable, intent(inout)         :: tuvxyz(:,:)
  integer,                    intent(inout)         :: ntuvxyz
  integer, allocatable                              :: ipts(:)
  integer                                           :: npts
  logical                                           :: overlap
  logical                                           :: interior(2)
  logical                                           :: separable
  real(kind=fp), dimension(3)                       :: tuv, xyz, tuv_subdiv
  integer                                           :: stat_newpoint
  integer                                           :: stat_subdiv
  integer                                           :: nchild(2)
  type(region), pointer                             :: newregion_c
  type(region), pointer                             :: newregion_s
  integer                                           :: i, j, k, ipt, jpt, ichild, jchild

  ! return if a degeneracy has been encountered previously
  ! (...)

  ! inherit from parents all the  points contained in the current regions
  if ( ntuvxyz > 0 ) then ! <-------------------------------------------------+
     call inherit_points( &                                                   !
          region_c, &                                                         !
          tuvxyz(1,1:ntuvxyz), &                                              !
          ntuvxyz )                                                           !
     call inherit_points( &                                                   !
          region_s, &                                                         !
          tuvxyz(2:3,1:ntuvxyz), &                                            !
          ntuvxyz )                                                           !
  end if ! <------------------------------------------------------------------+


  ! get list of already discovered points contained in both the current regions  
  npts = 0
  if ( region_c%npts > 0 .and. region_s%npts > 0 ) then ! <-------------------+
     call intersection_arrays( &                                              !
          region_c%ipts(1:region_c%npts), &                                   !
          region_s%ipts(1:region_s%npts), &                                   !
          ipts )                                                              !
     if ( allocated(ipts) ) npts = size(ipts)                                 !
  end if ! <------------------------------------------------------------------+


  ! if there are no already discovered points contained in both the current regions, check if 
  ! they intersect at endpoints/corners (2*4 tests)
  if ( npts < 1 ) then ! <--------------------------------------------------------------------+
     do k = 1,2 ! <----------------------------------------------------------------+          !
        do j = 1,2 ! <----------------------------------------------------------+  !          !
           xyz_s = region_s%poly(1)%ptr%coef( &                                 !  !          !
                1 + (j-1)*region_s%poly(1)%ptr%degr(1), &                       !  !          !
                1 + (k-1)*region_s%poly(1)%ptr%degr(2), &                       !  !          !
                1:3 )                                                           !  !          !
           do i = 1,2 ! <----------------------------------------------------+  !  !          !
              xyz_c = region_c%poly(1)%ptr%coef( &                           !  !  !          !
                   1 + (i-1)*region_c%poly(1)%ptr%degr(1), &                 !  !  !          !
                   1:3, &                                                    !  !  !          !
                   1 )                                                       !  !  !          !
              if ( sum( (xyz_s - xyz_c)**2 ) < EPSxyzsqr ) then ! <----+     !  !  !          !
                 xyz = 0.5_fp * ( xyz_c + xyz_s )                      !     !  !  !          !
                 tuv = real( [i,j,k]-1, kind=fp ) / real( &            !     !  !  !          !
                      [ region_c%poly(1)%ptr%degr(1), &                !     !  !  !          !
                      region_s%poly(1)%ptr%degr(1:2) ], &              !     !  !  !          !
                      kind=fp )                                        !     !  !  !          !
                 tuv = -1._fp + 2._fp * tuv                            !     !  !  !          !
                 !                                                     !     !  !  !          !
                 call append_vector( &                                 !     !  !  !          !
                      [tuv,xyz], &                                     !     !  !  !          !
                      6, &                                             !     !  !  !          !
                      tuvxyz, &                                        !     !  !  !          !
                      ntuvxyz )                                        !     !  !  !          !
                 !                                                     !     !  !  !          !
                 call append_n( &                                      !     !  !  !          !
                      ipts, &                                          !     !  !  !          !
                      npts, &                                          !     !  !  !          !
                      [ntuvxyz], &                                     !     !  !  !          !
                      1, &                                             !     !  !  !          !
                      unique=.true. )                                  !     !  !  !          !
                 !                                                     !     !  !  !          !
                 call add_point_bottom_up( region_c, ntuvxyz )         !     !  !  !          !
                 call add_point_bottom_up( region_s, ntuvxyz )         !     !  !  !          !
                 !                                                     !     !  !  !          !
              end if ! <-----------------------------------------------+     !  !  !          !
           end do ! <--------------------------------------------------------+  !  !          !
        end do ! <--------------------------------------------------------------+  !          !
     end do ! <--------------------------------------------------------------------+          !
  end if ! <----------------------------------------------------------------------------------+

  ! if there are still no intersection points, check if the bounding boxes overlap
  if ( npts < 1 ) then ! <--------------------------------------------------------------------+
     if ( .not.associated(region_c%xyzbox) ) then ! <--------+                                !
        allocate( region_c%xyzbox )                          !                                !
        call bernOBB1( &                                     !                                !
             region_c%poly(1)%ptr%coef( &                    !                                !
             1:region_c%poly(1)%ptr%degr(1)+1, &             !                                !
             1:3, &                                          !                                !
             1), &                                           !                                !
             region_c%poly(1)%ptr%degr(1), &                 !                                !
             region_c%xyzbox )                               !                                !
     end if ! <----------------------------------------------+                                !
     !                                                                                        !
     if ( .not.associated(region_s%xyzbox) ) then ! <--------+                                !
        allocate( region_s%xyzbox )                          !                                !
        call bernOBB2( &                                     !                                !
             region_s%poly(1)%ptr%coef( &                    !                                !
             1:region_s%poly(1)%ptr%degr(1)+1, &             !                                !
             1:region_s%poly(1)%ptr%degr(2)+1, &             !                                !
             1:3), &                                         !                                !
             region_s%poly(1)%ptr%degr, &                    !                                !
             region_s%xyzbox )                               !                                !
     end if ! <----------------------------------------------+                                !
     !                                                                                        !
     call overlap_OBBs( &                                                                     !
          region_c%xyzbox, &                                                                  !
          region_s%xyzbox, &                                                                  !
          overlap )                                                                           !
     !                                                                                        !
     if ( .not.overlap ) return ! disjoint bounding boxes => empty intersection               !
     !                                                                                        !
  else ! <------------------------------------------------------------------------------------+
     ! ( npts > 0 )                                                                           !
     ! are there any point interior the curve or surface region? ...                          !
     do jpt = 1,npts  ! <------------------------------------------+                          !
        ipt = ipts(jpt)                                            !                          !
        interior(1) = is_in_interval_strict( &                     !                          !
             tuvxyz(1,ipt), &                                      !                          !
             region_c%uvbox(1), region_c%uvbox(2) )                !                          !
        interior(2) = ( &                                          !                          !
             is_in_interval_strict( &                              !                          !
             tuvxyz(2,ipt), &                                      !                          !
             region_s%uvbox(1), region_s%uvbox(2) ) .and. &        !                          !
             is_in_interval_strict( &                              !                          !
             tuvxyz(3,ipt), &                                      !                          !
             region_s%uvbox(3), region_s%uvbox(4) ) )              !                          !
        !                                                          !                          !
        ! ... if so, subdivide both regions at that point          !                          !
        if ( any(interior) ) then ! <-------------------+          !                          !
           tuv_subdiv = tuvxyz(1:3,ipt)                 !          !                          !
           exit                                         !          !                          !
        end if ! <--------------------------------------+          !                          !
     end do ! <----------------------------------------------------+                          !
     !                                                                                        !
     if ( all(.not.interior) ) then  ! <--------------------------------------------------+   !
        ! if only one intersection point, check if the pair can intersect at other points !   !
        if ( npts == 1 ) then ! <-----------------------------------------------+         !   !
           call intersect_curve_surface_elsewhere( &                            !         !   !
                region_c%poly(1)%ptr, &                                         !         !   !
                region_s%poly(1)%ptr, &                                         !         !   !
                tuvxyz(4:6,ipts(1)), &                                          !         !   !
                separable, &                                                    !         !   !
                randomize=.true. )                                              !         !   !
           !                                                                    !         !   !
           if ( separable ) return ! the pair cannot intersect at other points  !         !   !
        end if ! <--------------------------------------------------------------+         !   !
     end if ! <---------------------------------------------------------------------------+   !
     !                                                                                        !
  end if ! <----------------------------------------------------------------------------------+


  ! search for a new intersection point
  if ( npts < 1 .or. all(.not.interior) ) then ! <--------------------------------------------+
     ! (...)                                                                                  !
     !                                                                                        !
     ! if a degeneracy has been encountered, report it                                        !
     ! (...)                                                                                  !
     !                                                                                        !
     ! if a point has been found, check whether it is a duplicate                             !
     if ( stat_newpoint == 0 ) then ! <---------------------------------------------------+   !
        do ipt = 1,ntuvxyz ! <----------------------------------------------+             !   !
           if ( sum( (xyz - tuvxyz(4:6,ipt))**2 ) < EPSxyzsqr ) then ! <--+ !             !   !
              stat_newpoint = 1                                           ! !             !   !
              call add_point_bottom_up( region_c, ipt )                   ! !             !   !
              call add_point_bottom_up( region_s, ipt )                   ! !             !   !
              exit                                                        ! !             !   !
           end if ! <-----------------------------------------------------+ !             !   !
        end do ! <----------------------------------------------------------+             !   !
        !                                                                                 !   !
        ! if this is actually a new intersection point, append it to the tuvxyz list      !   !
        call append_vector( &                                                             !   !
             [tuv,xyz], &                                                                 !   !
             6, &                                                                         !   !
             tuvxyz, &                                                                    !   !
             ntuvxyz )                                                                    !   !
        !                                                                                 !   !
        ! append that point to the current regions' lists and to that of their ascendants !   !
        call add_point_bottom_up( region_c, ntuvxyz )                                     !   !
        call add_point_bottom_up( region_s, ntuvxyz )                                     !   !
        !                                                                                 !   !
        ! subdivide both regions at that point                                            !   !
        tuv_subdiv = tuv                                                                  !   !
     end if  ! <--------------------------------------------------------------------------+   !
     !                                                                                        !
     ! if no new point has been discovered, subdivide both regions at their centerpoint       !
     if ( stat_newpoint /= 0 ) then ! <-----------------------------------------------+       !
        tuv_subdiv(1) = 0.5_fp * ( region_c%uvbox(1) + region_c%uvbox(2) )            !       !
        tuv_subdiv(2:3) = 0.5_fp * ( region_s%uvbox([1,3]) + region_s%uvbox([2,4]) )  !       !
     end if ! <-----------------------------------------------------------------------+       !
     !                                                                                        !
  end if ! <----------------------------------------------------------------------------------+


  !! subdivide the curve region
  call subdiv_region( &
       region_c, &
       tuv_subdiv(1), &
       stat_subdiv )
  ! linear change of variable -> local frame of Bézier subcurve
  tuv_subdiv(1) = 0.5_fp * ( ab2n1p1( tuv_subdiv(1), region_c%uvbox(1), region_c%uvbox(2) ) + 1._fp )

  if ( stat_subdiv == 1 ) then ! <--------------------------------------------------+
     ! the subdivision point is at one of the curve's endpoints                     !
     nchild(1) = 1                                                                  !
  else ! ---------------------------------------------------------------------------+
     ! the subdivision point is interior to the curve                               !
     nchild(1) = size(region_c%child)                                               !
     if ( stat_subdiv == 0 ) then ! <-------------------------+                     !
        ! the curve region has no children yet                !                     !
        do ichild = 1,size(region_c%child) ! <-------------+  !                     !
           allocate( region_c%child(ichild)%poly(1) )      !  !                     !
           allocate( region_c%child(ichild)%poly(1)%ptr )  !  !                     !
        end do ! <-----------------------------------------+  !                     !
        call subdiv_bezier1( &                                !                     !
             region_c%poly(1)%ptr, &                          !                     !
             tuv_subdiv(1), &                                 !                     !
             bl=region_c%child(1)%poly(1)%ptr, &              !                     !
             br=region_c%child(2)%poly(1)%ptr )               !                     !
     end if ! <-----------------------------------------------+                     !
  end if ! <------------------------------------------------------------------------+



  !! subdivide the surface region
  call subdiv_region( &
       region_s, &
       tuv_subdiv(2:3), &
       stat_subdiv )
  ! linear change of variables -> local frame of Bézier subsurface
  tuv_subdiv(2) = 0.5_fp * ( ab2n1p1( tuv_subdiv(2), region_s%uvbox(1), region_s%uvbox(2) ) + 1._fp )
  tuv_subdiv(3) = 0.5_fp * ( ab2n1p1( tuv_subdiv(3), region_s%uvbox(3), region_s%uvbox(4) ) + 1._fp )

  if ( stat_subdiv == 2 ) then ! <--------------------------------------------------+
     ! the subdivision point is at one of the surface's corners                     !
     nchild(2) = 1                                                                  !
  else ! ---------------------------------------------------------------------------+
     ! the subdivision point is interior to the surface                             !
     nchild(2) = size(region_s%child)                                               !
     if ( stat_subdiv >= 0 ) then ! <-------------------------+                     !
        ! the surface region has no children yet              !                     !
        do ichild = 1,size(region_s%child) ! <-------------+  !                     !
           allocate( region_s%child(ichild)%poly(1) )      !  !                     !
           allocate( region_s%child(ichild)%poly(1)%ptr )  !  !                     !
        end do ! <-----------------------------------------+  !                     !
     end if ! <-----------------------------------------------+                     !
  end if ! <------------------------------------------------------------------------+

  if ( stat_subdiv == 0 ) then ! <--------------------------------------------------+
     ! subdivide into 4 children                                                    !
     call subdiv_bezier2( &                                                         !
          region_s%poly(1)%ptr, &                                                   !
          tuv_subdiv(2:3), &                                                        !
          bsw=region_s%child(1)%poly(1)%ptr, &                                      !
          bse=region_s%child(2)%poly(1)%ptr, &                                      !
          bnw=region_s%child(3)%poly(1)%ptr, &                                      !
          bne=region_s%child(4)%poly(1)%ptr )                                       !
  elseif ( stat_subdiv == 1 ) then ! -----------------------------------------------+
     ! subdivide into 2 children                                                    !
     if ( region_s%child(2)%uvbox(1) <= region_s%uvbox(1) + EPSregion ) then ! <--+ !
        call subdiv_bezier2_only_v( &                                             ! !
             region_s%poly(1)%ptr, &                                              ! !
             v=tuv_subdiv(3), &                                                   ! !
             bs=region_s%child(1)%poly(1)%ptr, &                                  ! !
             bn=region_s%child(2)%poly(1)%ptr )                                   ! !
     else ! ----------------------------------------------------------------------+ !
        call subdiv_bezier2_only_u( &                                             ! !
             region_s%poly(1)%ptr, &                                              ! !
             u=tuv_subdiv(2), &                                                   ! !
             bw=region_s%child(1)%poly(1)%ptr, &                                  ! !
             be=region_s%child(2)%poly(1)%ptr )                                   ! !
     end if ! <-------------------------------------------------------------------+ !
  end if ! <------------------------------------------------------------------------+

  ! carry on the recursion with new pairs of regions
  if ( all(nchild < 2) ) then
     ! erreur
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
        ! (...)
        !                                                !   !
     end do ! <------------------------------------------+   !
     !                                                       !
  end do ! <-------------------------------------------------+

end subroutine intersect_curve_surface
