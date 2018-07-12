recursive subroutine intersect_surface_pair( &
     surfroot, &
     region, &
     interdata, &
     uvxyz, &
     nuvxyz, &
     stat_degeneracy )
  use mod_math
  use mod_bernstein2
  use mod_diffgeom2
  use mod_regiontree
  use mod_tolerances
  use mod_types_intersection
  implicit none
  type(ptr_surface),            intent(in)    :: surfroot(2)
  type(ptr_region),             intent(inout) :: region(2)
  type(type_intersection_data), intent(inout) :: interdata
  real(kind=fp), allocatable,   intent(inout) :: uvxyz(:,:)
  integer,                      intent(inout) :: nuvxyz
  integer,                      intent(inout) :: stat_degeneracy
  logical                                     :: overlap, separable
  integer                                     :: stat_collineal
  real(kind=fp)                               :: uv_collineal(2,2), xyz_collineal(3), n_collineal(3)
  real(kind=fp)                               :: uv_subdiv(2,2)
  type(ptr_polynomial)                        :: poly(2)
  real(kind=fp)                               :: param_vector(3)
  type(ptr_region)                            :: newregion(2)
  integer                                     :: isurf, ipt

  if ( stat_degeneracy > 1 ) return ! a degeneracy has been encountered

  ! inherit from parent regions all intersection points contained in current regions
  if ( nuvxyz > 0 ) then ! <-----------------+
     do isurf = 1,2 ! <----------------+     !
        call inherit_points( &         !     !
             region(isurf)%ptr, &      !     !
             uvxyz, &                  !     !
             nuvxyz )                  !     !
     end do ! <------------------------+     !
  end if ! <---------------------------------+
  
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
  n_collineal(:) = 0._fp
  call find_collineal_corners( &
       region, &
       stat_collineal, &
       uv_collineal, &
       xyz_collineal, &
       n_collineal )

  if ( stat_collineal <= 0 ) then ! <----------------------------+
     ! there is a pair of collineal corner points, we will       !
     ! subdivide both region at their parametric center point    !
     do isurf = 1,2 ! <----------------------------+             !
        uv_subdiv(:,isurf) = 0.5_fp * ( &          !             !
             region(isurf)%ptr%uvbox([1,3]) + &    !             !
             region(isurf)%ptr%uvbox([2,4]) )      !             !
     end do ! <------------------------------------+             !
     if ( stat_collineal < 0 ) then ! <-----------------------+  !
        ! there is a singular intersection point located at   !  !
        ! a corner in both surfaces, we add that intersection !  !
        ! point to the collection                             !  !
        call add_intersection_point( &                        !  ! *****
             uv, &                                            !  ! *****
             xyz, &                                           !  ! *****
             interdata, &                                     !  ! *****
             ipt )                                            !  ! *****
        do isurf = 1,2 ! <----------------+                   !  !
           call add_point_bottom_up( &    !                   !  !
                region(isurf)%ptr, &      !                   !  !
                ipt )                     !                   !  !
        end do ! <------------------------+                   !  !
        !                                                     !  !
        ! then we check if the pair of surfaces can intersect !  !
        ! at other points                                     !  !
        do isurf = 1,2 ! <----------------------------------+ !  !
           poly(isurf)%ptr => region(isurf)%ptr%poly(1)%ptr ! !  !
        end do ! <------------------------------------------+ !  !
        call intersect_elsewhere( &                           !  !
             poly, &                                          !  !
             xyz_collineal, &                                 !  !
             separable, &                                     !  !
             randomize=.true. )                               !  !
        if ( separable ) return                               !  !
        !                                                     !  !
     end if ! <-----------------------------------------------+  !
  end if ! <-----------------------------------------------------+



  ! check Hohmeyer's loop detection criterion
  do isurf = 1,2 ! <----------------------------------+
     poly(isurf)%ptr => region(isurf)%ptr%poly(1)%ptr !
  end do ! <------------------------------------------+
  call loop_detection_criterion( &
       poly, &
       param_vector, &
       stat_loopdetection, &
       ( stat_collineal <= 0 ), &
       n_collineal )
  
  if ( stat_loopdetection < 0 ) then
     ! a degeneracy has been encountered
     stat_degeneracy = stat_loopdetection
     return
  end if

  ! If the loop detection criterion IS satisfied, then the intersection of the 
  ! current surface regions is a (possibly empty) set of isolated points and/or 
  ! open curve branches. All the isolated points as well as the curves' endpoints
  ! are located on the boundary of a surface region.
  if ( stat_loopdetection == 0 ) then ! <---------------------------------------+
     ! first we start new temporary region trees, rooted at copies of current   !
     ! regions (all the data contained in the regions is copied, except the     !
     ! pointers to parents)                                                     !
     do isurf = 1,2 ! <-------------------+                                     !
        allocate( newregion(isurf)%ptr )  !                                     !
        call copy_region( &               !                                     !
             region(isurf)%ptr, &         !                                     !
             newregion(isurf)%ptr )       !                                     !
     end do ! <---------------------------+                                     !
     !                                                                          !
     call intersect_simple_surfaces( &                                          ! *****
          )                                                                     ! *****
     !                                                                          !
     ! copy the data in temporay regions back to the current regions            !
     ! (essentially the indices of newly discovered intersection points)...     !
     do isurf = 1,2 ! <------------------------------------------------------+  !
        if ( newregion(isurf)%ptr%npts > 0 ) then ! <---------------------+  !  !
           call append_n( &                                               !  !  !
                region(isurf)%ptr%ipts, &                                 !  !  !
                region(isurf)%ptr%npts, &                                 !  !  !
                newregion(isurf)%ptr%ipts(1:newregion(isurf)%ptr%npts), & !  !  !
                2, &                                                      !  !  !
                unique=.true. )                                           !  !  !
        end if ! <--------------------------------------------------------+  !  !
        ! ...then free the temporay region trees...                          !  !
        nullify( &                                                           !  !
             newregion(isurf)%ptr%xyzbox,      &                             !  !
             newregion(isurf)%ptr%poly(1)%ptr, &                             !  !
             newregion(isurf)%ptr%poly(2)%ptr )                              !  !
        deallocate( newregion(isurf)%ptr%poly )                              !  !
        call free_region_tree( newregion(isurf)%ptr )                        !  !
        deallocate( newregion(isurf)%ptr )                                   !  !
     end do ! <--------------------------------------------------------------+  !
     !                                                                          !
     ! ... and finally return (there cannot be other intersection points/curves ! 
     ! between the current surface regions)                                     !
     return                                                                     !
     !                                                                          !
  end if ! <--------------------------------------------------------------------+
  
end subroutine intersect_surface_pair
