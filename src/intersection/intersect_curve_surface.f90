! if the curve and surface do not intersect except at already discovered points, return

! else, if there are no already discovered points of intersection interior to the curve or surface :
!    check bounding boxes and perform numerical iteration to find an intersection point :
!        if a degenerate point of intersection is found :
!            ....
!        else, if a new intersection point is found :
!            subdivide the curve and surface at that point and carry on the recursion
!        else, id no point of intersection is found :
!            subdivide the curve and surface at their parametric midpoint and carry on the recursion
! else, if a 

recursive subroutine intersect_curve_surface( &
    rootc, &
    roots, &
    regionc, &
    regions, &
    tuvxyz )
  use mod_math
  !use mod_list
  use mod_util
  use mod_diffgeom
  use mod_types_intersection
  implicit none
  type(type_parametric_curve),    intent(in)    :: rootc
  type(type_parametric_surface),  intent(in)    :: roots
  type(type_curve_region),        intent(inout) :: regionc
  type(type_surface_region),      intent(inout) :: regions
  real(kind=MATHpr), allocatable, intent(inout) :: tuvxyz(:)
  integer, allocatable                          :: shared_pts(:), nshared_pts
  
  ! do the curve and surface intersect at already discovered points?
  nshared_pts = 0
  !if ( npoints > 0 ) then
  if ( allocated(tuvxyz) ) then
     if ( size(tuvxyz) > 0 ) then
        if ( allocated( regionc%ipts_cs ) .and. allocated( regions%ipts_cs ) ) then
           call intersection_arrays( regionc%ipts_cs, regions%ipts_cs, shared_pts )
           nshared_pts = size( shared_pts )
        end if
     end if
  end if
  
  if ( nshared_pts == 1 ) then
     
     
  elseif ( nshared_pts > 1 ) then
     
     
  else
     
     
  end if

end subroutine intersect_curve_surface
