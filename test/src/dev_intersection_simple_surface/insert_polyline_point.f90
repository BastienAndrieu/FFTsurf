subroutine insert_polyline_point( &
    uv, &
    xyz, &
    polyline, &
    i )
  use mod_math
  ! Inserts a uv-xyz point in an intersection_polyline after the i-th point.
  ! If 'i' is not provided, the point is inserted after the current last point.
  implicit none
  integer, parameter                              :: PARAM_xtra_np = 10
  real(kind=fp),                    intent(in)    :: uv(2,2)
  real(kind=fp),                    intent(in)    :: xyz(3)
  type(type_intersection_polyline), intent(inout) :: polyline
  integer, optional,                intent(in)    :: i
  integer                                         :: iprev, stat_alloc
  
  if ( present(i) ) then
     iprev = i
  else
     iprev = polyline%np
  end if

  if ( iprev >= size(polyline%xyz,2) ) then
     call reallocate_polyline( &
          polyline, &
          polyline%np + PARAM_xtra_np, &
          stat_alloc )
     if ( stat_alloc == 1 ) then
        STOP 'insert_polyline_point : could not reallocate polyline%uv'
     elseif ( stat_alloc == 2 ) then
        STOP 'insert_polyline_point : could not reallocate polyline%xyz'
     end if
  end if

  if ( iprev < polyline%np ) then
     polyline%uv(:,:,iprev+2:polyline%np+1) = polyline%uv(:,:,iprev+1:polyline%np)
     polyline%xyz(:,iprev+2:polyline%np+1) = polyline%xyz(:,iprev+1:polyline%np)
  end if
  polyline%uv(:,:,iprev+1) = uv
  polyline%xyz(:,iprev+1) = xyz
  polyline%np = polyline%np + 1

end subroutine insert_polyline_point
