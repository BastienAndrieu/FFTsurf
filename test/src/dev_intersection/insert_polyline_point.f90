subroutine insert_polyline_point( &
    uv, &
    xyz, &
    stat, &
    polyline, &
    i )
  use mod_math
  ! Inserts a uv-xyz point in an intersection_polyline after the i-th point.
  ! If 'i' is not provided, the point is inserted after the current last point.
  implicit none
  integer, parameter                              :: PARAM_xtra_np = 10
  real(kind=fp),                    intent(in)    :: uv(2,2)
  real(kind=fp),                    intent(in)    :: xyz(3)
  integer,                          intent(out)   :: stat
  type(type_intersection_polyline), intent(inout) :: polyline
  integer, optional,                intent(in)    :: i
  integer                                         :: iprev

  if ( present(i) ) then
     iprev = i
     iprev = min(iprev,polyline%np)
     iprev = max(iprev,0)
  else
     iprev = polyline%np
  end if

  if ( .not.allocated(polyline%uv)  ) allocate(polyline%uv(2,2,PARAM_xtra_np))
  if ( .not.allocated(polyline%xyz) ) allocate(polyline%xyz(3,PARAM_xtra_np) )

  if ( polyline%np + 1 > size(polyline%xyz,2) ) then
     call reallocate_polyline( &
          polyline, &
          polyline%np + PARAM_xtra_np, &
          stat )
     if ( stat > 0 ) return
  end if

  stat = 0
  if ( iprev < polyline%np ) then
     polyline%uv(:,:,iprev+2:polyline%np+1) = polyline%uv(:,:,iprev+1:polyline%np)
     polyline%xyz(:,iprev+2:polyline%np+1)  = polyline%xyz(:,iprev+1:polyline%np)
  end if
  polyline%uv(:,:,iprev+1) = uv
  polyline%xyz(:,iprev+1)  = xyz
  polyline%np = polyline%np + 1

end subroutine insert_polyline_point
