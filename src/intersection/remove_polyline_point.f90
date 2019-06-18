subroutine remove_polyline_point( &
    polyline, &
    i, &
    stat )
  use mod_math
  use mod_types_intersection
  implicit none
  type(type_intersection_polyline), intent(inout) :: polyline
  integer, optional,                intent(in)    :: i
  integer,                          intent(out)   :: stat

  if ( i < 1 .or. i > polyline%np ) then
     stat = 1
     return
  end if

  if ( allocated(polyline%uv) ) then
     polyline%uv(1:2,1:2,i:polyline%np-1) = polyline%uv(1:2,1:2,i+1:polyline%np)
  end if

  if ( allocated(polyline%xyz) ) then
     polyline%xyz(1:3,i:polyline%np-1) = polyline%xyz(1:3,i+1:polyline%np)
  end if

  if ( allocated(polyline%s) ) then
     polyline%s(i:polyline%np-1) = polyline%s(i+1:polyline%np)
  end if

  polyline%np = polyline%np - 1
  stat = 0
  return
  
end subroutine remove_polyline_point
