subroutine insert_polyline_point( &
     uv, &
     xyz, &
     pline, &
     i )
  use mod_math
  use mod_types_intersection
  ! Inserts a uv-xyz point in an intersection_polyline after the i-th point.
  ! If 'i' is not provided, the point is inserted after the last point.
  implicit none
  real(kind=MATHpr),                intent(in)           :: uv(2,2)
  real(kind=MATHpr),                intent(in)           :: xyz(3)
  type(type_intersection_polyline), intent(inout)        :: pline
  integer,                          intent(in), optional :: i
  integer                                                :: iprev, stat

  if ( present(i) ) then
     iprev = i
  else
     iprev = pline%np
  end if

  if ( iprev >= size(pline%xyz,2) ) then
     call reallocate_polyline( &
          pline, &
          pline%np + param_xtra_polyline_np, &
          stat )
     if ( stat > 0 ) STOP 'insert_polyline_point : could not reallocate polyline'
  end if

  if ( iprev < pline%np ) then
     pline%uv(:,:,iprev+2:pline%np+1) = pline%uv(:,:,iprev+1:pline%np)
     pline%xyz(:,iprev+2:pline%np+1) = pline%xyz(:,iprev+1:pline%np)
  end if
  pline%uv(:,:,iprev+1) = uv
  pline%xyz(:,iprev+1) = xyz
  pline%np = pline%np + 1

end subroutine insert_polyline_point
