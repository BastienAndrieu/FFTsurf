subroutine compare_points_on_surface( &
     p1, &
     p2, &
     stat )
  use mod_types_intersection
  use mod_tolerances
  ! Compares two instances of point_on_surface. 
  ! Returns stat=0 if the points are not coincident (within the "same-point tolerance").
  ! Else, stat=1 if the points lie ont the same surface, 
  ! and stat=2 if they do not.
  implicit none
  type(type_point_on_surface), intent(in)  :: p1, p2
  integer,                     intent(out) :: stat

  stat = 0
  if ( associated(p1%surface, p2%surface) ) then
     if ( maxval( abs(p1%uv - p2%uv) ) < p1%uvtol + p2%uvtol ) then
        stat = 1
        return
     end if
  end if

  if ( sum( (p1%xyz - p2%xyz)**2 ) < EPSxyzsqr ) then
     stat = 2
  end if

end subroutine compare_points_on_surface
