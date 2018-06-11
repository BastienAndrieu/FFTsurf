subroutine compute_xyz_intersection_point( &
     point )
  use mod_math
  use mod_types_intersection
  ! Computes the mean xyz coordinates of several instances of
  ! point_on_surface representing the same intersection_point.
  implicit none
  type(type_intersection_point), intent(inout) :: point
  type(type_point_on_surface), pointer         :: pos
  integer                                      :: n

  point%xyz(:) = 0._MATHpr
  pos => point%head
  do while ( associated(pos) ) 
     n = n + 1
     point%xyz = point%xyz + pos%xyz
     pos => pos%next
  end do
  point%xyz = point%xyz / real( n, kind=MATHpr )

end subroutine compute_xyz_intersection_point
