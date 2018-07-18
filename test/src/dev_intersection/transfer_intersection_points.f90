subroutine transfer_intersection_points( &
     from, &
     to )
  use mod_types_intersection
  implicit none
  type(type_intersection_data), intent(inout) :: from
  type(type_intersection_data), intent(inout) :: to
  integer                                     :: ip
  
  to%np = from%np

  do ip = 1,from%np
     to%points(ip)%xyz  =  from%points(ip)%xyz
     to%points(ip)%pos  => from%points(ip)%pos
     to%points(ip)%npos =  from%points(ip)%npos
     nullify(from%points(ip)%pos)
  end do
  deallocate(from%points)
  from%np = 0

end subroutine transfer_intersection_points
