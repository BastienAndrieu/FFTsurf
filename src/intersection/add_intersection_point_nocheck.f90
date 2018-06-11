subroutine add_intersection_point_nocheck( &
    surf, &
    uv, &
    xyz, &
    interdat )
  use mod_math
  use mod_types_intersection
  ! Adds a new intersection_point without checking for duplicates
  implicit none
  type(ptr_parametric_surface),  intent(in)    :: surf(2)
  real(kind=MATHpr),             intent(in)    :: uv(2,2)
  real(kind=MATHpr),             intent(in)    :: xyz(3)
  type(type_intersection_data),  intent(inout) :: interdat
  type(type_point_on_surface), pointer         :: pos
  integer                                      :: isurf

  interdat%np = interdat%np + 1

  if ( interdat%np > size(interdat%points) ) then
     STOP 'add_intersection_point_nocheck : interdat%np > size(interdat%points)'
     ! reallocate ...
  end if

  interdat%points(interdat%np)%xyz = xyz
  allocate( interdat%points(interdat%np)%head )
  pos => interdat%points(interdat%np)%head
  pos%surface => surf(1)%ptr
  pos%uv = uv(:,1)
  pos%xyz = xyz
  
  allocate( pos%next )
  pos => pos%next
  pos%surface => surf(2)%ptr
  pos%uv = uv(:,2)
  pos%xyz = xyz
  nullify( pos%next )

end subroutine add_intersection_point_nocheck
