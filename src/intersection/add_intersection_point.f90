subroutine add_intersection_point( &
     surf, &
     region, &
     uv, &
     interdat )
  use mod_math
  use mod_diffgeom
  use mod_types_intersection
  implicit none
  type(type_parametric_surface), intent(in), target :: surf(2)
  type(ptr_surface_region),      intent(in)         :: region(2)
  real(kind=MATHpr),             intent(in)         :: uv(2,2)
  type(type_intersection_data),  intent(inout)      :: interdat
























  !type(type_point_on_surface), target               :: pos(2)
  !type(type_point_on_surface), pointer              :: p
  !integer                                           :: stat(2)
  !integer                                           :: id
  !integer                                           :: isurf, i, j
  
  !id = interdat%np
  !do isurf = 1,2
  !   if ( .not.allocated( region(isurf)%ptr%ipoints ) ) cycle
  !   do j = 1,region(isurf)%ptr%npoints
  !      i = region(isurf)%ptr%ipoints(i)
  !      p => interdat%points(i)%head
  !   end do
  !end do





  !do isurf = 1,2
  !   pos(isurf)%surface => surf(isurf)
  !   pos(isurf)%uv = uv(:,isurf)
  !   call eval( pos(isurf)%xyz, surf(isurf), uv(:,isurf) )
  !end do
  !do i = 1,interdat%np
  !   p => interdat%points(i)%head
  !   do isurf = 1,2
  !      call compare_points_on_surface( &
  !           p, &
  !           pos(isurf), &
  !           stat(isurf) )
  !      if ( stat(isurf) == 0 ) continue ! distinct points
  !      if ( stat(isurf) == 1 ) then
  !         ! same point, same surface
  !      else
  !         ! same point, different surface
  !      end if
  !      id = i
  !      return
  !   end do
  !end do


  ! compare with existing points in interdat (only againts points within the intersected regions)
  ! add point if unique (and get id), else get id of same point
  
  ! add id of the point to list ipoints in region(1,2)%ptr
  ! increment counter npoints in region(1,2)%ptr
  ! if not allocated of size < id, (re)allocate ipoints

  


end subroutine add_intersection_point
