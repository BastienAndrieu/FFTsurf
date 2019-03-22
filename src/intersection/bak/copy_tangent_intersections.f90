subroutine copy_tangent_intersections( &
     from, &
     to )
  use mod_types_intersection
  implicit none
  type(type_intersection_data), intent(in), target :: from
  type(type_intersection_data), intent(inout)      :: to
  type(type_point_on_surface), pointer             :: pos => null()
  real(kind=fp), allocatable                       :: uv(:,:)
  type(ptr_surface), allocatable                   :: surf(:)
  integer                                          :: icurv, j, ipoint, iendpoints(2), isurf

  do icurv = 1,from%nc
     if ( .not.from%curves(icurv)%smooth ) cycle
     
     do j = 1,2
        if ( j == 1 ) then
           ipoint = from%curves(icurv)%isplit(1,1)
        else
           ipoint = from%curves(icurv)%isplit(1,from%curves(icurv)%nsplit)
        end if

        allocate(uv(2,from%points(ipoint)%npos), surf(from%points(ipoint)%npos))
        pos => from%points(ipoint)%pos
        do isurf = 1,from%points(ipoint)%npos
           uv(:,isurf) = pos%uv 
           surf(isurf)%ptr => pos%surf
           pos => pos%next
        end do        
        
        call add_intersection_point( &
             uv(1:2,1:from%points(ipoint)%npos), &
             from%points(ipoint)%xyz, &
             surf(1:from%points(ipoint)%npos), &
             from%points(ipoint)%npos, &
             to, &
             iendpoints(j) )
        deallocate(uv, surf)
        
     end do

     call add_intersection_curve( &
          to, &
          from%curves(icurv)%param_vector, &
          iendpoints, &
          from%curves(icurv)%uvbox )

     to%curves(to%nc)%dummy = from%curves(icurv)%dummy
     to%curves(to%nc)%smooth = from%curves(icurv)%smooth
     do isurf = 1,2
        to%curves(to%nc)%surf(isurf)%ptr => from%curves(icurv)%surf(isurf)%ptr
     end do
     to%curves(to%nc)%polyline => from%curves(icurv)%polyline
     if ( allocated(to%curves(to%nc)%iedge) ) deallocate(to%curves(to%nc)%iedge)
     to%curves(to%nc)%isplit(2,1:2) = from%curves(icurv)%isplit(2,1:2)
     
  end do
  
end subroutine copy_tangent_intersections
