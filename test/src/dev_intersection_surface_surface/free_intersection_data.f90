subroutine free_intersection_data( &
     interdat )
  implicit none
  type(type_intersection_data), intent(inout) :: interdat
  type(type_point_on_surface), pointer        :: current => null(), next => null()
  integer                                     :: ip, ic
  
  if ( allocated(interdat%points) ) then

     do ip = 1,interdat%np
        interdat%points(ip)%npos = 0
        current => interdat%points(ip)%pos
        do while ( associated(current) )
           nullify( current%surf )
           next => current%next
           deallocate( current )
           nullify( current )
           current => next
        end do
     end do

     deallocate( interdat%points )
     interdat%np = 0

  end if





  if ( allocated(interdat%curves) ) then

     do ic = 1,interdat%nc
        if ( allocated(interdat%curves(ic)%polyline%s) ) deallocate( interdat%curves(ic)%polyline%s )
        if ( allocated(interdat%curves(ic)%polyline%uv) ) deallocate( interdat%curves(ic)%polyline%uv )
        if ( allocated(interdat%curves(ic)%polyline%xyz) ) deallocate( interdat%curves(ic)%polyline%xyz )
        interdat%curves(ic)%polyline%np = 0
     end do

     deallocate( interdat%curves )
     interdat%nc = 0

  end if

end subroutine free_intersection_data
