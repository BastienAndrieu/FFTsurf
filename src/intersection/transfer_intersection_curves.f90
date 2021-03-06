subroutine transfer_intersection_curves( &
     from, &
     to )
  use mod_types_intersection
  implicit none
  type(type_intersection_data), intent(inout) :: from
  type(type_intersection_data), intent(inout) :: to
  integer                                     :: ic

  to%nc = from%nc

  do ic = 1,from%nc
     to%curves(ic)%dummy        =  from%curves(ic)%dummy
     to%curves(ic)%smooth       =  from%curves(ic)%smooth
     to%curves(ic)%surf(1)%ptr  => from%curves(ic)%surf(1)%ptr
     to%curves(ic)%surf(2)%ptr  => from%curves(ic)%surf(2)%ptr
     to%curves(ic)%uvbox        =  from%curves(ic)%uvbox
     to%curves(ic)%param_vector =  from%curves(ic)%param_vector
     to%curves(ic)%w0           =  from%curves(ic)%w0
     if ( allocated(from%curves(ic)%isplit) ) then
        call move_alloc(from=from%curves(ic)%isplit, to=to%curves(ic)%isplit)
     end if
     to%curves(ic)%nsplit       =  from%curves(ic)%nsplit
     if ( allocated(from%curves(ic)%iedge) ) then
        call move_alloc(from=from%curves(ic)%iedge, to=to%curves(ic)%iedge)
     end if
     to%curves(ic)%polyline     => from%curves(ic)%polyline
     nullify( &
          from%curves(ic)%surf(1)%ptr, &
          from%curves(ic)%surf(2)%ptr, &
          from%curves(ic)%polyline)
  end do

  deallocate(from%curves)
  from%nc = 0

end subroutine transfer_intersection_curves
