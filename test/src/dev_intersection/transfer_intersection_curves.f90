subroutine transfer_intersection_curves( &
     from, &
     to )
  implicit none
  type(type_intersection_data), intent(inout) :: from
  type(type_intersection_data), intent(inout) :: to
  integer                                     :: ic

  to%nc = from%nc

  do ic = 1,from%nc
     to%curves(ic)%smooth       =  from%curves(ic)%smooth
     to%curves(ic)%surf(1)%ptr  => from%curves(ic)%surf(1)%ptr
     to%curves(ic)%surf(2)%ptr  => from%curves(ic)%surf(2)%ptr
     to%curves(ic)%uvbox        =  from%curves(ic)%uvbox
     to%curves(ic)%param_vector =  from%curves(ic)%param_vector
     to%curves(ic)%root         => from%curves(ic)%root
     to%curves(ic)%polyline     => from%curves(ic)%polyline
     nullify(from%curves(ic)%surf(1)%ptr,from%curves(ic)%surf(2)%ptr)
     nullify(from%curves(ic)%root, from%curves(ic)%polyline)
  end do

  deallocate( from%curves )
  from%nc = 0

end subroutine transfer_intersection_curves
