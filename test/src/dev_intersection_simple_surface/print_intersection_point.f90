subroutine print_intersection_point( point )
  implicit none
  type(type_intersection_point), intent(in) :: point
  type(type_point_on_surface), pointer      :: pos
  integer                                   :: ipos

  print *,'--------------------'
  print *,'xyz = ',point%xyz
  print *,'npos =',point%npos
  if ( point%npos > 0 ) then
     pos => point%pos
     do ipos = 1,point%npos
        print *,'uv =',pos%uv
        pos => pos%next
     end do
  end if
  print *,'--------------------'

end subroutine print_intersection_point
