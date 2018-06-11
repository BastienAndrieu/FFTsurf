recursive subroutine write_surface_region( &
     region, &
     file_unit )
  use mod_types_intersection
  ! Writes the data contained in a (leaf) surface_region.
  implicit none
  type(type_surface_region), intent(in) :: region
  integer,                   intent(in) :: file_unit
  integer                               :: ichild

  if ( associated(region%child) ) then
     ! the current region has child regions, carry on the recursion
     do ichild = 1,size(region%child)
        call write_surface_region( &
             region%child(ichild), &
             file_unit )
     end do
  else
     ! the current region is a leaf, write its data
     write (file_unit,*) region%uvbox
  end if

end subroutine write_surface_region
