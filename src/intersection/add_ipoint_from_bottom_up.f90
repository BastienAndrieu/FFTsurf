recursive subroutine add_ipoint_from_bottom_up( &
     id, &
     region )
  use mod_types_intersection
  integer,                   intent(in)    :: id
  type(type_surface_region), intent(inout) :: region
  
  call append( region%ipts_ss, id )

  if ( associated( region%parent ) ) then
     call add_ipoint_from_bottom_up( &
          id, &
          region%parent )
  end if


end subroutine add_ipoint_from_bottom_up
