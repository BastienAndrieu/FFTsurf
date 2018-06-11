subroutine init_surface_region( &
     region, &
     uvbox )
  use mod_math
  use mod_types_intersection
  ! Initializes a surface_region region with no child, given a uvbox
  implicit none
  real(kind=MATHpr),         intent(in)  :: uvbox(2,2)
  type(type_surface_region), intent(out) :: region

  region%uvbox = uvbox
  if ( associated(region%xyzbox) ) deallocate( region%xyzbox )
  if ( associated(region%pnbox) ) deallocate( region%pnbox )
  if ( associated(region%child) ) deallocate( region%parent )
  if ( associated(region%child) ) deallocate( region%child )
  nullify( region%xyzbox, region%pnbox, region%parent, region%child )

  if ( allocated(region%ipts_ss) ) deallocate(region%ipts_ss)
  if ( allocated(region%ipts_cs) ) deallocate(region%ipts_cs)

end subroutine init_surface_region
