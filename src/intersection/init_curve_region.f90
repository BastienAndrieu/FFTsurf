subroutine init_curve_region( &
     region, &
     tbox )
  use mod_math
  use mod_types_intersection
  ! Initializes a curve_region region with no child, given an interval tbox
  implicit none
  real(kind=MATHpr),       intent(in)  :: tbox(2)
  type(type_curve_region), intent(out) :: region

  region%tbox = tbox
  if ( associated(region%xyzbox) ) deallocate( region%xyzbox )
  if ( associated(region%child) ) deallocate( region%child )
  nullify( region%xyzbox, region%child )
  
  if ( allocated(region%ipts_cs) ) deallocate(region%ipts_cs)

end subroutine init_curve_region
