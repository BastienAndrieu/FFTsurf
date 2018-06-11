recursive subroutine intersect_simple_surfaces( &
     surfroot, &
     surfc, &
     region, &
     param_vector, &
     interdat )
  use mod_math
  use mod_chebyshev
  use mod_diffgeom
  use mod_types_intersection
  use mod_tolerances
  ! ...
  implicit none
  type(ptr_parametric_surface), intent(in)    :: surfroot(2)
  type(ptr_chebyshev_series2),  intent(in)    :: surfc(2)
  real(kind=MATHpr),            intent(in)    :: param_vector(3)
  type(ptr_surface_region),     intent(inout) :: region(2)
  type(type_intersection_data), intent(inout) :: interdat
  integer                                     :: icurv, ivar, ival
  
  return
  ! intersect the 4 borders of one surface with the other, and vice versa
  do icurv = 1,2

     do ivar = 1,2

        do ival = 1,2
           
           call intersect_border_surface( &
                surfroot, &
                surfc, &
                region, &
                icurv, &
                ivar, &
                ival, &
                interdat )

        end do

     end do

  end do
  
  
end subroutine intersect_simple_surfaces
         
