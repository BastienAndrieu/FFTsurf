subroutine subdiv_curve_region( &
     region, &
     ts, &
     stat )
  use mod_math
  use mod_types_intersection
  ! Subdivides a curve_region region at a given parametric point 'ts'.
  ! The region is subidivided into 2 child regions unless the subdivision point 
  ! is close to one of the region's endpoints.
  ! Only the 'tbox' attribute of the children are assigned 
  ! ( the remaining attributes require additional data that are outside the scope 
  ! of this subroutine ).
  implicit none
  real(kind=MATHpr),       intent(in)            :: ts
  type(type_curve_region), intent(inout), target :: region
  integer,                 intent(out)           :: stat
  real(kind=MATHpr)                              :: t(3)
  integer                                        :: ichild

  if ( associated( region%child ) ) return ! the region already has children

  t([1,3]) = region%tbox
  t(2) = ts

  ! check if the subivision point is close to one of the region's endpoints
  if ( abs(t(2) - t(1)) < EPSregion .or.  abs(t(3) - t(2)) < EPSregion ) then
     stat = 1
     return

  else
     ! subdivide the region into 2 child regions
     stat = 0
     allocate( region%child(2) )
     do ichild = 1,2
        call init_curve_region( &
             region, &
             t(ichild + [0,1]) )
     end do
     
  end if

end subroutine subdiv_curve_region
