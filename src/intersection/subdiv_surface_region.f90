subroutine subdiv_surface_region( &
     region, &
     uvs, &
     stat )
  use mod_math
  use mod_types_intersection
  ! Subdivides a surface_region region at a given parametric point 'uvs'.
  ! If the subdivision point is close to one of the region's boundaries, 
  ! the region is subdivided into 2 child regions, else, it is subdivided 
  ! into 4 child regions.
  ! Only the 'uvbox' and 'parent' attributes of the children are assigned 
  ! ( the remaining attributes require additional data that are outside the scope 
  ! of this subroutine ).
  implicit none
  real(kind=MATHpr),         intent(in)            :: uvs(2)
  type(type_surface_region), intent(inout), target :: region
  integer,                   intent(out)           :: stat
  real(kind=MATHpr)                                :: uv(3,2)
  logical                                          :: degenerate(2)
  integer                                          :: idim, ichild, jchild

  if ( associated( region%child ) ) return ! the region already has children

  uv([1,3],:) = region%uvbox
  uv(2,:) = uvs

  ! check if the subivision point is close to the region's boundary
  do idim = 1,2
     degenerate(idim) = ( abs(uv(2,idim) - uv(1,idim)) < EPSregion ) .or. &
          ( abs(uv(3,idim) - uv(2,idim)) < EPSregion )
  end do

  if ( all(degenerate) ) then
     ! degenerate subdivision (corner)
     stat = 3
     return
     !STOP 'subdiv_surface_region : degenerate subdivision (corner)'

  elseif ( degenerate(1) ) then
     ! subdivision into 2 child regions along horizontal line v = vs
     stat = 1
     allocate( region%child(2) )
     do jchild = 1,2
        call init_surface_region( &
             region%child(jchild), &
             reshape( [uv([1,3],1), uv(jchild+[0,1],2)], [2,2] ) )
     end do

  elseif ( degenerate(2) ) then
     ! subdivision into 2 child regions along vertical line u = us
     stat = 2
     allocate( region%child(2) )
     do ichild = 1,2
        call init_surface_region( &
             region%child(ichild), &
             reshape( [uv(ichild+[0,1],1), uv([1,3],2)], [2,2] ) )
     end do

  else
     ! subdivision into 4 child regions at point (u,v) = (us,vs)
     stat = 0
     allocate( region%child(4) )
     do jchild = 1,2
        do ichild = 1,2
           call init_surface_region( &
                region%child( 2*(jchild-1) + ichild ), &
                reshape( [uv(ichild+[0,1],1), uv(jchild+[0,1],2)], [2,2] ) )
        end do
     end do

  end if

  ! make all children point to their parent, i.e. the current region
  do ichild = 1,size(region%child)
     region%child(ichild)%parent => region
  end do
  

end subroutine subdiv_surface_region
