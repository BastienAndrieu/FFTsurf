recursive subroutine free_curve_region_tree( &
     region )
  use mod_types_intersection
  ! Recursively frees all the curve_regions in the subtree
  ! rooted at the node 'region'
  implicit none
  type(type_curve_region), intent(inout) :: region
  integer                                :: ichild

  if ( associated(region%child) ) then
     do ichild = 1,size(region%child)
        call free_curve_region_tree( region%child(ichild) )
     end do
     deallocate( region%child )
     nullify( region%child )
  end if

  if ( associated(region%xyzbox) ) deallocate(region%xyzbox)
  if ( allocated(region%ipts_cs) ) deallocate(region%ipts_cs)

end subroutine free_curve_region_tree
