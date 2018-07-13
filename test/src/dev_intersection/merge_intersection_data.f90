subroutine merge_intersection_data( &
     surf, &
     uvxyz, &
     nuvxyz, &
     interdata_local, &
     interdata_global )
  use mod_math
  use mod_types_intersection
  ! Trace all intersection curves, intersect them and subidivide them accordingly and 
  implicit none
  type(ptr_surface),            intent(in)    :: surf(2)
  integer,                      intent(in)    :: nuvxyz
  real(kind=fp),                intent(in)    :: uvxyz(7,nuvxyz)
  type(type_intersection_data), intent(in)    :: interdata_local
  type(type_intersection_data), intent(inout) :: interdata_global
  integer                                     :: id_global(nuvxyz)
  integer                                     :: ip, ic

  ! add new intersection points
  do ip = 1,nuvxyz
     call add_intersection_point( &
          reshape(uvxyz(1:4,ip), [2,2]), &
          uvxyz(5:7,ip), &
          surf, &
          interdata_global, &
          id_global(ip) ) 
  end do

  do ic = 1,interdata_local%nc
     ! add curve
     call add_intersection_curve( &
          interdata_global, &
          interdata_local%curves(ic)%param_vector, &
          id_global(interdata_local%curves(ic)%root%endpoints), &
          interdata_local%curves(ic)%uvbox )

     ! trace polyline
     
  end do





end subroutine merge_intersection_data
