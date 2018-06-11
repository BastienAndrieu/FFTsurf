subroutine intersect_border_surface( &
     surfroot, &
     surfc, &
     region, &
     icurv, &
     ivar, &
     ival, &
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
  type(ptr_surface_region),     intent(inout) :: region(2)
  integer,                      intent(in)    :: icurv, ivar, ival
  type(type_intersection_data), intent(inout) :: interdat
  type(type_parametric_curve)                 :: curv
  type(type_curve_region)                     :: regionc
  real(kind=MATHpr), allocatable              :: tuvxyz(:,:)
  integer, allocatable                        :: union_pts(:)
  real(kind=MATHpr)                           :: uv(2,2)
  integer                                     :: i, j, isurf

  ! convert surface border to parametric curve
  call cs2edge2cs1( &
       surfc(icurv)%ptr, &
       ivar, &
       ival, &
       curv%c )
  call compute_first_derivatives( curv )

  ! initialize curve region tree
  call init_curve_region( &
       regionc, &
       real( [-1,1], kind=MATHpr ) )

  ! intersect the curve and surface
  call intersect_curve_surface( &
       curv, &
       surfroot( 1 + mod(icurv,2) )%ptr, &
       regionc, &
       region( 1 + mod(icurv,2) )%ptr, &
       tuvxyz )

  call free_curve_region_tree( regionc )


  if ( allocated(tuvxyz) ) then

     ! check if the curve-surface intersection points have already been discovered
     call union_arrays( &
          region(1)%ptr%ipts_ss, &
          region(2)%ptr%ipts_ss, &
          union_pts )
     outer : do i = 1,size(tuvxyz)
        do j = 1,size(union_pts)
           if ( sum( (tuvxyz(4:6,i) - interdat%points(union_pts(j))%xyz)**2 ) &
                < EPSxyzsqr ) cycle outer 
        end do

        ! this is a new point, its uv coordinates in both root surfaces are :
        uv(:,1+mod(icurv,2)) = tuvxyz(2:3,i)

        uv(ivar,icurv) = region(icurv)%ptr%uvbox(ival,ivar)
        uv(1+mod(ivar,2),icurv) = n1p12ab( &
             tuvxyz(1,i), &
             region(icurv)%ptr%uvbox(1+mod(ivar,2),1), &
             region(icurv)%ptr%uvbox(1+mod(ivar,2),2) )
        !uv(1+mod(ivar,2),icurv) = 0.5_MATHpr * ( &
        !     region(icurv)%ptr%uvbox(1+mod(ivar,2),1) * ( 1._MATHpr - tuvxyz(1,i) ) + &
        !     region(icurv)%ptr%uvbox(1+mod(ivar,2),2) * ( 1._MATHpr + tuvxyz(1,i) ) )


        ! add this new intersection point to global list (interdat) ...
        call add_intersection_point_nocheck( &
             surfroot, &
             uv, &
             tuvxyz(4:6,i), &
             interdat )

        ! ... and to the local lists of each surface
        do isurf = 1,2
           call add_ipoint_from_bottom_up( &
                interdat%np, &
                region(isurf) )
           
        end do

     end do outer

  end if
  

end subroutine intersect_border_surface
