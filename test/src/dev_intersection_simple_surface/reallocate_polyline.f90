subroutine reallocate_polyline( &
     polyline, &
     np, &
     stat_alloc )
  use mod_math
  ! Reallocates the uv and xyz arrays of an intersection_polyline
  ! to size (2,2,np) and (3,np), respectively.
  ! The actual length of the polyline (pline%np) is kept unchanged.
  implicit none
  type(type_intersection_polyline), intent(inout) :: polyline
  integer,                          intent(in)    :: np
  integer,                          intent(out)   :: stat_alloc
  real(kind=fp), allocatable                      :: uvtmp(:,:,:), xyztmp(:,:)

  if ( allocated(polyline%uv) ) then ! <--------------------------------------------+
     if ( size(polyline%uv,3) < np ) then ! <-----------------------------------+   !
        call move_alloc( polyline%uv, uvtmp )                                   !   !
        allocate( polyline%uv(2,2,np), stat=stat_alloc )                        !   !
        if ( stat_alloc == 0 ) polyline%uv(1:2,1:2,1:size(uvtmp,3)) = uvtmp     !   !
     else ! --------------------------------------------------------------------+   !
        stat_alloc = 0                                                          !   !
     end if ! <-----------------------------------------------------------------+   !
  else ! ---------------------------------------------------------------------------+
     allocate( polyline%uv(2,2,np), stat=stat_alloc )                               !
  end if ! <------------------------------------------------------------------------+

  if ( stat_alloc /= 0 ) then ! <-----------+
     stat_alloc = 1                         !
     return                                 !
  end if ! <--------------------------------+

  
  if ( allocated(polyline%xyz) ) then ! <-------------------------------------------+
     if ( size(polyline%xyz,2) < np ) then ! <----------------------------------+   !
        call move_alloc( polyline%xyz, xyztmp )                                 !   !
        allocate( polyline%xyz(3,np), stat=stat_alloc )                         !   !
        if ( stat_alloc == 0 ) polyline%xyz(1:3,1:size(xyztmp,2)) = xyztmp      !   !
     else ! --------------------------------------------------------------------+   !
        stat_alloc = 0                                                          !   !
     end if ! <-----------------------------------------------------------------+   !
  else ! ---------------------------------------------------------------------------+
     allocate( polyline%xyz(3,np), stat=stat_alloc )                                !
  end if ! <------------------------------------------------------------------------+

  if ( stat_alloc /= 0 ) then ! <-----------+
     stat_alloc = 2                         !
     return                                 !
  end if ! <--------------------------------+

end subroutine reallocate_polyline
