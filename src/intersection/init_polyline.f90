subroutine init_polyline( &
     pline, &
     stat )
  use mod_types_intersection
  ! Allocates the uv and xyz arrays of an intersection_polyline
  ! to size (2,2,param_init_polyline_np) and (3,param_init_polyline_np), 
  ! respectively. The actual length of the polyline (pline%np) is set to zero.
  implicit none
  type(type_intersection_polyline), intent(inout) :: pline
  integer,                          intent(out)   :: stat

  stat = 0

  if ( allocated(pline%uv) ) then
     deallocate( pline%uv, STAT=stat )
     if ( stat > 0 ) STOP 'init_polyline : could not deallocate uv'
  end if

  if ( allocated(pline%xyz) ) then
     deallocate( pline%xyz, STAT=stat )
     if ( stat > 0 ) STOP 'init_polyline : could not deallocate xyz'
  end if

  allocate( &
       pline%uv(2,2,param_init_polyline_np), &
       pline%xyz(3,param_init_polyline_np), &
       STAT=stat )
  if ( stat > 0 ) return

  pline%np = 0

end subroutine init_polyline
