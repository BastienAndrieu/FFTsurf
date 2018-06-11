subroutine reallocate_polyline( &
     pline, &
     np, &
     stat )
  use mod_math
  use mod_types_intersection
  ! Reallocates the uv and xyz arrays of an intersection_polyline
  ! to size (2,2,np) and (3,np), respectively.
  ! The actual length of the polyline (pline%np) is kept unchanged.
  implicit none
  type(type_intersection_polyline), intent(inout) :: pline
  integer,                          intent(in)    :: np
  integer,                          intent(out)   :: stat
  real(kind=MATHpr)                               :: uv(2,2,pline%np)
  real(kind=MATHpr)                               :: xyz(3,pline%np)
  logical                                         :: nonempty

  stat = 0

  if ( np <= pline%np ) return

  nonempty = ( pline%np > 0) 

  if ( allocated(pline%uv) ) then
     if ( nonempty ) uv = pline%uv(:,:,1:pline%np)
     deallocate( pline%uv, STAT=stat )
     if ( stat > 0 ) STOP 'reallocate_polyline : could not deallocate uv'
  end if

  if ( allocated(pline%xyz) ) then
     if ( nonempty ) xyz = pline%xyz(:,1:pline%np)
     deallocate( pline%xyz, STAT=stat )
     if ( stat > 0 ) STOP 'reallocate_polyline : could not deallocate xyz'
  end if

  allocate( pline%uv(2,2,np), pline%xyz(3,np), STAT=stat )
  if ( stat > 0 ) return

  if ( nonempty ) then
     pline%uv(:,:,1:size(uv,3)) = uv
     pline%xyz(:,1:size(xyz,2)) = xyz
  end if

end subroutine reallocate_polyline
