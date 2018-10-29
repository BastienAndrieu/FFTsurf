subroutine read_intersection_curves( &
     filename, &
     surf, &
     interdata )
  use mod_util
  use mod_diffgeom
  use mod_types_intersection
  implicit none
  character(*),                 intent(in)    :: filename
  type(type_surface), target,   intent(in)    :: surf(:)
  type(type_intersection_data), intent(inout) :: interdata
  integer                                     :: fid, ncurves, np
  real(kind=fp), allocatable                  :: xyz(:,:), uv(:,:,:)
  type(ptr_surface)                           :: surfpair(2)
  integer                                     :: ic, ip, pair(2), endpt(2)

  call get_free_unit(fid)
  open( &
       unit=fid, &
       file=filename, &
       action='read' )
  read (fid,*) ncurves
  do ic = 1,ncurves
     read (fid,*) pair
     surfpair(1)%ptr => surf(pair(1))
     surfpair(2)%ptr => surf(pair(2))
     read (fid,*) np

     allocate(xyz(3,np), uv(2,2,np))
     do ip = 1,np
        read (fid,*) xyz(:,ip)
     end do
     do ip = 1,np
        read (fid,*) uv(:,1,ip), uv(:,2,ip)
     end do

     call add_intersection_point( &
          uv(:,:,1), &
          xyz(:,1), &
          surfpair, &
          2, &
          interdata, &
          endpt(1) ) 
     call add_intersection_point( &
          uv(:,:,np), &
          xyz(:,np), &
          surfpair, &
          2, &
          interdata, &
          endpt(2) )

     call add_intersection_curve( &
          interdata, &
          [0._fp, 0._fp, 0._fp], &
          endpt, &
          spread(spread([-1._fp, 1._fp], 2, 2), 3, 2) )
     interdata%curves(interdata%nc)%surf(1)%ptr => surf(pair(1))
     interdata%curves(interdata%nc)%surf(2)%ptr => surf(pair(2))
     interdata%curves(interdata%nc)%isplit(2,:) = [1,np]

     allocate(interdata%curves(interdata%nc)%polyline)
     interdata%curves(interdata%nc)%polyline%np = np
     call move_alloc(from=xyz, to=interdata%curves(interdata%nc)%polyline%xyz)
     call move_alloc(from=uv , to=interdata%curves(interdata%nc)%polyline%uv )

  end do
  close(fid)

end subroutine read_intersection_curves
