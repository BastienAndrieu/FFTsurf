subroutine write_intersection_data( &
     interdata, &
     filepoints, &
     filecurves )
  use mod_util
  use mod_types_intersection
  implicit none
  type(type_intersection_data), intent(in) :: interdata
  character(*),                 intent(in) :: filepoints, filecurves
  integer                                  :: fileunit
  integer                                  :: ip, ic, is

  call get_free_unit( fileunit )

  ! intersection points
  open( &
       unit = fileunit, &
       file = filepoints, &
       action = 'write' )
  do ip = 1,interdata%np
     write ( fileunit, * ) interdata%points(ip)%xyz
  end do
  close( fileunit )

  ! intersection curves
  open( &
       unit = fileunit, &
       file = filecurves, &
       action = 'write' )
  write ( fileunit, * ) interdata%nc
  do ic = 1,interdata%nc
     write ( fileunit, * ) logic2int(interdata%curves(ic)%dummy)
     write ( fileunit, * ) logic2int(interdata%curves(ic)%smooth)
     write ( fileunit, * ) interdata%curves(ic)%uvbox(:,:,1)
     write ( fileunit, * ) interdata%curves(ic)%uvbox(:,:,2)
     write ( fileunit, * ) interdata%curves(ic)%nsplit
     do ip = 1,interdata%curves(ic)%nsplit
        write ( fileunit, * ) interdata%curves(ic)%isplit(:,ip)
     end do
     do is = 1,interdata%curves(ic)%nsplit-1
        write ( fileunit, * ) 1!class(is)
     end do
     if ( associated(interdata%curves(ic)%polyline) ) then
        write ( fileunit, * ) interdata%curves(ic)%polyline%np
        do ip = 1,interdata%curves(ic)%polyline%np
           write ( fileunit, * ) interdata%curves(ic)%polyline%uv(:,:,ip), interdata%curves(ic)%polyline%xyz(:,ip)
        end do
     else
        write ( fileunit, * ) 0
     end if
  end do
  close( fileunit )    

end subroutine write_intersection_data
