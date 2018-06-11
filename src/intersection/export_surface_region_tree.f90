subroutine export_surface_region_tree( &
     root, &
     filename )
  use mod_util
  use mod_types_intersection
  ! Exports a surface_region subtree to a given file.
  ! (essentially for debugging purposes)
  ! See the 'write_surface_region' subroutine for details.
  implicit none
  type(type_surface_region), intent(in) :: root
  character(*),              intent(in) :: filename
  integer                               :: file_unit

  call get_free_unit( file_unit )

  open( &
       unit = file_unit, &
       file = filename, &
       action = 'write' )

  call write_surface_region( &
       root, &
       file_unit )

  close( file_unit )

  PRINT *,'surface_region tree written in ',filename

end subroutine export_surface_region_tree
