program debugobb

  use mod_obb
  
  implicit none
  
  type(type_obb) :: box1, box2
  logical        :: overlap

  call read_obb( box1, '/d/bandrieu/GitHub/FFTsurf/test/dev_intersection_simple_surface/xyzbox_c.dat' )
  call read_obb( box2, '/d/bandrieu/GitHub/FFTsurf/test/dev_intersection_simple_surface/xyzbox_s.dat' )
    
  call overlap_OBBs( &
       box1, &
       box2, &
       overlap )

contains

  subroutine read_obb( box, filename )
    use mod_util
    implicit none
    type(type_obb), intent(out) :: box
    character(*),   intent(in)  :: filename
    integer                     :: fid, i
    call get_free_unit( fid )
    open( unit=fid, file=filename, action='read' )
    read (fid,*) box%ctr
    read (fid,*) box%rng
    do i = 1,3
       read (fid,*) box%axe(:,i)
    end do
    close(fid)
  end subroutine read_obb

end program debugobb
