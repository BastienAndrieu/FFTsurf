program separation

  use mod_math
  use mod_separation

  implicit none

  integer, parameter               :: nmax = 5000
  integer                          :: n1, n2
  real(kind=fp), dimension(nmax,3) :: sep1, sep2
  real(kind=fp)                    :: vec(3,2)
  logical                          :: separable(2)
  integer                          :: i


  call read_pts( &
       'separation/sep1.dat', &
       sep1, &
       n1, &
       nmax )

  call read_pts( &
       'separation/sep2.dat', &
       sep2, &
       n2, &
       nmax )



  
  call separating_plane( &
       sep1(1:n1,1:3), &
       sep2(1:n2,1:3), &
       n1, &
       n2, &
       vec(:,1), &
       separable(1) )

  call separating_plane( &
       sep1(1:n1,1:3), &
       -sep2(1:n2,1:3), &
       n1, &
       n2, &
       vec(:,2), &
       separable(2) )


  PRINT *,'SEPARABLE =',SEPARABLE
  DO I = 1,2
     IF ( SEPARABLE(I) ) PRINT *,'VEC',I,' =',VEC(:,I)
  END DO

contains

  subroutine read_pts( &
       filename, &
       pts, &
       n_pts, &
       n_pts_max )
  use mod_util
  implicit none
  integer,       intent(in)  :: n_pts_max
  character(*),  intent(in)  :: filename
  real(kind=fp), intent(out) :: pts(n_pts_max,3)
  integer,       intent(out) :: n_pts
  real(kind=fp)              :: point(3)
  integer                    :: fileunit, io

  call get_free_unit( fileunit )
  open( unit=fileunit, file=filename, action='read' )
  n_pts = 0
  do 
     read (fileunit, *, iostat=io) point
     if ( io /= 0 ) exit
     n_pts = n_pts + 1
     if (n_pts > n_pts_max) STOP 'read_pts : n_pts > max.'
     pts(n_pts,:) = point
  end do
  close(fileunit)

end subroutine read_pts


end program separation
