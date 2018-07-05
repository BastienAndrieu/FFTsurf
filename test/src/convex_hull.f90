program convex_hull

  use mod_util
  use mod_math
  use mod_tolerances
  use mod_geometry

  implicit none
  
  integer, parameter   :: nmax = 10000
  integer              :: n
  real(kind=fp)        :: xy(2,nmax), rng(2)
  integer, allocatable :: hull(:)
  integer              :: nhull
  integer              :: fid, narg, i, j
  character(100)       :: nametest, arg
  logical              :: normalize
  character            :: strj
  integer*8            :: tic, toc, countrate

  narg = command_argument_count()
  if ( narg < 1) then
     nametest = ''
  else
     call get_command_argument(1, nametest)
     read ( nametest, * ) i
     if ( i < 1 ) nametest = ''
  end if

  if ( narg < 2) then
     normalize = .false.
  else
     call get_command_argument(2, arg)
     read ( arg, * ) normalize
  end if
  print *,'normalize?',normalize

  !print *,'EPS =',EPSmath
  
  call read_points( 'convex_hull/xy' // trim(nametest) //'.dat', xy, n, 2, nmax )
  PRINT *,'N =',N

  rng = maxval( xy(:,1:n), dim=2 ) - minval( xy(:,1:n), dim=2 )
  PRINT *,'RANGES =',RNG
  if ( normalize ) then
     do i = 1,2
        !if ( rng(i) < real( 2._fp*EPSmath, kind=fp) .and. &
        !     rng(i) > 0._fp ) xy(i,:) = xy(i,:) / rng(i)
        if ( rng(i) > EPSfp ) xy(i,:) = xy(i,:) / rng(i)
     end do
  end if

  allocate( hull(n) )
  do j = 1,2
     print *,''
     call system_clock( tic, countrate )
     if ( j == 1 ) then
        print *,'Gift-wrapping (Jarvis march)'
        call convex_hull_2d( xy, n, hull, nhull )
     else
        print *,'Quickhull'
        call quickhull( xy, n, hull, nhull )
     end if
     call system_clock( toc )
     print *,'elapsed :', real( toc - tic ) / real( countrate )
     !nhull = min( nhull, n )
     PRINT *,'NHULL =',NHULL
     !PRINT *,'HULL =',HULL(1:NHULL)

     write ( strj, '(I1)' ) j
     call get_free_unit( fid )
     open( &
          unit = fid, &
          file = 'convex_hull/hull' // strj // '.dat', &
          action = 'write' )
     do i = 1,nhull
        write (fid,*) hull(i)
     end do
     close( fid )
  end do


  call get_free_unit( fid )
  open( &
       unit = fid, &
       file = 'convex_hull/xy.dat', &
       action = 'write' )
  do i = 1,n
     write (fid,*) xy(:,i)
  end do
  close( fid )

  

contains

  subroutine read_points( filename, x, n, dim, nmax )
    implicit none
    character(*),  intent(in)    :: filename
    integer,       intent(in)    :: dim, nmax
    real(kind=fp), intent(inout) :: x(dim,nmax)
    integer,       intent(out)   :: n
    integer                      :: fid, io
    real(kind=fp)                :: p(dim)

    call get_free_unit( fid )
    open( &
         unit = fid, &
         file = filename, &
         action = 'read' )
    n = 0
    do
       read (fid,*,iostat=io) p
       if ( io /= 0 ) exit
       n = n + 1
       if ( n > nmax ) then
          print *,'read_points : /!\ n > nmax, stop reading'
          exit
       end if
       x(:,n) = p
    end do
    close( fid )
    
  end subroutine read_points

end program convex_hull
