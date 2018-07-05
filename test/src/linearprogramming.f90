program linearprogramming

  use mod_util
  use mod_math
  use mod_linearprogramming

  implicit none
  
  character(100)             :: nametest
  integer                    :: fid
  integer                    :: dim, nc, i
  real(kind=fp), allocatable :: c(:), A(:,:), x(:)
  logical                    :: singular
  integer                    :: stat

  if ( command_argument_count() < 1) then
     nametest = '1'
  else
     call get_command_argument(1, nametest)
  end if

  call get_free_unit( fid )
  open( &
       unit=fid, &
       file='linearprogramming/lpp' // trim(nametest) // '.dat', &
       action='read' )
  read (fid,*) !'dimension'
  read (fid,*) dim
  allocate( c(dim), x(dim) )
  read (fid,*) !'objectif'
  read (fid,*) c
  read (fid,*) !'contraintes'
  read (fid,*) nc
  allocate( A(nc,dim+1) )
  do i = 1,nc
     read (fid,*) A(i,:)
  end do
  close( fid )
  
  call solve_NxN( x, A(1:dim,1:dim), -A(1:dim,dim+1), singular )
  call lp_solve( &
       x, &
       stat, &
       A, &
       c, &
       dim, &
       nc )

  PRINT *,'STAT =',STAT
  IF ( STAT == 0 ) PRINT *,'X =',X

  open( &
       unit=fid, &
       file='linearprogramming/result.dat', &
       action='write' )
  write (fid,*) 'lpp' // trim(nametest) // '.dat'
  write (fid,*) stat
  if ( stat == 0 ) then
     write (fid,*) x
  end if
  close( fid )

end program linearprogramming
