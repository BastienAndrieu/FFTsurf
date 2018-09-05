program graph

  use mod_util
  use mod_math
  use mod_graph
  
  implicit none

  character(99)              :: arg, fname
  character(2)               :: strnum
  integer                    :: num, fid
  integer                    :: narc, nnod, m, nloops
  integer, allocatable       :: arc2nod(:,:), looparc(:,:), loopnod(:,:), lenloop(:)
  real(kind=fp), allocatable :: arc_angles(:,:)
  integer                    :: iarc, i

  if ( command_argument_count() < 1 ) then
     fname = '/d/bandrieu/stck/Bureau/FFTsurf2/test/graph/data_test_loop.dat'
  else
     call get_command_argument(1, arg)
     read (arg,*) num
     write (strnum, '(I2.2)') num
     fname = '/d/bandrieu/GitHub/FFTsurf/test/graph/test'//strnum//'.dat'
  end if
  
  ! read embedded graph
  call get_free_unit(fid)
  open(unit=fid, file=trim(fname), action='read')
  
  read (fid,*) narc
  allocate(arc2nod(2,narc), arc_angles(2,narc))
  do iarc = 1,narc
     read (fid,*) arc2nod(:,iarc)
  end do
  do iarc = 1,narc
     read (fid,*) arc_angles(:,iarc)
  end do

  read (fid,*) nnod
  
  close(fid)

  
  
  ! make loops
  m = min(narc, nnod)
  allocate(looparc(m,m), loopnod(m,m), lenloop(m))
  call make_loops( &
       arc2nod, &
       arc_angles, &
       narc, &
       nnod, &
       nloops, &
       looparc, &
       !loopnod, &
       lenloop )

  do i = 1,nloops
     do iarc = 1,lenloop(i)
        loopnod(iarc,i) = arc2nod(1,looparc(iarc,i))
     end do
  end do


  
  open(unit=fid, file='graph/result.dat', action='write')
  write (fid,'(A)') trim(fname)
  write (fid,*) nloops
  do i = 1,nloops
     write (fid,*) looparc(1:lenloop(i),i)
     write (fid,*) loopnod(1:lenloop(i),i)
  end do
  close(fid)
  
end program graph
