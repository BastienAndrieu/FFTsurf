program decasteljau
  
  use mod_util
  use mod_math
  use mod_bernstein

  implicit none
  
  integer :: degr, dim

  real(kind=MATHpr), dimension(:,:), allocatable :: b, bl, br
  real(kind=MATHpr), dimension(:), allocatable   :: f
  real(kind=MATHpr)                              :: t
  integer                                        :: i, j, file_unit

  call get_free_unit( file_unit )
  open( unit=file_unit, file='decasteljau/b.dat', action='read' )
  read ( file_unit, * ) degr, dim
  degr = degr - 1
  allocate( b(degr+1,dim), bl(degr+1,dim), br(degr+1,dim), f(dim) )
  do j = 1,dim
     do i = 1,degr+1
        read ( file_unit, * ) b(i,j)
     end do
  end do
  close ( file_unit )

  open( unit=file_unit, file='decasteljau/t.dat', action='read' )
  read ( file_unit, * ) t
  close ( file_unit )


  call de_casteljau( &
       b, &
       t, &
       degr, &
       dim, &
       f, &
       bl, &
       br )
  !call de_casteljau( &
  !     b, &
  !     t, &
  !     degr, &
  !     dim, &
  !     bl=bl, &
  !     br=br )
  
  
  open( unit=file_unit, file='decasteljau/bl.dat', action='write' )
  do i = 1,degr+1
     write ( file_unit, * ) bl(i,:)
  end do
  close ( file_unit )

  open( unit=file_unit, file='decasteljau/br.dat', action='write' )
  do i = 1,degr+1
     write ( file_unit, * ) br(i,:)
  end do
  close ( file_unit )

  open( unit=file_unit, file='decasteljau/f.dat', action='write' )
  write ( file_unit, * ) f
  close ( file_unit )

end program decasteljau
