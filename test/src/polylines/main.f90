  
  use mod_util
  use mod_math
  
  character                  :: strnum
  type(type_matrix)          :: xy(2)
  integer                    :: npts
  integer, allocatable       :: isegments(:,:)
  real(kind=fp), allocatable :: lambda(:,:)
  integer                    :: fid, i

  do i = 1,2
     write (strnum,'(i1)') i
     call read_matrix('polylines/xy' // strnum // '.dat', xy(i))
     call transpose_matrix(xy(i))
  end do

  stat = 1
  call intersect_2Dpolylines( &
       xy, &
       [1,1], &
       [size(xy(1)%mat,2), size(xy(2)%mat,2)], &
       npts, &
       isegments, &
       lambda )

  call get_free_unit(fid)
  open(unit=fid, file='polylines/result.dat', action='write')
  write(fid,*) npts
  do i = 1,npts
     write(fid,*) isegments(:,i)
  end do
  do i = 1,npts
     write(fid,*) lambda(:,i)
  end do

contains


  subroutine read_matrix( filename, a )
    implicit none
    character(*),      intent(in)    :: filename
    type(type_matrix), intent(inout) :: a
    integer                          :: m, n, i, fid

    fid = 1
    open(unit=fid, file=filename, action='read')
    read (fid,*) m, n
    if ( allocated(a%mat) ) deallocate(a%mat)
    allocate(a%mat(m,n))
    do i = 1,m
       read (fid,*) a%mat(i,1:n)
    end do
    close(fid)

  end subroutine read_matrix



  subroutine transpose_matrix( a )
    implicit none
    type(type_matrix), intent(inout) :: a
    real(kind=fp)                    :: tmp(size(a%mat,2),size(a%mat,1))

    tmp = transpose(a%mat)
    deallocate(a%mat)
    allocate(a%mat(size(tmp,1),size(tmp,2)))
    a%mat = tmp

  end subroutine transpose_matrix
