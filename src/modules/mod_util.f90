! =======================================================
! Module mod_util
! 29/05/18
! =======================================================

module mod_util

  implicit none

  interface append
     module procedure append_integer, append_double
  end interface append

  interface append_vec
     module procedure append_vec_i, append_vec_d
  end interface append_vec


  type type_list
     integer              :: n = 0
     integer, allocatable :: list(:)
  end type type_list

contains

  subroutine get_free_unit( iunit )
    ! Returns a free file unit for I/O
    implicit none
    integer, intent(out) :: iunit
    integer              :: i, stat
    logical              :: isopen
    iunit = 0
    do i = 1,99
       if ( i == 5 .or. i == 6 ) cycle

       inquire( unit=i, opened=isopen, iostat=stat )
       if ( stat /= 0 ) cycle
       if ( .not.isopen ) then
          iunit = i
          return
       end if
    end do
  end subroutine get_free_unit




  function int2logic(i)
    ! Converts an integer to a logical
    implicit none
    integer :: i
    logical :: int2logic

    if (i == 0) then
       int2logic = .false.
    else
       int2logic = .true.
    end if

  end function int2logic




  function logic2int(l)
    ! Converts a logical to an integer
    implicit none
    logical :: l
    integer :: logic2int

    if (l) then
       logic2int = 1
    else
       logic2int = 0
    end if

  end function logic2int




  subroutine is_in_list( &
       list, &
       val, &
       id )
    implicit none
    integer, intent(in)  :: list(:)
    integer, intent(in)  :: val
    integer, intent(out) :: id
    integer              :: i

    id = 0
    do i = 1,size(list)
       if ( list(i) == val ) then
          id = i
          return
       end if
    end do

  end subroutine is_in_list



  subroutine generate_name_from_date( &
       name )
    implicit none
    character(*), intent(out) :: name
    character(8)              :: date
    character(10)             :: time

    call date_and_time( date, time )
    name = date(7:8)//'_'//date(5:6)//'_'//date(1:4)//'_'//time(1:2)//'h'//time(3:4)

  end subroutine generate_name_from_date




  subroutine append_n( &
       list, &
       nlist, &
       elem, &
       nelem, &
       unique )
    implicit none
    integer, allocatable, intent(inout) :: list(:)
    integer,              intent(inout) :: nlist
    integer,              intent(in)    :: nelem
    integer,              intent(in)    :: elem(nelem)
    logical, optional,    intent(in)    :: unique
    logical                             :: mask(nelem)
    integer, allocatable                :: tmp(:)
    integer                             :: i, nnew


    if ( .not.allocated(list) ) then

       allocate( list(nelem) )
       list = elem
       nlist = nelem

    else

       mask(1:nelem) = .true.
       if ( present(unique) ) then
          if ( unique ) then
             do i = 1,nelem
                if ( any(list(1:nlist) == elem(i)) ) then
                   mask(i) = .false.
                   cycle
                end if
             end do
          end if
       end if

       nnew = count(mask)
       if ( nlist + nnew > size(list) ) then
          call move_alloc( from=list, to=tmp )
          allocate( list(nlist + nnew) )
          list(1:size(tmp)) = tmp
          deallocate( tmp )
       end if
       do i = 1,nelem
          if ( mask(i) ) then
             nlist = nlist + 1
             list(nlist) = elem(i)
          end if
       end do

    end if

    !PRINT *,'NLSIT=',NLIST
    !PRINT *,LIST
    !PRINT *,'------------------'
    !PRINT *,''

  end subroutine append_n





  subroutine append_list( &
       array, &
       length, &
       element, &
       unique )
    implicit none
    integer, allocatable, intent(inout) :: array(:)
    integer,              intent(inout) :: length
    integer,              intent(in)    :: element
    logical, optional,    intent(in)    :: unique
    integer, allocatable                :: tmp(:)

    if ( .not.allocated(array) ) then

       allocate( array(1) )
       array = element
       length = 1

    else

       if ( present(unique) ) then
          if ( unique ) then
             if ( any(array(1:length) == element) ) return
          end if
       end if

       length = length + 1
       if ( length > size(array) ) call move_alloc( from=array, to=tmp )
       allocate( array(length) )
       array(1:size(tmp)) = tmp
       array(length) = element
       deallocate(tmp)

    end if

  end subroutine append_list





  subroutine append_integer( &
       array, &
       element, &
       noduplicates, &
       newlength )
    ! Appends an integer element to a rank-1, dynamic array
    implicit none
    integer,              intent(in)    :: element
    integer, allocatable, intent(inout) :: array(:)
    logical, optional,    intent(in)    :: noduplicates
    integer, optional,    intent(out)   :: newlength
    integer, allocatable                :: tmp(:)
    integer                             :: length, i




    if ( allocated(array) ) then

       length = size(array)

       if ( present(noduplicates) ) then
          if ( noduplicates ) then
             do i = 1,length
                if ( array(i) == element ) then
                   if ( present(newlength) ) newlength = length
                   return
                end if
             end do
          end if
       end if

       allocate( tmp(length+1) )
       tmp(1:length) = array(1:length)
       tmp(length+1) = element

       deallocate( array )
       call move_alloc( FROM=tmp, TO=array )

    else

       allocate( array(1) )
       array(1) = element

    end if

    if ( present(newlength) ) newlength = size(array)

  end subroutine append_integer



  subroutine append_double( &
       array, &
       element )
    ! Appends a double precision real element to a rank-1, dynamic array
    implicit none
    real*8,              intent(in)    :: element
    real*8, allocatable, intent(inout) :: array(:)
    real*8, allocatable                :: tmp(:)
    integer                            :: length

    if ( allocated(array) ) then
       length = size(array)
       allocate( tmp(length+1) )
       tmp(1:length) = array(1:length)
       tmp(length+1) = element

       deallocate( array )
       call move_alloc( FROM=tmp, TO=array )
    else
       allocate( array(1) )
       array(1) = element
    end if

  end subroutine append_double



  subroutine union_arrays( &
       array1, &
       array2, &
       union )
    implicit none
    integer,              intent(in)  :: array1(:), array2(:)
    integer, allocatable, intent(out) :: union(:)
    integer, dimension(size(array1))  :: array1sorted
    integer, dimension(size(array2))  :: array2sorted
    integer                           :: i, j

    ! first, sort both arrays
    array1sorted = array1
    array2sorted = array2
    call bubblesort_integers( array1sorted )
    call bubblesort_integers( array2sorted )

    ! then merge them incrementally
    i = 1
    j = 1
    do while( i <= size(array1) .and. j <= size(array2) )

       if ( array1sorted(i) < array2sorted(j) ) then
          call append( union, array1sorted(i) )
          i = i + 1
       elseif ( array1sorted(i) > array2sorted(j) ) then
          call append( union, array2sorted(j) )
          j = j + 1
       else
          call append( union, array1sorted(i) )
          i = i + 1
          j = j + 1
       end if

    end do

    do while( i <= size(array1) )
       call append( union, array1sorted(i) )
       i = i + 1
    end do

    do while( j <= size(array2) )
       call append( union, array2sorted(j) )
       j = j + 1
    end do

  end subroutine union_arrays




  subroutine intersection_arrays( &
       array1, &
       array2, &
       intersection )
    implicit none
    integer,              intent(in)  :: array1(:), array2(:)
    integer, allocatable, intent(out) :: intersection(:)
    integer                           :: i

    call copy_array( &
         pack( array1, [ ( any( array1(i) == array2 ), i = 1, size(array1) ) ] ), &
         intersection )

  end subroutine intersection_arrays


  subroutine copy_array( &
       source_array, &
       target_array )
    implicit none
    integer,              intent(in)  :: source_array(:)
    integer, allocatable, intent(out) :: target_array(:)

    if ( size(source_array) > 0 ) then
       allocate( target_array(size(source_array)) )
       target_array = source_array
    end if

  end subroutine copy_array





  subroutine bubblesort_integers( a )
    implicit none
    integer, intent(inout) :: a(:)
    integer                :: tmp, i, j
    logical                :: swapped

    do j = size(a)-1,1,-1
       swapped = .false.
       do i = 1, j
          if ( a(i) > a(i+1) ) then
             tmp = a(i)
             a(i) = a(i+1)
             a(i+1) = tmp
             swapped = .true.
          end if
       end do
       if (.not. swapped) return
    end do

  end subroutine bubblesort_integers




  subroutine bubblesort_double( a )
    implicit none
    real*8, intent(inout) :: a(:)
    real*8                :: tmp
    integer               :: i, j
    logical               :: swapped

    do j = size(a)-1,1,-1
       swapped = .false.
       do i = 1, j
          if ( a(i) > a(i+1) ) then
             tmp = a(i)
             a(i) = a(i+1)
             a(i+1) = tmp
             swapped = .true.
          end if
       end do
       if (.not. swapped) return
    end do

  end subroutine bubblesort_double






  subroutine write_matrix( a, m, n, filename )
    implicit none
    integer,      intent(in) :: m, n
    real*8,       intent(in) :: a(m,n)
    character(*), intent(in) :: filename
    integer                  :: fileunit, i, j

    call get_free_unit( fileunit )
    open( unit=fileunit, file=filename, action='write' )
    do i = 1,m
       !write ( unit=fileunit, fmt=* ) a(i,:)
       do j = 1,n
          write ( unit=fileunit, fmt='(ES22.15,1x)', advance='no' ) a(i,j)
       end do
       write ( unit=fileunit, fmt=* )
    end do
    close( fileunit )
  end subroutine write_matrix





  subroutine randperm( perm, n )
    implicit none
    integer, intent(in)  :: n
    integer, intent(out) :: perm(n)
    integer              :: i, j, t

    perm = [ ( i, i=1,n ) ]

    do i = n,2,-1
       call random_integer( j, i-1 )
       j = j + 1
       t = perm(j)
       perm(j) = perm(i)
       perm(i) = t
    end do

  end subroutine randperm



  subroutine random_integer( i, imax )
    implicit none
    integer, intent(in)  :: imax
    integer, intent(out) :: i
    real                 :: r 

    call random_number( r )
    i = nint( real(imax) * r )

  end subroutine random_integer




  subroutine append_vec_i( &
       vec, &
       dim, &
       array, &
       n )
    use mod_math
    implicit none
    integer,              intent(in)    :: dim
    integer,              intent(in)    :: vec(dim)
    integer, allocatable, intent(inout) :: array(:,:)
    integer,              intent(inout) :: n
    integer, allocatable                :: tmp(:,:)

    if ( .not.allocated(array) ) allocate(array(dim,1))
    if ( dim > size(array,1) ) STOP 'append_vec_i : dim > size(array,1)'

    if ( n + 1 > size(array,2) ) then
       call move_alloc(from=array, to=tmp)
       allocate(array(dim,n+1))
       array(:,1:n) = tmp(:,1:n)
       deallocate(tmp)
    end if
    n = n + 1
    array(1:dim,n) = vec(1:dim)

  end subroutine append_vec_i




  subroutine append_vec_d( &
       vec, &
       dim, &
       array, &
       n )
    use mod_math
    implicit none
    integer,             intent(in)    :: dim
    real*8,              intent(in)    :: vec(dim)
    real*8, allocatable, intent(inout) :: array(:,:)
    integer,             intent(inout) :: n
    real*8, allocatable                :: tmp(:,:)

    if ( .not.allocated(array) ) allocate(array(dim,1))
    if ( dim > size(array,1) ) STOP 'append_vec_d : dim > size(array,1)'

    if ( n + 1 > size(array,2) ) then
       call move_alloc(from=array, to=tmp)
       allocate(array(dim,n+1))
       array(:,1:n) = tmp(:,1:n)
       deallocate(tmp)
    end if
    n = n + 1
    array(1:dim,n) = vec(1:dim)

  end subroutine append_vec_d





  subroutine insert_after( &
      list, &
      nlist, &
      elem, &
      iprev )
    implicit none
    integer, allocatable, intent(inout) :: list(:)
    integer,              intent(inout) :: nlist
    integer,              intent(in)    :: elem
    integer,              intent(inout) :: iprev
    integer, allocatable                :: tmp(:)

    if ( .not.allocated(list) ) then
       allocate(list(1))
       iprev = 0
       nlist = 1
       list(1) = elem
    else
       iprev = min(iprev, nlist)
       iprev = max(iprev, 0)
       if ( nlist + 1 > size(list) ) then
          call move_alloc( from=list, to=tmp )
          allocate(list(nlist+1))
          list(1:nlist) = tmp
          deallocate(tmp)
       end if
       list(iprev+2:nlist+1) = list(iprev+1:nlist)
       list(iprev+1) = elem
       nlist = nlist + 1
    end if

  end subroutine insert_after


  
  subroutine insert_n_after( &
       list, &
       nlist, &
       n, &
       elem, &
       iprev )
    implicit none
    integer, allocatable, intent(inout) :: list(:)
    integer,              intent(inout) :: nlist
    integer,              intent(in)    :: n
    integer,              intent(in)    :: elem(n)
    integer,              intent(inout) :: iprev
    integer, allocatable                :: tmp(:)

    if ( .not.allocated(list) ) then
       allocate(list(n))
       iprev = 0
       nlist = n
       list(1:n) = elem(1:n)
    else
       iprev = min(iprev, nlist)
       iprev = max(iprev, 0)
       if ( nlist + n > size(list) ) then
          call move_alloc( from=list, to=tmp )
          allocate(list(nlist+n))
          list(1:nlist) = tmp(1:nlist)
          deallocate(tmp)
       end if
       list(iprev+n+1:nlist+n) = list(iprev+1:nlist)
       list(iprev+1:iprev+n) = elem(1:n)
       nlist = nlist + n
    end if

  end subroutine insert_n_after





  subroutine insert_column_after( &
      array, &
      m, &
      n, &
      col, &
      iprev )
    implicit none
    integer, parameter                  :: n_xtra = 1
    integer, allocatable, intent(inout) :: array(:,:)
    integer,              intent(in)    :: m
    integer,              intent(inout) :: n
    integer,              intent(in)    :: col(m)
    integer,              intent(inout) :: iprev
    integer, allocatable                :: tmp(:,:)

    if ( .not.allocated(array) ) then
       allocate(array(m,n_xtra))
       iprev = 0
       n = 1
       array(1:m,1) = col(1:m)
    else
       iprev = min(iprev, n)
       iprev = max(iprev, 0)
       if ( n + 1 > size(array,2) ) then
          call move_alloc( from=array, to=tmp )
          allocate(array(max(m,size(tmp,1)),n + n_xtra))
          array(1:size(tmp,1),1:n) = tmp(1:size(tmp,1),1:n)
          deallocate(tmp)
       end if
       array(1:m,iprev+2:n+1) = array(1:m,iprev+1:n)
       array(1:m,iprev+1)     = col(1:m)
       n = n + 1
    end if

  end subroutine insert_column_after


  subroutine remove_from_list( &
       elem, &
       list, &
       nlist )
    implicit none
    integer, intent(in)    :: elem
    integer, intent(inout) :: list(:)
    integer, intent(inout) :: nlist
    integer                :: i

    i = 1
    do while ( i <= nlist )
       if ( list(i) == elem ) then
          !PRINT *,'REMOVING',ELEM
          list(i:nlist-1) = list(i+1:nlist)
          nlist = nlist - 1
       else
          i = i + 1
       end if
    end do
    
  end subroutine remove_from_list




  subroutine reallocate_list( &
       list, &
       n )
    implicit none
    integer,              intent(in)    :: n
    integer, allocatable, intent(inout) :: list(:)
    integer, allocatable                :: tmp(:)

    if ( .not.allocated(list) ) then
       allocate(list(n))
       return
    else
       call move_alloc(from=list, to=tmp)
       allocate(list(n))
       list(1:min(n,size(tmp))) = tmp(1:min(n,size(tmp)))
       deallocate(tmp)
    end if
    
  end subroutine reallocate_list
  
end module mod_util
