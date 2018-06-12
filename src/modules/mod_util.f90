! =======================================================
! Module mod_util
! 29/05/18
! =======================================================

module mod_util

  implicit none

  interface append
     module procedure append_integer, append_double
  end interface append

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




  !subroutine intersection_lists( &
  !     list1, &
  !     list2, &
  !     inter, &
  !     ninter )
  !  implicit none
  !  integer, intent(in)  :: list1(:), list2(:)
  !  integer, intent(out) :: inter(:), ninter
  !  integer              :: i, j
  !  ninter = 0
  !  do i = 1,size(list1)
  !     do j = 1,size(list2)
  !       if ( list1(i) == list2(j) ) then
  !          ninter = ninter + 1
  !          if ( ninter > size(inter) ) then
  !             ninter = -1
  !             return
  !          end if
  !          inter(ninter) = list1(i)
  !       end if
  !     end do
  !  end do
  !end subroutine intersection_lists


  
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



end module mod_util
