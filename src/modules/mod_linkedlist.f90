! Integer linked-lists
module mod_linkedlist

  implicit none

  type type_linkedlist
     integer                        :: i
     type(type_linkedlist), pointer :: next => null()
  end type type_linkedlist

contains

  subroutine init_linkedlist( &
       self, &
       i )
    implicit none
    type(type_linkedlist), pointer              :: self
    integer,               intent(in), optional :: i
    
    allocate( self )
    nullify( self%next )

    if ( present(i) ) self%i = i

  end subroutine init_linkedlist




  subroutine free_linkedlist( &
       self )
    implicit none
    type(type_linkedlist), pointer :: self
    type(type_linkedlist), pointer :: current, next
    
    current => self
    do while ( associated(current) )
       next => current%next
       deallocate( current )
       nullify( current )
       current => next
    end do

  end subroutine free_linkedlist
  


  
  subroutine insert_after( &
       self, &
       i )
    implicit none
    type(type_linkedlist), pointer              :: self
    integer,               intent(in), optional :: i
    type(type_linkedlist), pointer              :: next

    allocate( next )

    if ( present(i) ) next%i = i

    next%next => self%next
    self%next => next

  end subroutine insert_after



  
  recursive subroutine linkedlist_length( &
       self, &
       length )
    implicit none
     type(type_linkedlist), pointer       :: self
     integer,               intent(inout) :: length

     if ( .not.associated(self) ) return
     length = length + 1
     call linkedlist_length( &
          self%next, &
          length )
     
   end subroutine linkedlist_length



  subroutine linkedlist2array( &
       head, &
       array, &
       n_opt )
    implicit none
    type(type_linkedlist), pointer     :: head
    integer, allocatable,  intent(out) :: array(:)
    integer, optional,     intent(in)  :: n_opt
    integer                            :: n, i
    type(type_linkedlist), pointer     :: current

    n = 0
    if ( present(n_opt) ) then
       n = n_opt
    else
       call linkedlist_length( head, n )
    end if
    
    allocate( array(n) )
    current => head
    do i = 1,n
       array(i) = current%i
       current => current%next
    end do

  end subroutine linkedlist2array



  subroutine array2linkedlist( &
       array, &
       head )
    implicit none
    integer,               intent(in) :: array(:)
    type(type_linkedlist), pointer    :: head
    type(type_linkedlist), pointer    :: current
    integer                           :: i

    call init_linkedlist( &
         head, &
         array(1) )

    current => head
    do i = 2,size(array)
       call insert_after( &
            current, &
            array(i) )
       current => current%next
    end do

  end subroutine array2linkedlist


end module mod_linkedlist
