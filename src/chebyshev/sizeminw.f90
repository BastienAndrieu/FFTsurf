function sizeminw( n ) 
  ! Minimal size of work array (see dfftpack documentation, dcosti subroutine)
  implicit none
  integer, intent(in) :: n
  integer             :: sizeminw
  sizeminw = 3*n + 15
end function sizeminw
