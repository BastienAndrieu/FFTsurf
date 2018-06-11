program test_union_intersection_arrays

  use mod_util

  integer                             :: n1, n2
  integer, dimension(:),  allocatable :: a1, a2, u, i


  print *,'n1 ?'
  read (*,*) n1
  allocate( a1(n1) )
  if ( n1 > 0 ) then
     print *,'a1 ?'
     read (*,*) a1
  end if

  print *,'n2 ?'
  read (*,*) n2
  allocate( a2(n2) )
  if ( n2 > 0 ) then
     print *,'a2 ?'
     read (*,*) a2
  end if

  call union_arrays( a1, a2, u )
  print *,'       union =', u
  
  call intersection_arrays( a1, a2, i )
  print *,'intersection =', i

end program test_union_intersection_arrays
