program util

  use mod_util

  integer              :: n, i
  integer, allocatable :: list(:)

  print *,'n = ?'
  read (*,*) n
  allocate(list(n))

  list(1:n) = [(i, i=1,n)]
  print *,list

  print *,'ifirst = ?'
  read (*,*) i
  
  call circular_permutation( &
       list(1:n), &
       n, &
       i )
  print *,list
  
end program util
