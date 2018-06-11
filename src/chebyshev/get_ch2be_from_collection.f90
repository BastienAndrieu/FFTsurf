subroutine get_ch2be_from_collection( &
     collection, &
     N, &
     c2b )
  use mod_math
  implicit none
  integer,            intent(in)    :: N
  type(type_matrix),  intent(inout) :: collection(:)
  real(kind=MATHpr)                 :: c2b(N,N)

  if ( N > size(collection) ) then
     call ch2be_matrix( c2b, N-1 )
  else
     if ( .not.allocated(collection(N)%mat) ) then
        allocate( collection(N)%mat(N,N) )
        call ch2be_matrix( collection(N)%mat, N-1 )
     end if
     c2b = collection(N)%mat
  end if

end subroutine get_ch2be_from_collection
