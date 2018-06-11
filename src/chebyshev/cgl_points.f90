subroutine cgl_points( &
     x, &
     a, &
     b, &
     n )
    ! Returns the n Chebyshev-Gauss-Lobatto points
    ! mapped to the interval [a,b]
    implicit none
    integer,           intent(in)  :: n
    real(kind=MATHpr), intent(in)  :: a, b
    real(kind=MATHpr), intent(out) :: x(n)
    integer                        :: i = 0

    if (n == 1) then
       x = 0.5_MATHpr*(a+b)
    else
       x = a + (b-a) * 0.5_MATHpr * ( 1._MATHpr + &
            cos(MATHpi* [( real(i,kind=MATHpr)/real(n-1,kind=MATHpr) , &
            i = n-1,0,-1 )]) )
    end if

  end subroutine cgl_points
