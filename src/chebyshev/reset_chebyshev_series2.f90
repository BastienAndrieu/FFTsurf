subroutine reset_chebyshev_series2( &
     c, &
     degr, &
     dim )
  implicit none
  integer,                      intent(in)  :: degr(2)
  integer,                      intent(in)  :: dim
  type(type_chebyshev_series2), intent(out) :: c

  if ( allocated(c%coef) ) then
     !if ( any(c%degr < degr) .or. c%dim < dim ) deallocate(c%coef)
     if ( size(c%coef,1) <= degr(1) .or. &
          size(c%coef,2) <= degr(2) .or. &
          size(c%coef,3) < dim ) deallocate(c%coef)
  end if
  if ( .not.allocated(c%coef) ) allocate(c%coef(1:degr(1)+1, 1:degr(2)+1, dim))

  c%degr = degr
  c%dim = dim
  c%coef(:,:,:) = 0._MATHpr

end subroutine reset_chebyshev_series2
