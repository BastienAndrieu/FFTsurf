subroutine reset_bernstein_series1( &
     b, &
     degr, &
     dim )
  implicit none
  integer,                      intent(in)  :: degr
  integer,                      intent(in)  :: dim
  type(type_bernstein_series1), intent(out) :: b

  if ( allocated(b%coef) ) then
     if ( size(b%coef,1) <= degr .or. &
          size(b%coef,2) < dim ) deallocate(b%coef)
  end if
  if ( .not.allocated(b%coef) ) allocate(b%coef(1:degr+1, dim))

  b%degr = degr
  b%dim = dim
  b%coef(:,:) = 0._MATHpr

end subroutine reset_bernstein_series1
