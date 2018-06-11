subroutine reset_bernstein_series2( &
     b, &
     degr, &
     dim )
  implicit none
  integer,                      intent(in)  :: degr(2)
  integer,                      intent(in)  :: dim
  type(type_bernstein_series2), intent(out) :: b

  if ( allocated(b%coef) ) then
     if ( size(b%coef,1) <= degr(1) .or. &
          size(b%coef,2) <= degr(2) .or. &
          size(b%coef,3) < dim ) deallocate(b%coef)
  end if
  if ( .not.allocated(b%coef) ) allocate( b%coef(1:degr(1)+1, 1:degr(2)+1, dim) )

  b%degr = degr
  b%dim = dim
  b%coef(:,:,:) = 0._MATHpr

end subroutine reset_bernstein_series2
