subroutine write_bernstein_series2( &
     b, &
     filename )
  use mod_util
  implicit none
  character(*),                 intent(in) :: filename
  type(type_bernstein_series2), intent(in) :: b
  integer                                  :: fileunit, icoef, jcoef, idim

  call get_free_unit( fileunit )
  if ( fileunit == 0 ) STOP "write_bernstein_series2 : could not find free unit"

  open( unit = fileunit, file = filename, action = "write" )
  write (fileunit,*) b%degr+1, b%dim

  do idim = 1,b%dim
     do jcoef = 1,b%degr(2)+1
        do icoef = 1,b%degr(1)+1
           write (fileunit,"(ES22.15)") b%coef(icoef,jcoef,idim)
        end do
     end do
  end do
  close(fileunit)
end subroutine write_bernstein_series2
