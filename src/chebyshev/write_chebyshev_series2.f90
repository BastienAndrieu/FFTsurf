subroutine write_chebyshev_series2( &
     c, &
     filename )
  use mod_util
  implicit none
  character(*),                 intent(in) :: filename
  type(type_chebyshev_series2), intent(in) :: c
  integer                                  :: fileunit, icoef, jcoef, idim

  call get_free_unit( fileunit )
  if ( fileunit == 0 ) STOP "write_chebyshev_series2 : could not find free unit"

  open( unit = fileunit, file = filename, action = "write" )
  write (fileunit,*) c%degr+1, c%dim

  do idim = 1,c%dim
     do jcoef = 1,c%degr(2)+1
        do icoef = 1,c%degr(1)+1
           write (fileunit,"(ES22.15)") c%coef(icoef,jcoef,idim)
        end do
     end do
  end do
  close(fileunit)
end subroutine write_chebyshev_series2
