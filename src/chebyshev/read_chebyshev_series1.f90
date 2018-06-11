subroutine read_chebyshev_series1( &
     c, &
     filename )
  use mod_util
  implicit none
  character(*),                 intent(in)  :: filename
  type(type_chebyshev_series1), intent(out) :: c
  integer                                   :: degr, dim
  integer                                   :: fileunit, icoef, idim

  call get_free_unit( fileunit )
  if ( fileunit == 0 ) STOP "read_chebyshev_series1 : could not find free unit"

  open( unit = fileunit, file = filename, action = "read" )
  read (fileunit,*) degr, dim
  degr = degr - 1

  if ( degr < 0 ) STOP "read_chebyshev_series1 : degr < 0"

  call reset_chebyshev_series1( &
       c, &
       degr, &
       dim)

  do idim = 1,dim
     do icoef = 1,degr+1
        read (fileunit,*) c%coef(icoef,idim)
     end do
  end do
  close(fileunit)

end subroutine read_chebyshev_series1
