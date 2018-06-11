subroutine read_chebyshev_series2( &
     c, &
     filename )
  use mod_util
  implicit none
  character(*),                 intent(in)  :: filename
  type(type_chebyshev_series2), intent(out) :: c
  integer                                   :: degr(2), dim
  integer                                   :: fileunit, icoef, jcoef, idim

  call get_free_unit( fileunit )
  if ( fileunit == 0 ) STOP "read_chebyshev_series2 : could not find free unit"

  open( unit = fileunit, file = filename, action = "read" )
  read (fileunit,*) degr, dim
  degr = degr - 1

  if ( any(degr < 0) ) STOP "read_chebyshev_series2 : degr < 0"

  call reset_chebyshev_series2( &
       c, &
       degr, &
       dim)

  do idim = 1,dim
     do jcoef = 1,degr(2)+1
        do icoef = 1,degr(1)+1
           read (fileunit,*) c%coef(icoef,jcoef,idim)
        end do
     end do
  end do
  close(fileunit)

end subroutine read_chebyshev_series2
