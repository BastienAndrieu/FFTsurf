program trace_border_adaptive

  use mod_math
  use mod_options
  use mod_util
  use mod_diffgeom
  use mod_eos
  
  implicit none

  type(type_options)         :: options
  type(type_surface)         :: surf
  integer                    :: iborder
  real(kind=fp), allocatable :: xyz(:,:), uv(:,:)
  integer                    :: np, stat
  integer                    :: fid
  character                  :: strnum
  integer                    :: i

  !! read options file
  call read_options( &
       'trace_border_adaptive/options.opt', &
       options )

  !! import surface
  call read_polynomial( &
       surf%x, &
       'trace_border_adaptive/surf_x.cheb', &
       nvar=2, &
       base=1 )
  call economize2( &
       surf%x, &
       EPSmath )
  call compute_deriv1(surf)
  call compute_deriv2(surf)


  !! Go!
  call get_free_unit(fid)

  do iborder = 1,4
     call trace_border_polyline_adaptive( &
          surf, &
          iborder, &
          options%chord_err, &
          options%hmin, &
          options%hmax, &
          stat, &
          uv, &
          xyz, &
          np )
     if ( stat > 0 ) then
        PRINT *,'!!! STAT =', STAT
        STOP
     end if

     write (strnum, '(i1)') iborder
     open( &
          unit=fid, &
          file='trace_border_adaptive/border_'//strnum//'.dat', &
          action='write' )
     do i = 1,np
        write (fid,*) xyz(1:3,i), uv(1:2,i)
     end do
     close(fid)
  end do

end program trace_border_adaptive
