!program dev_intersection

  use mod_util
  use mod_math
  use mod_polynomial
  use mod_diffgeom2  
  use mod_types_intersection

  implicit none

  LOGICAL, PARAMETER :: GLOBALDEBUG = .TRUE.
  
  integer                    :: narg
  character(100)             :: arg
  integer                    :: numtest
  character                  :: strnum
  character(2)               :: strnum2
  
  logical, parameter                      :: ECONOMIZE = .true.
  integer                                 :: nsurf
  type(type_surface), target, allocatable :: surf(:)
  logical, allocatable                    :: mask(:)
  type(type_intersection_data)            :: interdata_global
  integer                                 :: isurf

  ! =================================================================================
  ! Read command line arguments
  narg = command_argument_count()
  if (narg < 1) then
     numtest = 1
  else
     call get_command_argument(1, arg)
     read (arg,*) numtest
  end if

  PRINT *,'**********************************************'
  PRINT *,'NUMTEST =',NUMTEST
  PRINT *,'**********************************************'

  write (strnum2,'(I2.2)') numtest
  ! =================================================================================

  
  ! =================================================================================
  ! Import surfaces
  nsurf = 2
  allocate( surf(nsurf) )
  do isurf = 1,2
      write (strnum,'(I1)') isurf
      call read_polynomial( &
           surf(isurf)%x, &
           'coeffstest/C' // strnum // '_test' // strnum2 // '.txt', &
           nvar=2, &
           base=1 )

      if ( ECONOMIZE ) call economize2( surf(isurf)%x, EPSmath )

      call compute_deriv1( surf(isurf) )
      call compute_deriv2( surf(isurf) )
      call compute_pseudonormal( surf(isurf) )
      if ( ECONOMIZE ) call economize2( surf(isurf)%pn, EPSmath )


      CALL WRITE_POLYNOMIAL( surf(isurf)%x,  'dev_intersection/surfroot' // strnum // '_x.cheb'  )
      CALL WRITE_POLYNOMIAL( surf(isurf)%pn, 'dev_intersection/surfroot' // strnum // '_pn.cheb' )
  end do
  ! =================================================================================



  
  ! =================================================================================
  ! Compute all intersection
  allocate(mask(nsurf))
  mask(:) = .true.
  call intersect_all_surfaces( &
       surf, &
       nsurf, &
       interdata_global, &
       mask )
  ! =================================================================================



  ! =================================================================================
  ! Free memory
  deallocate(surf, mask)
  ! =================================================================================
  



contains
!end program dev_intersection





  subroutine write_interdat( &
       interdat, &
       filepoints, &
       filecurves )
    use mod_util
    use mod_types_intersection
    implicit none
    type(type_intersection_data), intent(in) :: interdat
    character(*),                 intent(in) :: filepoints, filecurves
    integer                                  :: fileunit
    integer                                  :: ip, ic

    call get_free_unit( fileunit )

    open( &
         unit = fileunit, &
         file = filepoints, &
         action = 'write' )
    do ip = 1,interdat%np
       write ( fileunit, * ) interdat%points(ip)%xyz!, interdat%points(ip)%pos%uv, interdat%points(ip)%pos%next%uv
    end do
    close( fileunit )

    open( &
         unit = fileunit, &
         file = filecurves, &
         action = 'write' )
    write ( fileunit, * ) interdat%nc
    do ic = 1,interdat%nc
       write ( fileunit, * ) interdat%curves(ic)%root%endpoints
       write ( fileunit, * ) interdat%curves(ic)%uvbox(:,:,1)
       write ( fileunit, * ) interdat%curves(ic)%uvbox(:,:,2)
       if ( associated(interdat%curves(ic)%polyline) ) then
          write ( fileunit, * ) interdat%curves(ic)%polyline%np
          do ip = 1,interdat%curves(ic)%polyline%np
             write ( fileunit, * ) interdat%curves(ic)%polyline%uv(:,:,ip), interdat%curves(ic)%polyline%xyz(:,ip)
          end do
       else
          write ( fileunit, * ) 0
       end if
    end do
    close( fileunit )

  end subroutine write_interdat

