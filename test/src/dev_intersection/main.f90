!program dev_intersection

  use mod_util
  use mod_math
  use mod_polynomial
  use mod_diffgeom2  
  use mod_types_intersection

  implicit none

  LOGICAL, PARAMETER :: GLOBALDEBUG = .true.
  
  integer                    :: narg
  character(100)             :: arg
  integer, allocatable       :: numtest(:), numsurf(:)
  !integer                    :: numtest
  character                  :: strnum
  character(2)               :: strnum2
  integer*8                  :: tic, toc, count_rate

  logical, parameter                      :: ECONOMIZE = .true.
  integer                                 :: nsurf
  type(type_surface), target, allocatable :: surf(:)
  logical, allocatable                    :: mask(:)
  type(type_intersection_data)            :: interdata_global
  integer                                 :: isurf

  ! =================================================================================
  ! Read command line arguments
  narg = command_argument_count()
  if ( narg < 1 ) then
     nsurf = 2
     allocate( numtest(nsurf), numsurf(nsurf) )
     numtest(:) = 1
     numsurf = [1,2]
  elseif ( mod(narg,2) /= 0 ) then
     nsurf = 2
     allocate( numtest(nsurf), numsurf(nsurf) )
     call get_command_argument(1, arg)
     read (arg,*) numtest(1)
     numtest(2) = numtest(1)
     numsurf = [1,2]
  else
     nsurf = narg / 2
     allocate( numtest(nsurf), numsurf(nsurf) )
     do isurf = 1,nsurf
        call get_command_argument(2*isurf-1, arg)
        read (arg,*) numtest(isurf)
        call get_command_argument(2*isurf, arg)
        read (arg,*) numsurf(isurf)
     end do
  end if


  PRINT *,'NUMTEST =',NUMTEST
  PRINT *,'NUMSURF =',NUMSURF

  !narg = command_argument_count()
  !if ( narg < 1 ) then
  !   numtest = 1
  !else
  !   call get_command_argument(1, arg)
  !   read (arg,*) numtest
  !end if
  !PRINT *,'**********************************************'
  !PRINT *,'NUMTEST =',NUMTEST
  !PRINT *,'**********************************************
  !write (strnum2,'(I2.2)') numtest
  ! =================================================================================

  
  ! =================================================================================
  ! Import surfaces
  allocate( surf(nsurf) )
  !do isurf = 1,2
  do isurf = 1,nsurf
     write (strnum2,'(I2.2)') numtest(isurf)
     write (strnum,'(I1)') numsurf(isurf)
     !write (strnum,'(I1)') isurf
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

     write (strnum,'(I1)') isurf
     CALL WRITE_POLYNOMIAL( surf(isurf)%x,  'dev_intersection/surfroot' // strnum // '_x.cheb'  )
     CALL WRITE_POLYNOMIAL( surf(isurf)%pn, 'dev_intersection/surfroot' // strnum // '_pn.cheb' )
  end do
  ! =================================================================================



  
  ! =================================================================================
  ! Compute all intersection
  allocate(mask(nsurf))
  mask(:) = .true.
  call system_clock( tic, count_rate )
  call intersect_all_surfaces( &
       surf, &
       nsurf, &
       interdata_global, &
       mask )
  call system_clock( toc )
  PRINT *,''; PRINT *,''; PRINT *,''
  PRINT *,'ELAPSED =',REAL( TOC - TIC ) / REAL( COUNT_RATE )
  ! =================================================================================



  ! =================================================================================
  ! Free memory
  deallocate(numtest, numsurf)
  deallocate(surf, mask)
  call free_transformation_matrices()
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

