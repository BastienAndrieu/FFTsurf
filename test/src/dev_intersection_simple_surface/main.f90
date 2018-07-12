!program dev_intersection_simple_surface

  use mod_util
  use mod_math
  use mod_regiontree
  use mod_polynomial
  use mod_diffgeom2
  use mod_separation  

  implicit none

  type type_point_on_surface
     type(type_surface), pointer                :: surf => null()
     real(kind=fp)                              :: uv(2)
     type(type_point_on_surface), pointer       :: next => null()
  end type type_point_on_surface


  type type_intersection_point
     real(kind=fp)                              :: xyz(3)
     type(type_point_on_surface), pointer       :: pos => null()
     integer                                    :: npos = 0
  end type type_intersection_point


  type type_intersection_segment
     integer                                    :: endpoints(2)
     type(type_intersection_segment), pointer   :: parent => null()
     type(type_intersection_segment), pointer   :: child(:) => null()
  end type type_intersection_segment

  type type_intersection_polyline                                                                                               !
     integer                                    :: np = 0    ! actual nb. of poyline points                                     !
     real(kind=fp), allocatable                 :: s(:)      ! approximate curvilinear abscissa                                 !
     real(kind=fp), allocatable                 :: uv(:,:,:) ! uv coords.  (u/v, #surface, #point)                              !
     real(kind=fp), allocatable                 :: xyz(:,:)  ! xyz coords. (x/y/z, #point)                                      !
  end type type_intersection_polyline                                                                                           !

  type type_intersection_curve
     logical                                    :: smooth = .false.
     type(ptr_surface)                          :: surf(2)
     !type(ptr_region)                           :: region(2)
     real(kind=fp)                              :: uvbox(2,2,2) ! umin/max, vmin/max, #surf
     real(kind=fp)                              :: param_vector(3)
     type(type_intersection_segment)            :: root
     type(type_intersection_polyline)           :: polyline
  end type type_intersection_curve


  type type_intersection_data
     integer                                    :: np = 0
     type(type_intersection_point), allocatable :: points(:)
     integer                                    :: nc = 0
     type(type_intersection_curve), allocatable :: curves(:)
  end type type_intersection_data





  type type_proto_curve
     real(kind=fp) :: uvxyz(7,2)
  end type type_proto_curve
  

  type type_listcurves
     integer                       :: nc = 0
     type(type_proto_curve) :: curve(100)
  end type type_listcurves

  ! =================================================================================
  logical, parameter         :: ECONOMIZE = .true.

  integer                    :: narg, numtest
  character(100)             :: arg
  character                  :: strnum
  character(2)               :: strnum2
  integer                    :: fileunit

  type(type_surface), target :: surf(2)
  type(type_region), target  :: root(2)
  type(ptr_surface)          :: surfroot(2)
  type(ptr_region)           :: region(2)
  !type(type_listcurves)      :: listcurv
  type(type_intersection_data) :: interdat

  real(kind=fp), allocatable :: uvxyz(:,:)
  integer                    :: nuvxyz
  integer                    :: stat_degeneracy
  integer                    :: isurf
  integer*8                  :: tic, toc, count_rate, i
  ! =================================================================================


  ! Lecture arguments
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
  do isurf = 1,2
     write (strnum,'(I1)') isurf
     call read_polynomial( &
          surf(isurf)%x, &
          'coeffstest/C' // strnum // '_test' // strnum2 // '.txt', &
          nvar=2, &
          base=1 )

     if ( ECONOMIZE ) then
        !PRINT *,'     S DEGR =',SURF(ISURF)%S%DEGR
        call economize2( surf(isurf)%x, EPSmath )
        !PRINT *,' S_eco DEGR =',SURF(ISURF)%S%DEGR
     end if
     call write_polynomial( surf(isurf)%x, 'dev_intersection_simple_surface/c_'//strnum//'.cheb' )

     call compute_deriv1( surf(isurf) )
     call compute_deriv2( surf(isurf) )

     call write_polynomial( surf(isurf)%xu, 'dev_intersection_simple_surface/du_'//strnum//'.cheb' )
     call write_polynomial( surf(isurf)%xv, 'dev_intersection_simple_surface/dv_'//strnum//'.cheb' )

     call init_region( &
          root(isurf), &
          2, &
          [( [-1._fp, 1._fp], i=1,2 )] ) 
  end do
  ! =================================================================================




  allocate( ch2be_matrices( &
       max( maxval(surf(1)%x%degr), maxval(surf(2)%x%degr) ) + 1 &
       ) )

  do isurf = 1,2
     region(isurf)%ptr => root(isurf)
     surfroot(isurf)%ptr => surf(isurf)

     allocate( region(isurf)%ptr%poly(1) )
     allocate( region(isurf)%ptr%poly(1)%ptr )
     call cheb2bern_poly( &
          surf(isurf)%x, &
          region(isurf)%ptr%poly(1)%ptr )
  end do


  allocate( interdat%curves(100) )
  nuvxyz = 0
  stat_degeneracy = 0
  call system_clock( tic, count_rate )
  call intersect_simple_surfaces( &
       surfroot, &
       region, &
       [1._fp, 0._fp, 0._fp], &
       interdat, &
       !listcurv, &
       uvxyz, &
       nuvxyz, &
       stat_degeneracy )
  call system_clock( toc )
  PRINT *,''; PRINT *,''; PRINT *,''
  PRINT *,'ELAPSED =',REAL( TOC - TIC ) / REAL( COUNT_RATE )

  PRINT *,'STAT_DEGENERACY =',STAT_DEGENERACY

  !call get_free_unit( fileunit )
  !open( &
  !     unit=fileunit, &
  !     file='dev_intersection_simple_surface/curves.dat', &
  !     action='write' )
  !write ( fileunit, * ) listcurv%nc
  !do i = 1,listcurv%nc
  !   write ( fileunit, * ) listcurv%curve(i)%uvxyz(:,1)
  !   write ( fileunit, * ) listcurv%curve(i)%uvxyz(:,2)
  !end do
  !close( fileunit )




  if ( allocated( ch2be_matrices ) ) then
     do i = 1,size(ch2be_matrices)
        if ( allocated(ch2be_matrices(i)%mat) ) deallocate( ch2be_matrices(i)%mat )
     end do
     deallocate( ch2be_matrices )
  end if

  if ( allocated(uvxyz) ) deallocate( uvxyz )



  PRINT *,'-------------- INTERSECTION DATA --------------'
  PRINT *,INTERDAT%NP, ' POINTS'
  !DO I = 1,INTERDAT%NP
  !   CALL print_intersection_point( interdat%points(i) )
  !END DO
  
  
  call get_free_unit( fileunit )
  open( &
       unit=fileunit, &
       file='dev_intersection_simple_surface/intersection_points.dat', &
       action='write' )
  write ( fileunit, * ) interdat%np
  do i = 1,interdat%np
     !write ( fileunit, * ) interdat%points(i)%xyz
     write ( fileunit, * ) interdat%points(i)%xyz, interdat%points(i)%pos%uv, interdat%points(i)%pos%next%uv
  end do
  close( fileunit )

  PRINT *,INTERDAT%NC, ' CURVES'
  open( &
       unit=fileunit, &
       file='dev_intersection_simple_surface/intersection_curves.dat', &
       action='write' )
  write ( fileunit, * ) interdat%nc
  do i = 1,interdat%nc
     write ( fileunit, * ) interdat%curves(i)%root%endpoints
     !write ( fileunit, * ) interdat%curves(i)%region(1)%ptr%uvbox
     !write ( fileunit, * ) interdat%curves(i)%region(2)%ptr%uvbox
     write ( fileunit, * ) interdat%curves(i)%uvbox(:,:,1)
     write ( fileunit, * ) interdat%curves(i)%uvbox(:,:,2)
  end do
  close( fileunit )



  
  do isurf = 1,2
     write (strnum,'(I1)') isurf
     call export_region_tree( &
          region(isurf)%ptr, &
          'dev_intersection_simple_surface/tree_' // strnum // '.dat' )

     call free_polynomial( region(isurf)%ptr%poly(1)%ptr )
     deallocate( region(isurf)%ptr%poly(1)%ptr )
     deallocate( region(isurf)%ptr%poly )
     call free_region_tree( region(isurf)%ptr )

     call free_polynomial( surf(isurf)%x )
     call free_polynomial( surf(isurf)%xu )
     call free_polynomial( surf(isurf)%xv )
     call free_polynomial( surf(isurf)%xuu )
     call free_polynomial( surf(isurf)%xuv )
     call free_polynomial( surf(isurf)%xvv )

     !deallocate( root(isurf)%uvbox )
  end do

contains




  subroutine cheb2bern_poly( c, b )
    implicit none
    type(type_polynomial), intent(in)  :: c
    type(type_polynomial), intent(out) :: b
    real(kind=fp)                      :: au(c%degr(1)+1,c%degr(1)+1)
    real(kind=fp)                      :: av(c%degr(2)+1,c%degr(2)+1)
    integer                            :: i

    if ( c%base /= 1 ) STOP 'cheb2bern_poly : input polynomial not in Chebyshev basis'

    !PRINT *,'FOO'
    !PRINT *,'NVAR=',C%NVAR,', DEGR=',
    call reset_polynomial( poly=b, nvar=c%nvar, base=2, degr=c%degr(1:c%nvar), dim=c%dim )
    !PRINT *,'BAR'

    call get_cheb2bern_mat_from_collection( &
         ch2be_matrices, &
         c%degr(1)+1, &
         au )

    select case (c%nvar)
    case (1)
       b%coef(1:c%degr(1)+1,1:b%dim,1) = matmul( au, c%coef(1:c%degr(1)+1,1:c%dim,1) )

    case (2)
       call get_cheb2bern_mat_from_collection( &
            ch2be_matrices, &
            c%degr(2)+1, &
            av )
       av = transpose( av )

       do i = 1,c%dim
          b%coef(1:c%degr(1)+1,1:c%degr(2)+1,i) = matmul( &
               matmul( au, c%coef(1:c%degr(1)+1,1:c%degr(2)+1,i) ), &
               av )
       end do

    case default
       STOP 'cheb2bern_poly : nvar /= 1,2'

    end select

  end subroutine cheb2bern_poly

!end program dev_intersection_simple_surface
