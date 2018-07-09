!program dev_intersection_surface_surface

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
  !integer                    :: fileunit

  type(type_surface), target :: surf(2)
  type(type_region), target  :: root(2)
  type(ptr_surface)          :: surfroot(2)
  type(ptr_region)           :: region(2)
  type(type_intersection_data) :: interdat

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

     if ( ECONOMIZE ) call economize2( surf(isurf)%x, EPSmath )
     call write_polynomial( surf(isurf)%x, 'dev_intersection_surface_surface/c_'//strnum//'.cheb' )

     call compute_deriv1( surf(isurf) )
     call compute_deriv2( surf(isurf) )
     call compute_pseudonormal( surf(isurf) )
     if ( ECONOMIZE ) call economize2( surf(isurf)%pn, EPSmath )

     call write_polynomial( surf(isurf)%xu, 'dev_intersection_surface_surface/du_'//strnum//'.cheb' )
     call write_polynomial( surf(isurf)%xv, 'dev_intersection_surface_surface/dv_'//strnum//'.cheb' )
     call write_polynomial( surf(isurf)%pn, 'dev_intersection_surface_surface/pn_'//strnum//'.cheb' )

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

     allocate( region(isurf)%ptr%poly(2) )
     allocate( region(isurf)%ptr%poly(1)%ptr, region(isurf)%ptr%poly(2)%ptr )
     call cheb2bern_poly( &
          surf(isurf)%x, &
          region(isurf)%ptr%poly(1)%ptr )
     call cheb2bern_poly( &
          surf(isurf)%pn, &
          region(isurf)%ptr%poly(2)%ptr )
  end do

  allocate( interdat%curves(200) )
  stat_degeneracy = 0
  call system_clock( tic, count_rate )
  call intersect_surface_surface( &
       surfroot, &
       region, &
       interdat, &
       stat_degeneracy )!, &
       !.false., &
       ![0._fp, 0._fp, 0._fp] ) 
  call system_clock( toc )
  PRINT *,''; PRINT *,''; PRINT *,''
  PRINT *,'ELAPSED =',REAL( TOC - TIC ) / REAL( COUNT_RATE )

  PRINT *,'STAT_DEGENERACY =',STAT_DEGENERACY

  PRINT *,INTERDAT%NP,' INTERSECTION POINT(S)'
  PRINT *,INTERDAT%NC,' INTERSECTION CURVE(S)'

  call write_interdat( &
       interdat, &
       'dev_intersection_surface_surface/interdat_points.dat', &
       'dev_intersection_surface_surface/interdat_curves.dat' )

  call free_intersection_data( &
       interdat )
  
  if ( allocated( ch2be_matrices ) ) then
     do i = 1,size(ch2be_matrices)
        if ( allocated(ch2be_matrices(i)%mat) ) deallocate( ch2be_matrices(i)%mat )
     end do
     deallocate( ch2be_matrices )
  end if

  
  do isurf = 1,2
     write (strnum,'(I1)') isurf
     call export_region_tree( &
          region(isurf)%ptr, &
          'dev_intersection_surface_surface/tree_' // strnum // '.dat' )

     IF (.true.) THEN
        call free_polynomial( region(isurf)%ptr%poly(1)%ptr )
        call free_polynomial( region(isurf)%ptr%poly(2)%ptr )
        deallocate( region(isurf)%ptr%poly(1)%ptr, region(isurf)%ptr%poly(2)%ptr )
        deallocate( region(isurf)%ptr%poly )
        call free_region_tree( region(isurf)%ptr )

        call free_polynomial( surf(isurf)%x )
        call free_polynomial( surf(isurf)%xu )
        call free_polynomial( surf(isurf)%xv )
        call free_polynomial( surf(isurf)%xuu )
        call free_polynomial( surf(isurf)%xuv )
        call free_polynomial( surf(isurf)%xvv )

        !deallocate( root(isurf)%uvbox )
     END IF
  end do




contains




  subroutine cheb2bern_poly( c, b )
    implicit none
    type(type_polynomial), intent(in)    :: c
    type(type_polynomial), intent(inout) :: b
    real(kind=fp)                        :: au(c%degr(1)+1,c%degr(1)+1)
    real(kind=fp)                        :: av(c%degr(2)+1,c%degr(2)+1)
    integer                              :: i

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





  subroutine write_interdat( &
       interdat, &
       filepoints, &
       filecurves )
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
       write ( fileunit, * ) interdat%curves(ic)%polyline%np
       do ip = 1,interdat%curves(ic)%polyline%np
          write ( fileunit, * ) interdat%curves(ic)%polyline%uv(:,:,ip), interdat%curves(ic)%polyline%xyz(:,ip)
       end do
    end do
    close( fileunit )

  end subroutine write_interdat

!end program dev_intersection_surface_surface
