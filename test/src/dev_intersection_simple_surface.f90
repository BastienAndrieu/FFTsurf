program dev_intersection_simple_surface

  use mod_math
  use mod_bernstein
  use mod_diffgeom
  use mod_obb

  implicit none

  real(kind=MATHpr), parameter             :: EPSregion = real( 1.e-8, kind=MATHpr )

  type type_curve_region
     real(kind=MATHpr)                     :: tbox(2)
     type(type_obb), pointer               :: xyzbox => null()
     integer                               :: npts = 0
     integer, allocatable                  :: ipts(:)
     type(type_curve_region), pointer      :: parent => null()
     type(type_curve_region), pointer      :: child(:) => null()
     type(type_bernstein_series1), pointer :: bezier => null()
  end type type_curve_region


  type type_surface_region
     real(kind=MATHpr)                     :: uvbox(2,2)
     type(type_obb), pointer               :: xyzbox => null()
     integer                               :: npts = 0
     integer, allocatable                  :: ipts(:)
     type(type_surface_region), pointer    :: parent => null()
     type(type_surface_region), pointer    :: child(:) => null()
     type(type_bernstein_series2), pointer :: bezier => null()
  end type type_surface_region


  type ptr_surface_region
     type(type_surface_region), pointer    :: ptr => null()
  end type ptr_surface_region



  !type type_intersection_data                                                                                                   !
  !   integer                                    :: np = 0    ! actual nb. of intersection_points                                !
  !   integer                                    :: nc = 0    ! actual nb. of intersection_curves                                !
  !   type(type_intersection_point), allocatable :: points(:) ! list of intersection_points                                      !
  !   type(type_intersection_curve), allocatable :: curves(:) ! list of intersection_curves                                      !
  !end type type_intersection_data                                                                                               !


  ! =============================================================================

  logical, parameter             :: ECONOMIZE = .true.

  integer                        :: narg, numtest, icurv, ivar, ival
  character(100)                 :: arg
  character                      :: strnum
  character(2)                   :: strnum2

  type(type_parametric_surface), target :: surf(2)
  type(type_surface_region), target     :: root(2)
  type(ptr_surface_region)              :: region(2)
  type(ptr_parametric_surface)          :: surfroot(2)

  !type(type_parametric_curve)    :: curv
  !type(type_curve_region)        :: region_c
  !type(type_surface_region)      :: region_s
  real(kind=MATHpr), allocatable :: uvxyz(:,:)
  integer                        :: nuvxyz
  integer                        :: stat_degeneracy
  integer                        :: isurf, ipt
  integer*8                      :: tic, toc, count_rate

  !PRINT *,MAX(5,CEILING(-LOG10( 1.D-9 )))
  !STOP


  ! =================================================================================
  ! Lecture argument (nom du fichier configuration)
  narg = command_argument_count()
  if (narg < 1) then
     numtest = 1
  else
     call get_command_argument(1, arg)
     read (arg,*) numtest
  end if

  if (narg < 2) then
     icurv = 1
  else
     call get_command_argument(2, arg)
     read (arg,*) icurv
  end if

  if (narg < 3) then
     ivar = 1
  else
     call get_command_argument(3, arg)
     read (arg,*) ivar
  end if

  if (narg < 4) then
     ival = 1
  else
     call get_command_argument(4, arg)
     read (arg,*) ival
  end if

  PRINT *,'**********************************************'
  PRINT *,'NUMTEST =',NUMTEST
  !PRINT *,'  ICURV =',ICURV
  !PRINT *,'   IVAR =',IVAR
  !PRINT *,'   IVAL =',IVAL
  PRINT *,'**********************************************'

  write (strnum2,'(I2.2)') numtest
  ! =================================================================================


  ! =================================================================================
  do isurf = 1,2
     write (strnum,'(I1)') isurf
     call read_chebyshev_series2( &
          surf(isurf)%s, &
          '/stck/bandrieu/Bureau/coeffstest/C' // strnum // '_test' // strnum2 // '.txt' )
     !'/home/bastien/Bureau/coeffstest/C' // strnum // '_test' // strnum2 // '.txt' )

     if ( ECONOMIZE ) then
        !PRINT *,'     S DEGR =',SURF(ISURF)%S%DEGR
        call economization2( surf(isurf)%s, MATHeps )
        !PRINT *,' S_eco DEGR =',SURF(ISURF)%S%DEGR
     end if
     call write_chebyshev_series2( surf(isurf)%s, 'dev_intersection_simple_surface/c_'//strnum//'.cheb' )

     call compute_first_derivatives( surf(isurf) )
     call compute_second_derivatives( surf(isurf) )

     call init_surface_region( &
          root(isurf), &
          spread( real( [-1,1], kind=MATHpr ), 2, 2 ) )
  end do
  ! =================================================================================



  do isurf = 1,2
     region(isurf)%ptr => root(isurf)
     surfroot(isurf)%ptr => surf(isurf)

     allocate( region(isurf)%ptr%bezier )
     call cs2bs2( &
          surf(isurf)%s, &
          region(isurf)%ptr%bezier )
  end do




  allocate( uvxyz(7,10) )
  nuvxyz = 0
  stat_degeneracy = 0
  call system_clock( tic, count_rate )
  outer : do icurv = 1,2
     do ivar = 1,2
        do ival = 1,2
           call intersect_border_surface( &
                surfroot, &
                region, &
                icurv, &
                ivar, &
                ival, &
                uvxyz, &
                nuvxyz, &
                stat_degeneracy )
           if ( stat_degeneracy > 1 ) exit outer
        end do
     end do
  end do outer
  call system_clock( toc )
  PRINT *,'ELAPSED =',REAL( TOC - TIC ) / REAL( COUNT_RATE )

  open( unit=13, file='dev_intersection_simple_surface/uv_xyz.dat', action='write' )
  if ( nuvxyz < 1 ) then
     write (13,*) ''
  else
     do ipt = 1,nuvxyz
        write (13,*) uvxyz(:,ipt)
     end do
  end if
  close(13)
  
  DO ISURF = 1,2
     DEALLOCATE( REGION(ISURF)%PTR%BEZIER%COEF )
     DEALLOCATE( REGION(ISURF)%PTR%BEZIER )
  END DO

  

contains

  subroutine cs2bs1( c, b )
    implicit none
    type(type_chebyshev_series1), intent(in)  :: c
    type(type_bernstein_series1), intent(out) :: b
    real(kind=MATHpr)                         :: a(c%degr+1,c%degr+1)

    call reset_bernstein_series1( b, c%degr, c%dim )

    call ch2be_matrix( a, c%degr )

    b%coef(1:c%degr+1,:) = matmul( a, c%coef(1:c%degr+1,:) )

  end subroutine cs2bs1





  subroutine cs2bs2( c, b )
    implicit none
    type(type_chebyshev_series2), intent(in)  :: c
    type(type_bernstein_series2), intent(out) :: b
    real(kind=MATHpr)                         :: au(c%degr(1)+1,c%degr(1)+1)
    real(kind=MATHpr)                         :: av(c%degr(2)+1,c%degr(2)+1)
    integer                                   :: i

    call reset_bernstein_series2( b, c%degr, c%dim )

    call ch2be_matrix( au, c%degr(1) )
    call ch2be_matrix( av, c%degr(2) )
    av = transpose( av )

    do i = 1,c%dim
       b%coef(1:c%degr(1)+1,1:c%degr(2)+1,i) = matmul( matmul( au, c%coef(1:c%degr(1)+1,1:c%degr(2)+1,i) ), av )
    end do

  end subroutine cs2bs2









  subroutine intersect_simple_surfaces( &
      surfroot, &
      region, &
      param_vector, &
      interdat )
    implicit none
    integer, parameter                          :: nuvxyz_init = 10
    type(ptr_parametric_surface), intent(in)    :: surfroot(2)
    type(ptr_surface_region),     intent(inout) :: region(2)
    real(kind=MATHpr),            intent(in)    :: param_vector(3)
    integer,                      intent(inout) :: interdat ! *********
    real(kind=MATHpr), allocatable              :: uvxyz(:,:)
    integer                                     :: nuvxyz
    integer, allocatable                        :: sharedpts(:)
    integer                                     :: isurf, icurv, ivar, ival, jpt, ipt
    
    do isurf = 1,2
       if ( .not.associated(region(isurf)%ptr%parent) ) cycle
       do jpt = 1,region(isurf)%ptr%parent%npts
          ipt = region(isurf)%ptr%parent%ipts(jpt)
          if ( &
               is_in_interval( uvxyz(2*isurf-1,ipt), region(isurf)%ptr%uvbox(1,1), region(isurf)%ptr%uvbox(2,1) ) .and. &
               is_in_interval( uvxyz(2*isurf,ipt), region(isurf)%ptr%uvbox(1,2), region(isurf)%ptr%uvbox(2,2) ) ) then
             !call append( &
             !     region(isurf)%ptr%ipts, &
             !     ipt, &
             !     noduplicates=.true., &
             !     newlength=region(isurf)%ptr%npts )
          end if
       end do
    end do

    allocate( uvxyz(7,nuvxyz_init) )
    nuvxyz = 0
    if ( region(1)%ptr%npts > 0 .and. region(2)%ptr%npts > 0 ) then
       call intersection_arrays( &
            region(1)%ptr%ipts(1:region(1)%ptr%npts), &
            region(2)%ptr%ipts(1:region(2)%ptr%npts), &
            sharedpts )
       if ( allocated(sharedpts) ) then
          nuvxyz = size(sharedpts)
          ! ***** 
       end if
    end if
    
    
    ! intersect the 4 borders of each surface with the other surface
    do icurv = 1,2
       do ivar = 1,2
          do ival = 1,2

             call intersect_border_surface( &
                  surfroot, &
                  region, &
                  icurv, &
                  ivar, &
                  ival, &
                  uvxyz, &
                  nuvxyz, &
                  stat_degeneracy )
             if ( stat_degeneracy > 1 ) then
                ! propagation erreur vers routine appelante ...
                return
             end if

          end do
       end do
    end do


    if ( nuvxyz < 1 ) return



  end subroutine intersect_simple_surfaces



























  subroutine intersect_border_surface( &
       root_s, &
       region, &
       icurv, &
       ivar, &
       ival, &
       uvxyz, &
       nuvxyz, &
       stat_degeneracy )
    use mod_math
    use mod_chebyshev
    use mod_bernstein
    use mod_diffgeom
    use mod_tolerances
    implicit none
    integer, parameter                            :: npts_init = 10
    type(ptr_parametric_surface),   intent(in)    :: root_s(2)
    type(ptr_surface_region),       intent(inout) :: region(2)
    integer,                        intent(in)    :: icurv, ivar, ival
    real(kind=MATHpr), allocatable, intent(inout) :: uvxyz(:,:)
    integer,                        intent(inout) :: nuvxyz
    integer,                        intent(inout) :: stat_degeneracy
    type(type_parametric_curve)                   :: root_c
    type(type_curve_region)                       :: region_c
    type(type_surface_region)                     :: region_s
    real(kind=MATHpr), allocatable                :: tuvxyz(:,:)
    integer                                       :: ntuvxyz
    real(kind=MATHpr)                             :: uv(2,2)
    integer                                       :: isurf, ipt, jpt

    !PRINT *,''; PRINT *,''; PRINT *,''; PRINT *,''
    !PRINT *,'ICURV, IVAR, IVAL =', ICURV, IVAR, IVAL

    isurf = 1 + mod(icurv,2)

    ! convert surface border to parametric curve
    call cs2edge2cs1( &
         root_s(icurv)%ptr%s, &
         ivar, &
         ival, &
         root_c%c )
    call economization1( root_c%c, MATHeps )
    call compute_first_derivatives( root_c )
    call compute_second_derivatives( root_c )


    ! initialize curve region tree
    call init_curve_region( &
         region_c, &
         real( [-1,1], kind=MATHpr ) )

    ! initialize surface region tree
    call init_surface_region( &
         region_s, &
         spread( real( [-1,1], kind=MATHpr ), 2, 2 ) )

    ! compute curve Bezier control points
    allocate( region_c%bezier )
    call cs2bs1( &
         root_c%c, &
         region_c%bezier )

    ! copy surface Bezier control points
    region_s%bezier => region(isurf)%ptr%bezier

    
    !CALL WRITE_BERNSTEIN_SERIES1( REGION_C%BEZIER, 'dev_intersection_simple_surface/root_c_bezier.bern' )
    !CALL WRITE_BERNSTEIN_SERIES2( REGION_S%BEZIER, 'dev_intersection_simple_surface/root_s_bezier.bern' )


    allocate( tuvxyz(6,npts_init) )
    ntuvxyz = 0
    stat_degeneracy = 0
    call intersect_curve_surface( &
         root_c, &
         root_s(isurf)%ptr, &
         region_c, &
         region_s, &
         tuvxyz, &
         ntuvxyz, &
         stat_degeneracy )
    
    !open( unit=13, file='dev_intersection_simple_surface/tuv_xyz.dat', action='write' )
    !if ( ntuvxyz < 1 ) then
    !   write (13,*) ''
    !else
    !   do ipt = 1,ntuvxyz
    !      write (13,*) tuvxyz(:,ipt)
    !   end do
    !end if
    !close(13)

    
    !call export_surface_region_tree( &
    !     region_s, &
    !     'dev_intersection_simple_surface/tree_s.dat' )
    !call export_curve_region_tree( &
    !     region_c, &
    !     'dev_intersection_simple_surface/tree_c.dat' )
    IF ( STAT_DEGENERACY > 10 ) THEN
       PRINT *,'STAT_DEGENERACY =',STAT_DEGENERACY
       RETURN
       !STOP '********************'
    END IF

    !PRINT  *,'';PRINT  *,'';PRINT  *,'';PRINT  *,'';PRINT  *,'';PRINT  *,'';
    outer : do ipt = 1,ntuvxyz
       do jpt = 1,nuvxyz
          if ( sum( ( tuvxyz(4:6,ipt) - uvxyz(5:7,jpt) )**2 ) < EPSxyzsqr ) cycle outer
       end do

       uv(ivar,icurv) = region(icurv)%ptr%uvbox(ival,ivar)
       uv(1+mod(ivar,2),icurv) = n1p12ab( &
            tuvxyz(1,ipt), &
            region(icurv)%ptr%uvbox(1,1+mod(ivar,2)), &
            region(icurv)%ptr%uvbox(2,1+mod(ivar,2)) )
       uv(:,isurf) = tuvxyz(2:3,ipt)
       !PRINT *,uv,tuvxyz(4:6,ipt)

       call append_vector( &
            [ uv(:,1), uv(:,2), tuvxyz(4:6,ipt) ], &
            7, &
            uvxyz, &
            nuvxyz )

    end do outer

    call free_surface_region_tree( region_s, no_dealloc_bezier=.true. )
    call free_curve_region_tree( region_c )

  end subroutine intersect_border_surface

























































































































  recursive subroutine intersect_curve_surface( &
       root_c, & ! 3a
       root_s, & ! 3a
       region_c, &
       region_s, &
       coords, &   ! 1d, 2b
       npts, &
       stat_newpoint )
    use mod_util    
    use mod_math
    use mod_diffgeom
    use mod_obb
    use mod_tolerances
    implicit none
    type(type_parametric_curve),    intent(in)            :: root_c
    type(type_parametric_surface),  intent(in)            :: root_s
    type(type_curve_region),        intent(inout), target :: region_c
    type(type_surface_region),      intent(inout), target :: region_s
    real(kind=MATHpr), allocatable, intent(inout)         :: coords(:,:) ! t,u,v,x,y,z
    integer,                        intent(inout)         :: npts
    integer,                        intent(inout)         :: stat_newpoint
    integer, allocatable                                  :: sharedpts(:)
    integer                                               :: n_sharedpts
    logical                                               :: separable, overlap
    real(kind=MATHpr)                                     :: tuv_subdiv(3)
    logical                                               :: interior(2)
    real(kind=MATHpr)                                     :: dist, dist_ijk
    integer                                               :: ijkm(3)
    real(kind=MATHpr)                                     :: tuv(3), xyz(3)
    integer                                               :: stat_subdiv, nchild(2)
    type(type_curve_region), pointer                      :: newregion_c
    type(type_surface_region), pointer                    :: newregion_s
    integer                                               :: ipt, jpt, i, j, k, ichild, jchild

    if ( stat_degeneracy > 1 ) return

    !PRINT *,''
    !PRINT *,''
    !PRINT *,' TBOX =',REGION_C%TBOX
    !PRINT *,'UVBOX =',REGION_S%UVBOX
    !IF ( REGION_C%TBOX(1) < -0.99D0 .AND. REGION_C%TBOX(2) > 0.02D0 .AND. &
    !     REGION_S%UVBOX(1,1) > 0.50D0 .AND. REGION_S%UVBOX(2,1) > 0.99D0 .AND. &
    !     REGION_S%UVBOX(1,2) < -0.99D0 .AND. REGION_S%UVBOX(2,2) > 0.55D0 ) THEN
    !   PRINT *,''
    !   PRINT *,''
    !   PRINT *,''
    !   PRINT *,''
    !   PRINT *,' TBOX =',REGION_C%TBOX
    !   PRINT *,'UVBOX =',REGION_S%UVBOX
    !END IF

    !!IF ( REGION_C%TBOX(2) - REGION_C%TBOX(1) < 1.D-3 ) RETURN !EPSREGION ) RETURN
    !!IF ( REGION_S%UVBOX(2,1) - REGION_S%UVBOX(1,1) < 1.D-3 ) RETURN !EPSREGION ) RETURN
    !!IF ( REGION_S%UVBOX(2,2) - REGION_S%UVBOX(1,2) < 1.D-3 ) RETURN !EPSREGION ) RETURN


    !
    if ( associated(region_c%parent) ) then
       do jpt = 1,region_c%parent%npts
          ipt = region_c%parent%ipts(jpt)
          if ( is_in_interval( coords(1,ipt), region_c%tbox(1), region_c%tbox(2) ) ) then
             call append( &
                  region_c%ipts, &
                  ipt, &
                  noduplicates=.true., &
                  newlength=region_c%npts )
          end if
       end do
    end if


    if ( associated(region_s%parent) ) then
       do jpt = 1,region_c%parent%npts
          ipt = region_c%parent%ipts(jpt)
          if ( &
               is_in_interval( coords(2,ipt), region_s%uvbox(1,1), region_s%uvbox(2,1) ) .and. &
               is_in_interval( coords(3,ipt), region_s%uvbox(1,2), region_s%uvbox(2,2) ) ) then
             call append( &
                  region_s%ipts, &
                  ipt, &
                  noduplicates=.true., &
                  newlength=region_s%npts )
          end if
       end do
    end if
  




    ! get list of already discovered points common to both the current curve and surface regions
    n_sharedpts = 0
    if ( region_c%npts > 0 .and. region_s%npts > 0 ) then
       call intersection_arrays( &
            region_c%ipts(1:region_c%npts), &
            region_s%ipts(1:region_s%npts), &
            sharedpts )
       if ( allocated(sharedpts) ) n_sharedpts = size(sharedpts)
    end if
    
    !IF ( NPTS >= 2 .AND. N_SHAREDPTS > 0 ) THEN
    !   PRINT *,'';PRINT *,'';PRINT *,''
    !   PRINT *,' TBOX =',REGION_C%TBOX
    !   PRINT *,'UVBOX =',REGION_S%UVBOX
    !   PRINT *,N_SHAREDPTS,' SHARED POINT(S)'
    !   IF ( N_SHAREDPTS > 0 ) THEN
    !      DO IPT = 1,N_SHAREDPTS
    !         PRINT *,COORDS(:,SHAREDPTS(IPT))
    !      END DO
    !   END IF
    !END IF

    if ( n_sharedpts == 1 ) then
       ! check if the curve and surface regions can intersect at other (not yet discovered) points
       call intersect_curve_surface_elsewhere( &
            region_c%bezier, &
            region_s%bezier, &
            coords(4:6,sharedpts(1)), &
            separable, &
            randomize=.true. )
       !PRINT *,'MAY INTERSECT AT OTHER POINTS?', .NOT.SEPARABLE
       
       !IF ( REGION_C%TBOX(1) > 0._MATHPR .AND. &
       !     REGION_S%UVBOX(1,1) < 0._MATHPR .AND. &
       !     REGION_S%UVBOX(2,1) < 0._MATHPR .AND. &
       !     REGION_S%UVBOX(1,2) > 0._MATHPR .AND. &
       !     .NOT.SEPARABLE ) THEN
       !IF (.FALSE.) THEN!( .NOT.SEPARABLE .AND. NPTS == 2 ) THEN
       !   PRINT *,' TBOX =',REGION_C%TBOX
       !   PRINT *,'UVBOX =',REGION_S%UVBOX
       !   CALL WRITE_BERNSTEIN_SERIES1( REGION_C%BEZIER, 'dev_intersection_simple_surface/region_c_bezier.bern' )
       !   CALL WRITE_BERNSTEIN_SERIES2( REGION_S%BEZIER, 'dev_intersection_simple_surface/region_s_bezier.bern' )
       !   OPEN(13,FILE='dev_intersection_simple_surface/debug_separation_xyz.dat', ACTION='WRITE')
       !   WRITE(13,*) COORDS(4:6,SHAREDPTS(1))
       !   CLOSE(13)
       !   PRINT *,'---> DEBUG SEPARATION'
       !   !if ( separable ) STOP
       !   !STOP
       !   STAT_NEWPOINT = 66
       !   RETURN
       !END IF


       !IF ( REGION_C%TBOX(1) < -0.99D0 .AND. REGION_C%TBOX(2) > 0.02D0 .AND. &
       !     REGION_S%UVBOX(1,1) > 0.50D0 .AND. REGION_S%UVBOX(2,1) > 0.99D0 .AND. &
       !     REGION_S%UVBOX(1,2) < -0.99D0 .AND. REGION_S%UVBOX(2,2) > 0.55D0 ) THEN
       !   PRINT *,'MAY INTERSECT AT OTHER POINTS?', .NOT.SEPARABLE
       !   PRINT *,''
       !   PRINT *,''
       !   PRINT *,''
       !END IF

       if ( separable ) return
       !PRINT *,'*** MAY INTERSECT AT OTHER POINTS ***'
       
       
    end if

    if ( n_sharedpts > 0 ) then
       ! are any of the already discovered points interior to the curve/surface region
       do jpt = 1,n_sharedpts
          ipt = sharedpts(jpt)
          interior(1) = is_in_interval_strict( coords(1,ipt), region_c%tbox(1), region_c%tbox(2) )
          interior(2) = ( &
               is_in_interval_strict( coords(2,ipt), region_s%uvbox(1,1), region_s%uvbox(2,1) ) .and. &
               is_in_interval_strict( coords(3,ipt), region_s%uvbox(1,2), region_s%uvbox(2,2) ) )
          if ( any(interior) ) then
             tuv_subdiv(1) = coords(1,ipt)
             tuv_subdiv(2:3) = coords(2:3,ipt)
             exit
          end if
       end do
       !PRINT *,'INTERIOR POINTS?',INTERIOR

    else
       interior(:) = .false.

       if ( .not.associated(region_c%xyzbox) ) then
          allocate( region_c%xyzbox )
          call bernOBB1( &
               region_c%bezier%coef(1:region_c%bezier%degr+1,:), &
               region_c%bezier%degr, &
               region_c%xyzbox )
       end if

       if ( .not.associated(region_s%xyzbox) ) then
          allocate( region_s%xyzbox )
          call bernOBB2( &
               region_s%bezier%coef(1:region_s%bezier%degr(1)+1,1:region_s%bezier%degr(2)+1,:), &
               region_s%bezier%degr, &
               region_s%xyzbox )
       end if

       call overlap_OBBs( &
            region_c%xyzbox, &
            region_s%xyzbox, &
            overlap )

       if ( .not.overlap ) then
          !PRINT *,'DISJOINT OBBs'
          return
       end if

    end if

    !PRINT *,'INTERIOR?',INTERIOR
    !PRINT *,'ALL(.NOT.INTERIOR) =',ALL(.NOT.INTERIOR)
    if ( all(.not.interior) ) then ! <------------------------------------------------------------------------------------------+
       !! Search for a new intersection point                                                                                   !
       if ( .FALSE. ) THEN !( n_sharedpts == 0 ) then ! <-----------------------------------------------------------------------------------+     !
          ! start from the tuv point corresponding to the closest pair of control points                                  !     !
          dist = huge(1._MATHpr)                                                                                          !     !
          do k = 1,region_s%bezier%degr(2)+1 ! <-----------------------------------------------------------+              !     !
             do j = 1,region_s%bezier%degr(1)+1 ! <------------------------------------------------------+ !              !     !
                do i = 1,region_c%bezier%degr+1 ! <----------------------------------------------------+ ! !              !     !
                   dist_ijk = &                                                                        ! ! !              !     !
                        ( region_c%bezier%coef(i,1) - region_s%bezier%coef(j,k,1) )**2 + &             ! ! !              !     !
                        ( region_c%bezier%coef(i,2) - region_s%bezier%coef(j,k,2) )**2 + &             ! ! !              !     !
                        ( region_c%bezier%coef(i,3) - region_s%bezier%coef(j,k,3) )**2                 ! ! !              !     !
                   if ( dist > dist_ijk ) then ! <---- +                                               ! ! !              !     !
                      dist = dist_ijk                  !                                               ! ! !              !     !
                      ijkm = [i,j,k]                   !                                               ! ! !              !     !
                   end if ! <--------------------------+                                               ! ! !              !     !
                end do ! <-----------------------------------------------------------------------------+ ! !              !     !
             end do ! <----------------------------------------------------------------------------------+ !              !     !
          end do ! <---------------------------------------------------------------------------------------+              !     !
          tuv(1) = real( ijkm(1)-1, kind=MATHpr ) / real( region_c%bezier%degr,    kind=MATHpr )                          !     !
          tuv(2) = real( ijkm(2)-1, kind=MATHpr ) / real( region_s%bezier%degr(1), kind=MATHpr )                          !     !
          tuv(3) = real( ijkm(3)-1, kind=MATHpr ) / real( region_s%bezier%degr(2), kind=MATHpr )                          !     !
          tuv = 2._MATHpr * tuv - 1._MATHpr                                                                               !     !
          tuv(1) = n1p12ab( tuv(1), region_c%tbox(1), region_c%tbox(2) )                                                  !     !
          tuv(2) = n1p12ab( tuv(2), region_s%uvbox(1,1), region_s%uvbox(2,1) )                                            !     !
          tuv(3) = n1p12ab( tuv(3), region_s%uvbox(1,2), region_s%uvbox(2,2) )                                            !     !
          !PRINT *,':::::::::: SHOULD NOT BE HERE ::::::::::'
       else ! ------------------------------------------------------------------------------------------------------------+     !
          tuv = 0.5_MATHpr * [ sum( region_c%tbox ), sum( region_s%uvbox, dim=1 ) ]                                       !     !
          !PRINT *,''; PRINT *,''
          !PRINT *,' TBOX=',REAL(REGION_C%TBOX)
          !PRINT *,'UVBOX=',REAL(REGION_S%UVBOX)
          !tuv(1) = 0.5_MATHpr * sum( region_c%tbox )
          !tuv(2:3) = 0.5_MATHpr * sum( region_s%uvbox, dim=1 )
          !PRINT *,'MIDPOINT =',REAL(TUV)
       end if ! <---------------------------------------------------------------------------------------------------------+     !
       !PRINT *,'NEWTON <--- TUV0 =',REAL(TUV)
       !PRINT *,'------- NEWTON -------'
       !PRINT *,'N_SHAREDPTS =',N_SHAREDPTS
       !PRINT *,'INTERIOR?    ',INTERIOR
       !IF ( N_SHAREDPTS < 1 ) THEN
       !   CALL WRITE_OBB( REGION_C%XYZBOX, 'dev_intersection_simple_surface/xyzbox_c.dat' )
       !   CALL WRITE_OBB( REGION_S%XYZBOX, 'dev_intersection_simple_surface/xyzbox_s.dat' )
       !   CALL WRITE_BERNSTEIN_SERIES1( REGION_C%BEZIER, 'dev_intersection_simple_surface/region_c_bezier.bern' )
       !   CALL WRITE_BERNSTEIN_SERIES2( REGION_S%BEZIER, 'dev_intersection_simple_surface/region_s_bezier.bern' )
       !   STOP '--> VERIF OBBs'
       !END IF
       call newton_curve_surface( &                                                                                             !
            root_c, &                                                                                                           !
            root_s, &                                                                                                           !
            region_c%tbox, &                                                                                                    !
            region_s%uvbox, &                                                                                                   !
            tuv, &                                                                                                              !
            stat_newpoint, &                                                                                                    !
            xyz )      
       !PRINT *,'----------------------'                                                                                      !
       !PRINT *,'STAT_NEWPOINT =',STAT_NEWPOINT                                                                                 !
       !
       ! if we just found a degenerate point, return and report degeneracy                                                      !
       if ( stat_newpoint > 1 ) return                                                                                          ! 
       !
       if ( stat_newpoint == 0 ) then ! <--------------------------------------------------------------------------------+      !
          !PRINT *,'TUVXYZ =',REAL(TUV),REAL(XYZ)                                                                        !      !
          !PRINT *,'COORDS ='                                                                                            !      !
          !CALL PRINT_MAT( REAL(COORDS(1:6,1:NPTS)) )                                                                    !      !
          ! check if the "new" point has not already been discovered                                                     !      !
          do ipt = 1,npts ! <-------------------------------------------------------------------+                        !      !
             if ( sum( ( xyz - coords(4:6,ipt) )**2 ) < EPSxyzsqr ) then ! <---------+          !                        !      ! 
                stat_newpoint = 1                                                    !          !                        !      !
                exit                                                                 !          !                        !      !
             end if ! <--------------------------------------------------------------+          !                        !      !
          end do ! <----------------------------------------------------------------------------+                        !      !
       end if ! <--------------------------------------------------------------------------------------------------------+      !
       ! 
       if ( stat_newpoint == 0 ) then ! <--------------------------------------------------------------------------------+      !
          ! if we just found a new intersection point, add it to the lists, and subdivide at that point                  !      !
          !PRINT *,'NEW POINT =',TUV, XYZ
          !PRINT *,'IN   TBOX =',REGION_C%TBOX
          !PRINT *,'IN  UVBOX =',REGION_S%UVBOX
          !PRINT *,'';PRINT *,'';PRINT *,'';
          
          !call append_tuvxyz_point( &                                                                                   !      !
          !     coords, &                                                                                                !      !
          !     npts, &                                                                                                  !      !
          !     [tuv,xyz] )                                                                                              !      !
          call append_vector( &                                                                                          !      !
               [tuv,xyz], &                                                                                              !      !
               6, &                                                                                                      !      !
               coords, &                                                                                                 !      !
               npts )                                                                                                    !      !
          !PRINT *,'COORDS ='
          !CALL PRINT_MAT( TRANSPOSE(COORDS(:,1:NPTS)) )
          call add_to_parent_curve( region_c, npts )                                                                     !      !
          call add_to_parent_surface( region_s, npts )                                                                   !      !
          !      !
          !PRINT *,'    REGION_C%IPTS =',REGION_C%IPTS                                                                   !      !
          !PRINT *,'    REGION_S%IPTS =',REGION_S%IPTS                                                                   !      !
          !      !
          tuv_subdiv = tuv                                                                                               !      !
          !IF (NPTS > 2) THEN
          !   STAT_NEWPOINT = 77
          !   RETURN
          !END IF
       else ! -----------------------------------------------------------------------------------------------------------+      !
          ! else, subdivide at parametric midpoint                                                                       !      !
          tuv_subdiv(1) = 0.5_MATHpr * sum( region_c%tbox )                                                              !      !
          tuv_subdiv(2:3) = 0.5_MATHpr * sum( region_s%uvbox, dim=1 )                                                    !      !
       end if ! <--------------------------------------------------------------------------------------------------------+      !
    else ! ---------------------------------------------------------------------------------------------------------------------+
       !PRINT *,'---> HERE <---'                                                                                                 !
       IF ( N_SHAREDPTS < 1 ) STOP 'INTERIOR BUT N_SHAREDPTS < 1'                                                               !
       tuv_subdiv = coords(1:3,sharedpts(1))                                                                                    !
       !PRINT *,'TUV_SUBDIV =',REAL( TUV_SUBDIV )                                                                                !
    end if !  <-----------------------------------------------------------------------------------------------------------------+



    !PRINT *,'TUV_SUBDIV =',REAL(TUV_SUBDIV)
    !! Subdivide the curve region
    call subdiv_curve_region( &
         region_c, &
         tuv_subdiv(1), &
         stat_subdiv )
    !PRINT *,'STAT_SUBDIV C =',STAT_SUBDIV
    tuv_subdiv(1) = 0.5_MATHpr * ( ab2n1p1( tuv_subdiv(1), region_c%tbox(1), region_c%tbox(2) ) + 1._MATHpr )
    !PRINT *,'STAT_SUBDIV (CURVE) =',STAT_SUBDIV
    !PRINT *,'SIZE(CHILD) =',SIZE(REGION_C%CHILD)

    if ( stat_subdiv > 0) then
       !PRINT *,'DO NOT SUBDIVIDE (CURVE ENDPOINT)'
       nchild(1) = 1
    else
       nchild(1) = size(region_c%child)
       IF ( NCHILD(1) /= 2 ) THEN
          PRINT *,'ERROR : NCHILD(1) =',NCHILD(1)
          STOP
       END IF

       !do ichild = 1,nchild(1)
       !   do jpt = 1,region_c%npts
       !      ipt = region_c%ipts(jpt)
       !      !PRINT *,REAL(COORDS(1,IPT)), ' IN ',REAL(REGION_C%CHILD(ICHILD)%TBOX), ' ? ', &
       !      !     is_in_interval( coords(1,ipt), region_c%child(ichild)%tbox(1), region_c%child(ichild)%tbox(2) )
       !      if ( is_in_interval( coords(1,ipt), region_c%child(ichild)%tbox(1), region_c%child(ichild)%tbox(2) ) ) then
       !         !IF ( REGION_C%CHILD(ICHILD)%NPTS > 0 ) PRINT *,'==>>>>> APPEND',IPT,'TO',REGION_C%CHILD(ICHILD)%IPTS
       !         call append( region_c%child(ichild)%ipts, ipt, noduplicates=.true., newlength=region_c%child(ichild)%npts )
       !         !region_c%child(ichild)%npts = region_c%child(ichild)%npts + 1
       !      end if
       !   end do
       !   !PRINT *,'CHILD #',ICHILD,', NPTS =', REGION_C%CHILD(ICHILD)%NPTS
       !end do

       if ( stat_subdiv == 0 ) then
          allocate( region_c%child(1)%bezier, region_c%child(2)%bezier )
          call subdiv1( &
               region_c%bezier, &
               tuv_subdiv(1), &
               bl=region_c%child(1)%bezier, &
               br=region_c%child(2)%bezier )
       end if

    end if


    !! Subdivide the surface region
    call subdiv_surface_region( &
         region_s, &
         tuv_subdiv(2:3), &
         stat_subdiv )
    !PRINT *,'STAT_SUBDIV S =',STAT_SUBDIV
    tuv_subdiv(2) = 0.5_MATHpr * ( ab2n1p1( tuv_subdiv(2), region_s%uvbox(1,1), region_s%uvbox(2,1) ) + 1._MATHpr )
    tuv_subdiv(3) = 0.5_MATHpr * ( ab2n1p1( tuv_subdiv(3), region_s%uvbox(1,2), region_s%uvbox(2,2) ) + 1._MATHpr )

    if ( stat_subdiv == 3 ) then
       !PRINT *,'DO NOT SUBDIVIDE (SURFACE CORNER)'
       nchild(2) = 1
       !stat_newpoint = 44
       !RETURN !STOP '??? SURFACE'
    else
       !PRINT *,'NCHILD_S =',NCHILD(2)
       nchild(2) = size(region_s%child)
       !do ichild = 1,nchild(2)
       !   do jpt = 1,region_s%npts
       !      ipt = region_s%ipts(jpt)
       !      if ( &
       !           is_in_interval( coords(2,ipt), region_s%child(ichild)%uvbox(1,1), region_s%child(ichild)%uvbox(2,1) ) .and. &
       !           is_in_interval( coords(3,ipt), region_s%child(ichild)%uvbox(1,2), region_s%child(ichild)%uvbox(2,2) ) ) then
       !         !IF ( REGION_S%CHILD(ICHILD)%NPTS > 0 ) PRINT *,'==>>>>> APPEND',IPT,'TO',REGION_S%CHILD(ICHILD)%IPTS
       !         call append( region_s%child(ichild)%ipts, ipt, noduplicates=.true., newlength=region_s%child(ichild)%npts )
       !         !region_s%child(ichild)%npts = region_s%child(ichild)%npts + 1
       !      end if
       !   end do
       !end do
    end if

    if ( stat_subdiv == 0 ) then
       allocate( &
            region_s%child(1)%bezier, region_s%child(2)%bezier, &
            region_s%child(3)%bezier, region_s%child(4)%bezier )
       call subdiv2( &
            region_s%bezier, &
            tuv_subdiv(2:3), &
            bsw=region_s%child(1)%bezier, &
            bse=region_s%child(2)%bezier, &
            bnw=region_s%child(3)%bezier, &
            bne=region_s%child(4)%bezier )
    elseif ( stat_subdiv == 1 ) then
       allocate( region_s%child(1)%bezier, region_s%child(2)%bezier )
       call subdiv2_along_v( &
            region_s%bezier, &
            tuv_subdiv(3), &
            bs=region_s%child(1)%bezier, &
            bn=region_s%child(2)%bezier )
    elseif ( stat_subdiv == 2 ) then
       allocate( region_s%child(1)%bezier, region_s%child(2)%bezier )
       call subdiv2_along_u( &
            region_s%bezier, &
            tuv_subdiv(2), &
            bw=region_s%child(1)%bezier, &
            be=region_s%child(2)%bezier )
    end if

    !end if



    !! Carry on the recursion with the children
    if ( all(nchild < 2) ) then  
       PRINT *,' TBOX =',REGION_C%TBOX
       PRINT *,'UVBOX =',REGION_S%UVBOX
       PRINT *,N_SHAREDPTS,' SHARED POINTS'
       IF ( N_SHAREDPTS > 0) THEN
          DO I = 1,N_SHAREDPTS
             PRINT *,COORDS(:,SHAREDPTS(I))
          END DO
       END IF
       PRINT *,'INTERIOR?',INTERIOR
       PRINT *,'STAT_NEWPOINT =',STAT_NEWPOINT
       STAT_NEWPOINT = 99
       PRINT *,'NO MORE SUBDIVISION !!!!!'
       CALL WRITE_BERNSTEIN_SERIES1( REGION_C%BEZIER, 'dev_intersection_simple_surface/region_c_bezier.bern' )
       CALL WRITE_BERNSTEIN_SERIES2( REGION_S%BEZIER, 'dev_intersection_simple_surface/region_s_bezier.bern' )
              
       !CALL WRITE_OBB( REGION_C%XYZBOX, 'dev_intersection_simple_surface/xyzbox_c.dat' )
       !CALL WRITE_OBB( REGION_S%XYZBOX, 'dev_intersection_simple_surface/xyzbox_s.dat' )
       RETURN
       !STOP 'NO MORE SUBDIVISION !!!!!'
    end if

    !PRINT *,'NCHILD =',NCHILD
    do jchild = 1,nchild(2)
       if ( nchild(2) == 1 ) then
          newregion_s => region_s
       else
          newregion_s => region_s%child(jchild)
       end if

       do ichild = 1,nchild(1)
          if ( nchild(1) == 1 ) then
             newregion_c => region_c
          else
             newregion_c => region_c%child(ichild)
          end if


          call intersect_curve_surface( &
               root_c, &
               root_s, &
               newregion_c, &
               newregion_s, &
               coords, &
               npts, &
               stat_degeneracy )
       end do
    end do

  end subroutine intersect_curve_surface

























































































  subroutine intersect_curve_surface_elsewhere( &
       b_c, &
       b_s, &
       xyzinter, &
       separable, &
       randomize )
    use mod_math
    use mod_bernstein
    use mod_geometry
    use mod_separation
    implicit none
    type(type_bernstein_series1), intent(in)  :: b_c
    type(type_bernstein_series2), intent(in)  :: b_s
    real(kind=MATHpr),            intent(in)  :: xyzinter(3)
    logical,                      intent(out) :: separable
    logical, optional,            intent(in)  :: randomize
    real(kind=MATHpr)                         :: sep_c(b_c%degr+1,3)
    real(kind=MATHpr)                         :: sep_s((b_s%degr(1)+1)*(b_s%degr(2)+1),3)
    real(kind=MATHpr)                         :: vec(3), rot(3,3)
    integer                                   :: nbcps, nc, ns

    call rearrange_for_separability_test( &
         b_c%coef(1:b_c%degr+1,:), &
         b_c%degr+1, &
         xyzinter, &
         sep_c, &
         nc )

    nbcps = ( b_s%degr(1) + 1 ) * ( b_s%degr(2) + 1 )
    call rearrange_for_separability_test( &
         reshape( b_s%coef(1:b_s%degr(1)+1,1:b_s%degr(2)+1,:), [nbcps,3] ), &
         nbcps, &
         xyzinter, &
         sep_s, &
         ns )

    if ( nc < 1 .or. ns < 1 ) then
       separable = .true.
    else
       if ( present(randomize) ) then
          if ( randomize ) then
             call random_rotation_matrix3d( rot )
             IF ( .false. ) THEN
                PRINT *,'';PRINT *,'';PRINT *,''
                PRINT *,'RANDOM ROTATION MATRIX ='
                CALL PRINT_MAT( ROT )
             END IF
             sep_c(1:nc,:) = matmul( sep_c(1:nc,:), rot )
             sep_s(1:ns,:) = matmul( sep_s(1:ns,:), rot )
          end if
       end if

       call separating_plane( &
            sep_c(1:nc,1:3), &
            sep_s(1:ns,1:3), &
            nc, &
            ns, &
            vec, &
            separable )

       IF ( .FALSE. ) THEN!.NOT.SEPARABLE ) THEN
          PRINT *,'';PRINT *,'';PRINT *,''
          PRINT *,'XYZ_SEP (C) ='
          CALL PRINT_MAT( SEP_C(1:NC,:) )
          PRINT *,'XYZ_SEP (S) ='
          CALL PRINT_MAT( SEP_S(1:NS,:) )     
       END IF
    end if

  end subroutine intersect_curve_surface_elsewhere








  subroutine rearrange_for_separability_test( &
       bcp, &
       nbcp, &
       origin, &
       sep, &
       nsep )
    use mod_math
    use mod_tolerances
    implicit none
    integer,           intent(in)  :: nbcp
    real(kind=MATHpr), intent(in)  :: bcp(nbcp,3)
    real(kind=MATHpr), intent(in)  :: origin(3)
    real(kind=MATHpr), intent(out) :: sep(nbcp,3)
    integer,           intent(out) :: nsep
    real(kind=MATHpr)              :: xyzi(3)
    integer                        :: i

    nsep = 0
    do i = 1,nbcp
       xyzi = bcp(i,:) - origin
       if ( sum(xyzi**2) > EPSxyzsqr ) then
          nsep = nsep + 1
          sep(nsep,:) = xyzi
       end if
    end do

  end subroutine rearrange_for_separability_test








  recursive subroutine add_to_parent_curve( &
       region, &
       id )
    use mod_util
    implicit none
    type(type_curve_region), intent(inout) :: region
    integer,                 intent(in)    :: id

    !IF ( REGION%NPTS > 0 ) PRINT *,'==>>>>> APPEND',ID,'TO',REGION%IPTS
    call append( &
         region%ipts, &
         id )
    region%npts = region%npts + 1

    if ( associated(region%parent) ) then
       call add_to_parent_curve( &
            region%parent, &
            id )
    end if

  end subroutine add_to_parent_curve








  recursive subroutine add_to_parent_surface( &
       region, &
       id )
    use mod_util
    implicit none
    type(type_surface_region), intent(inout) :: region
    integer,                   intent(in)    :: id

    !IF ( REGION%NPTS > 0 ) PRINT *,'==>>>>> APPEND',ID,'TO',REGION%IPTS
    call append( &
         region%ipts, &
         id )
    region%npts = region%npts + 1

    if ( associated(region%parent) ) then
       call add_to_parent_surface( &
            region%parent, &
            id )
    end if

  end subroutine add_to_parent_surface







  subroutine append_tuvxyz_point( &
       coords, &
       npts, &
       tuvxyz )
    use mod_math
    implicit none
    real(kind=MATHpr),              intent(in)    :: tuvxyz(6)
    real(kind=MATHpr), allocatable, intent(inout) :: coords(:,:)
    integer,                        intent(inout) :: npts
    real(kind=MATHpr), allocatable                :: tmp(:,:)

    npts = npts + 1
    if ( .not.allocated(coords) ) allocate( coords(6,1) )
    if ( npts > size(coords,2) ) then
       call move_alloc( from=coords, to=tmp )
       allocate( coords(6,npts) )
       coords(:,1:npts) = tmp(:,1:npts)
       deallocate(tmp)
    end if
    !PRINT *,'APPENDING TUVXYZ=',REAL(TUVXYZ)
    coords(1:6,npts) = tuvxyz(1:6)

  end subroutine append_tuvxyz_point





  subroutine append_vector( &
       vec, &
       dim, &
       array, &
       n )
    use mod_math
    implicit none
    integer,                        intent(in)    :: dim
    real(kind=MATHpr),              intent(in)    :: vec(dim)
    real(kind=MATHpr), allocatable, intent(inout) :: array(:,:)
    integer,                        intent(inout) :: n
    real(kind=MATHpr), allocatable                :: tmp(:,:)

    if ( .not.allocated(array) ) allocate( array(dim,1) )
    if ( dim > size(array,1) ) STOP 'append_vector : dim > size(array,1)'
    if ( n + 1 > size(array,2) ) then
       call move_alloc( from=array, to=tmp )
       allocate( array(dim,n+1) )
       array(:,1:n) = tmp(:,1:n)
       deallocate(tmp)
    end if
    !PRINT *,'APPENDING VEC=',REAL(VEC)
    n = n + 1
    array(1:dim,n) = vec(1:dim)

  end subroutine append_vector






  subroutine subdiv_curve_region( &
       region, &
       ts, &
       stat )
    use mod_math
    ! Subdivides a curve_region region at a given parametric point 'ts'.
    ! The region is subidivided into 2 child regions unless the subdivision point 
    ! is close to one of the region's endpoints.
    ! Only the 'tbox' attribute of the children are assigned 
    ! ( the remaining attributes require additional data that are outside the scope 
    ! of this subroutine ).
    implicit none
    real(kind=MATHpr),       intent(in)            :: ts
    type(type_curve_region), intent(inout), target :: region
    integer,                 intent(out)           :: stat
    real(kind=MATHpr)                              :: t(3)
    integer                                        :: ichild

    if ( associated( region%child ) ) then
       ! the region already has children
       stat = -1
       return
    end if

    t([1,3]) = region%tbox
    t(2) = ts

    ! check if the subivision point is close to one of the region's endpoints
    if ( abs(t(2) - t(1)) < EPSregion .or.  abs(t(3) - t(2)) < EPSregion ) then
       stat = 1

    else
       ! subdivide the region into 2 child regions
       stat = 0
       allocate( region%child(2) )
       do ichild = 1,2
          call init_curve_region( &
               region%child(ichild), &
               t(ichild + [0,1]) )
          region%child(ichild)%parent => region
       end do



    end if

  end subroutine subdiv_curve_region







  subroutine subdiv_surface_region( &
       region, &
       uvs, &
       stat )
    use mod_math
    ! Subdivides a surface_region region at a given parametric point 'uvs'.
    ! If the subdivision point is close to one of the region's boundaries, 
    ! the region is subdivided into 2 child regions, else, it is subdivided 
    ! into 4 child regions.
    ! Only the 'uvbox' and 'parent' attributes of the children are assigned 
    ! ( the remaining attributes require additional data that are outside the scope 
    ! of this subroutine ).
    implicit none
    real(kind=MATHpr),         intent(in)            :: uvs(2)
    type(type_surface_region), intent(inout), target :: region
    integer,                   intent(out)           :: stat
    real(kind=MATHpr)                                :: uv(3,2)
    logical                                          :: degenerate(2)
    integer                                          :: idim, ichild, jchild

    if ( associated( region%child ) ) then
       ! the region already has children
       stat = -1
       return
    end if

    uv([1,3],:) = region%uvbox
    uv(2,:) = uvs

    ! check if the subivision point is close to the region's boundary
    do idim = 1,2
       degenerate(idim) = ( abs(uv(2,idim) - uv(1,idim)) < EPSregion ) .or. &
            ( abs(uv(3,idim) - uv(2,idim)) < EPSregion )
    end do

    if ( all(degenerate) ) then
       ! degenerate subdivision (corner)
       stat = 3
       return
       !STOP 'subdiv_surface_region : degenerate subdivision (corner)'

    elseif ( degenerate(1) ) then
       ! subdivision into 2 child regions along horizontal line v = vs
       stat = 1
       allocate( region%child(2) )
       do jchild = 1,2
          call init_surface_region( &
               region%child(jchild), &
               reshape( [uv([1,3],1), uv(jchild+[0,1],2)], [2,2] ) )
       end do

    elseif ( degenerate(2) ) then
       ! subdivision into 2 child regions along vertical line u = us
       stat = 2
       allocate( region%child(2) )
       do ichild = 1,2
          call init_surface_region( &
               region%child(ichild), &
               reshape( [uv(ichild+[0,1],1), uv([1,3],2)], [2,2] ) )
       end do

    else
       ! subdivision into 4 child regions at point (u,v) = (us,vs)
       stat = 0
       allocate( region%child(4) )
       do jchild = 1,2
          do ichild = 1,2
             call init_surface_region( &
                  region%child( 2*(jchild-1) + ichild ), &
                  reshape( [uv(ichild+[0,1],1), uv(jchild+[0,1],2)], [2,2] ) )
          end do
       end do

    end if

    ! make all children point to their parent, i.e. the current region
    do ichild = 1,size(region%child)
       region%child(ichild)%parent => region
    end do


  end subroutine subdiv_surface_region







  subroutine init_curve_region( &
       region, &
       tbox )
    use mod_math
    ! Initializes a curve_region region with no child, given an interval tbox
    implicit none
    real(kind=MATHpr),       intent(in)  :: tbox(2)
    type(type_curve_region), intent(out) :: region

    region%tbox = tbox

    if ( associated(region%xyzbox) ) deallocate( region%xyzbox )
    !nullify( region%xyzbox )

    if ( associated(region%parent) ) deallocate( region%parent )
    !nullify( region%parent )

    if ( associated(region%child) ) deallocate( region%child )
    !nullify( region%child )

    if ( associated(region%bezier) ) then 
       if ( allocated(region%bezier%coef) ) deallocate( region%bezier%coef )
       deallocate( region%bezier )
    end if
    !nullify( region%bezier )

    if ( allocated(region%ipts) ) deallocate( region%ipts )
    region%npts = 0

  end subroutine init_curve_region







  subroutine init_surface_region( &
       region, &
       uvbox )
    use mod_math
    ! Initializes a surface_region region with no child, given a uvbox
    implicit none
    real(kind=MATHpr),         intent(in)  :: uvbox(2,2)
    type(type_surface_region), intent(out) :: region

    region%uvbox = uvbox

    if ( associated(region%xyzbox) ) deallocate( region%xyzbox )
    !nullify( region%xyzbox )

    if ( associated(region%parent) ) deallocate( region%parent )
    !nullify( region%parent )

    if ( associated(region%child) ) deallocate( region%child )
    !nullify( region%child )

    if ( associated(region%bezier) ) then
       if ( allocated(region%bezier%coef) ) deallocate( region%bezier%coef )
       deallocate( region%bezier )
    end if
    !nullify( region%bezier )

    if ( allocated(region%ipts) ) deallocate( region%ipts )
    region%npts = 0

  end subroutine init_surface_region






  subroutine newton_curve_surface_old( &
       curv, &
       surf, &
       tbox, &
       uvbox, &
       tuv, &
       stat, &
       xyz )
    use mod_math
    use mod_diffgeom
    use mod_tolerances
    ! stat = 0 : converged
    !        1 : not converged
    !        2 : degeneracy
    implicit none
    integer, parameter                           :: nitmax = 1 + max( 5, ceiling(-log10(EPSxyz)) )  !10
    type(type_parametric_curve),   intent(in)    :: curv
    type(type_parametric_surface), intent(in)    :: surf
    real(kind=MATHpr),             intent(in)    :: tbox(2)
    real(kind=MATHpr),             intent(in)    :: uvbox(2,2)
    real(kind=MATHpr),             intent(inout) :: tuv(3)
    integer,                       intent(out)   :: stat
    real(kind=MATHpr),             intent(out)   :: xyz(3)
    real(kind=MATHpr), dimension(3)              :: lowerb, upperb, rng
    real(kind=MATHpr), dimension(3)              :: xyz_c, xyz_s, r, dtuv
    real(kind=MATHpr)                            :: jac(3,3), lambda, ressqr
    logical                                      :: singular
    integer                                      :: rank
    integer                                      :: it

    lowerb = [ tbox(1), uvbox(1,:) ]
    upperb = [ tbox(2), uvbox(2,:) ]
    rng = upperb - lowerb
    lowerb = lowerb - EPsuv*rng
    upperb = upperb + EPsuv*rng
    !PRINT *,'NEWTON :'
    !PRINT *,'TBOX =',REAL(TBOX)
    !PRINT *,'UVBOX =',REAL(UVBOX)
    !PRINT *,'LOWERB =',REAL(LOWERB)
    !PRINT *,'UPPERB =',REAL(UPPERB)

    stat = 1

    do it = 1,nitmax
       call eval( xyz_c, curv, tuv(1) )
       call eval( xyz_s, surf, tuv(2:3) )

       r = xyz_s - xyz_c ! residual vector
       ressqr = sum( r**2 ) ! squared residual 
       !PRINT *,'NEWTON, IT#',IT,', TUV =',REAL(TUV),', RES =',REAL(NORM2(R))
       !PRINT *,'NEWTON, IT#',IT,', RES =',REAL(NORM2(R))
       if ( stat == 0 ) then
          PRINT *,'NEWTON, IT#',IT,', RES =',REAL(NORM2(R))
       end if

       ! convergence criterion
       if ( ressqr < EPSxyzsqr ) then
          !if ( ressqr < MATHeps**2 ) stat = 0

          if ( stat == 0 ) then
             PRINT *,'NEWTON CONVERGED, RES=',SQRT(RESSQR)
             PRINT *,'TUVXYZ =',TUV,XYZ
             xyz = 0.5_MATHpr * ( xyz_s + xyz_c )
             return
          end if

          stat = 0
          !PRINT *,'NEWTON => XYZ =',REAL(XYZ)
          PRINT *,'+++++++ NEWTON CONVERGE #IT',IT,' RES=',SQRT(RESSQR)
          xyz = 0.5_MATHpr * ( xyz_s + xyz_c )
          return
       end if

       ! Jacobian matrix
       call evald( jac(:,1), curv, tuv(1) )
       jac(:,1) = -jac(:,1)
       call evald( jac(:,2), surf, tuv(2:3), 1 )
       call evald( jac(:,3), surf, tuv(2:3), 2 )

       ! solve for Newton step
       call linsolve_QR( &
            dtuv, &
            jac, &
            -r, &
            3, &
            3, &
            singular, &
            rank )

       !PRINT *,'JAC='
       !CALL PRINT_MAT( JAC )

       !if ( singular ) then ! degeneracy
       if ( rank < 3 ) then ! degeneracy
          PRINT *,''
          PRINT *,'IT #',IT
          PRINT *,'TUV =',TUV
          PRINT *,'RES =',SQRT(RESSQR)
          PRINT *,'JACOBIAN ='
          CALL PRINT_MAT( JAC )
          PRINT *,'RANK =',RANK
          PRINT *,'DTUV=',DTUV
          if ( stat == 0 ) then
             stat = 2
             return
          end if
       end if

       ! scale down Newton step to keep the solution inside feasible region
       call nd_box_constraint( &
            tuv, &
            lowerb, &
            upperb, &
            dtuv, &
            lambda )

       if ( lambda < -MATHeps ) then
          ! negative damped factor
          return
       end if

       dtuv = lambda * dtuv
       if ( abs(dtuv(1)) < EPSuv .or. sum(dtuv(2:3)**2) < EPSuvsqr ) then
          ! damped Newton step is too small
          PRINT *,'|DT| =',ABS(DTUV(1)),' |DUV| =',NORM2( DTUV(2:3) )
          return
       end if

       ! update solution
       tuv = tuv + lambda * dtuv

    end do


  end subroutine newton_curve_surface_old











  recursive subroutine free_curve_region_tree( &
       region )
    ! Recursively frees all the curve_regions in the subtree
    ! rooted at the node 'region'
    implicit none
    type(type_curve_region), intent(inout) :: region
    integer                                :: ichild

    if ( associated(region%child) ) then
       do ichild = 1,size(region%child)
          call free_curve_region_tree( region%child(ichild) )
       end do
       deallocate( region%child )
    end if

    if ( associated(region%xyzbox) ) deallocate( region%xyzbox )
    if ( allocated(region%ipts) ) deallocate( region%ipts )
    if ( associated(region%bezier) ) then
       if ( allocated(region%bezier%coef) ) deallocate( region%bezier%coef )
       deallocate( region%bezier )
    end if

  end subroutine free_curve_region_tree






  recursive subroutine free_surface_region_tree( &
       region, &
       no_dealloc_bezier )
    ! Recursively frees all the surface_regions in the subtree
    ! rooted at the node 'region'
    implicit none
    type(type_surface_region), intent(inout) :: region
    logical, optional,         intent(in)    :: no_dealloc_bezier
    integer                                  :: ichild

    if ( associated(region%child) ) then
       do ichild = 1,size(region%child)
          call free_surface_region_tree( region%child(ichild) )
       end do
       deallocate( region%child )
    end if

    if ( associated(region%xyzbox) )     deallocate( region%xyzbox )
    if ( allocated(region%ipts) )        deallocate( region%ipts )

    if ( present(no_dealloc_bezier) ) then
       if ( no_dealloc_bezier ) then
          if ( .not.associated(region%parent) ) return
       end if
    end if

    if ( associated(region%bezier) ) then
       if ( allocated(region%bezier%coef) ) deallocate( region%bezier%coef )
       deallocate( region%bezier )
    end if


  end subroutine free_surface_region_tree








  subroutine export_curve_region_tree( &
       root, &
       filename )
    use mod_util
    ! Exports a curve_region subtree to a given file.
    ! (essentially for debugging purposes)
    ! See the 'write_curve_region' subroutine for details.
    implicit none
    type(type_curve_region), intent(in) :: root
    character(*),            intent(in) :: filename
    integer                             :: file_unit

    call get_free_unit( file_unit )

    open( &
         unit = file_unit, &
         file = filename, &
         action = 'write' )

    call write_curve_region( &
         root, &
         file_unit )

    close( file_unit )

    PRINT *,'curve_region tree written in ',filename

  end subroutine export_curve_region_tree






  subroutine export_surface_region_tree( &
       root, &
       filename )
    use mod_util
    ! Exports a surface_region subtree to a given file.
    ! (essentially for debugging purposes)
    ! See the 'write_surface_region' subroutine for details.
    implicit none
    type(type_surface_region), intent(in) :: root
    character(*),              intent(in) :: filename
    integer                               :: file_unit

    call get_free_unit( file_unit )

    open( &
         unit = file_unit, &
         file = filename, &
         action = 'write' )

    call write_surface_region( &
         root, &
         file_unit )

    close( file_unit )

    PRINT *,'surface_region tree written in ',filename

  end subroutine export_surface_region_tree






  recursive subroutine write_curve_region( &
       region, &
       file_unit )
    ! Writes the data contained in a (leaf) curve_region.
    implicit none
    type(type_curve_region), intent(in) :: region
    integer,                 intent(in) :: file_unit
    integer                             :: ichild

    if ( associated(region%child) ) then
       ! the current region has child regions, carry on the recursion
       do ichild = 1,size(region%child)
          call write_curve_region( &
               region%child(ichild), &
               file_unit )
       end do
    else
       ! the current region is a leaf, write its data
       write (file_unit,*) region%tbox
    end if

  end subroutine write_curve_region





  recursive subroutine write_surface_region( &
       region, &
       file_unit )
    ! Writes the data contained in a (leaf) surface_region.
    implicit none
    type(type_surface_region), intent(in) :: region
    integer,                   intent(in) :: file_unit
    integer                               :: ichild

    if ( associated(region%child) ) then
       ! the current region has child regions, carry on the recursion
       do ichild = 1,size(region%child)
          call write_surface_region( &
               region%child(ichild), &
               file_unit )
       end do
    else
       ! the current region is a leaf, write its data
       write (file_unit,*) region%uvbox(:,1), region%uvbox(:,2)
    end if

  end subroutine write_surface_region










  
  
  
  subroutine newton_curve_surface( &
       curv, &
       surf, &
       tbox, &
       uvbox, &
       tuv, &
       stat, &
       xyz )
    use mod_math
    use mod_diffgeom
    use mod_tolerances    
    ! stat = 0 : converged
    !        1 : not converged
    !        2 : degeneracy
    implicit none
    real(kind=MATHpr), parameter                 :: THRESHOLD = real(1.e-2, kind=MATHpr)
    real(kind=MATHpr), parameter                 :: EPS = real(1.e-12, kind=MATHpr)
    real(kind=MATHpr), parameter                 :: TOL = EPSxyz
    real(kind=MATHpr), parameter                 :: EPSsqr = EPS**2
    real(kind=MATHpr), parameter                 :: TOLsqr = TOL**2
    integer, parameter                           :: nitmax = ceiling(-log10(EPS))
    integer, parameter                           :: nitcheck = 5

    type(type_parametric_curve),   intent(in)    :: curv
    type(type_parametric_surface), intent(in)    :: surf
    real(kind=MATHpr),             intent(in)    :: tbox(2)
    real(kind=MATHpr),             intent(in)    :: uvbox(2,2)
    real(kind=MATHpr),             intent(inout) :: tuv(3)
    integer,                       intent(out)   :: stat
    real(kind=MATHpr),             intent(out)   :: xyz(3)
    real(kind=MATHpr), dimension(3)              :: lowerb, upperb, rng
    real(kind=MATHpr), dimension(3)              :: xyz_c, xyz_s, r
    real(kind=MATHpr), dimension(3)              :: tuvtmp, dtuv
    real(kind=MATHpr)                            :: rescheck, res, restmp
    real(kind=MATHpr)                            :: jac(3,3), lambda
    logical                                      :: singular
    integer                                      :: rank
    integer                                      :: it

    ! Feasible t,u,v-domain
    lowerb = [ tbox(1), uvbox(1,:) ]
    upperb = [ tbox(2), uvbox(2,:) ]
    rng = upperb - lowerb
    lowerb = lowerb - EPsuv*rng
    upperb = upperb + EPsuv*rng


    stat = 1
    restmp = huge(1._MATHpr)
    newton_iteration : do it = 1,nitmax

       ! position vector
       call eval( xyz_c, curv, tuv(1) )   ! curve
       call eval( xyz_s, surf, tuv(2:3) ) ! surface

       r = xyz_s - xyz_c ! residual vector
       res = sum( r**2 ) ! squared norm of residual vector
       !PRINT *,'NEWTON, IT#',IT,', RES =',REAL(NORM2(R))
       !PRINT *,NORM2(R)
       
       ! check signs of convergence
       if ( it == 1 ) rescheck = THRESHOLD * res
       if ( it > nitcheck .and. res > rescheck ) then
       !   PRINT *,'NO SIGN OF CONVERGENCE, STOP NEWTON ITERATION'
          return
       end if

       ! convergence criterion
       if ( res < TOLsqr ) then
          stat = 0
          if ( res < restmp ) then
             restmp = res
             tuvtmp = tuv
             xyz = 0.5_MATHpr * ( xyz_s + xyz_c )
          end if
          if ( restmp < EPSsqr ) then
             tuv = tuvtmp
             return
          end if
       end if

       ! Jacobian matrix
       call evald( jac(:,1), curv, tuv(1) )
       jac(:,1) = -jac(:,1)
       call evald( jac(:,2), surf, tuv(2:3), 1 )
       call evald( jac(:,3), surf, tuv(2:3), 2 )

       ! solve for Newton step
       call linsolve_QR( &
            dtuv, &
            jac, &
            -r, &
            3, &
            3, &
            singular, &
            rank )

       if ( rank < 3 ) then ! degeneracy
          PRINT *,''
          PRINT *,'IT #',IT
          PRINT *,'TUV =',TUV
          PRINT *,'RES =',SQRT(RES)
          PRINT *,'JACOBIAN ='
          CALL PRINT_MAT( JAC )
          PRINT *,'RANK =',RANK
          PRINT *,'DTUV=',DTUV
          if ( stat == 0 ) then
             stat = 2
             return
          end if
       end if

       ! scale down Newton step to keep the solution inside feasible region
       call nd_box_constraint( &
            tuv, &
            lowerb, &
            upperb, &
            dtuv, &
            lambda )

       if ( lambda < -MATHeps ) return ! negative damped factor

       dtuv = lambda * dtuv
       if ( abs(dtuv(1)) < EPSuv .or. sum(dtuv(2:3)**2) < EPSuvsqr ) then
          ! damped Newton step is too small
          !PRINT *,'|DT| =',ABS(DTUV(1)),' |DUV| =',NORM2( DTUV(2:3) )
          return
       end if

       ! update solution
       tuv = tuv + lambda * dtuv

    end do newton_iteration

  end subroutine newton_curve_surface






end program dev_intersection_simple_surface
