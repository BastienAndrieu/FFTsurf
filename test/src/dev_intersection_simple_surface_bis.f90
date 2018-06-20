program dev_intersection_simple_surface

  use mod_util
  use mod_math
  use mod_regiontree
  use mod_polynomial
  use mod_diffgeom2
  use mod_separation  

  implicit none

  type type_intersection_curve
     real(kind=fp) :: uvxyz(7,2)
  end type type_intersection_curve

  type type_listcurves
     integer                       :: nc = 0
     type(type_intersection_curve) :: curve(100)
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
  type(type_listcurves)      :: interdat

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
          '/stck/bandrieu/Bureau/coeffstest/C' // strnum // '_test' // strnum2 // '.txt', &
                                !'/home/bastien/Bureau/coeffstest/C' // strnum // '_test' // strnum2 // '.txt', &
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

     allocate( region(isurf)%ptr%poly )
     call cheb2bern_poly( &
          surf(isurf)%x, &
          region(isurf)%ptr%poly )
  end do


  
  nuvxyz = 0
  stat_degeneracy = 0
  call system_clock( tic, count_rate )
  call intersect_simple_surfaces( &
       surfroot, &
       region, &
       [1._fp, 0._fp, 0._fp], &
       interdat, &
       uvxyz, &
       nuvxyz, &
       stat_degeneracy )
  call system_clock( toc )
  PRINT *,''; PRINT *,''; PRINT *,''
  PRINT *,'ELAPSED =',REAL( TOC - TIC ) / REAL( COUNT_RATE )

  PRINT *,'STAT_DEGENERACY =',STAT_DEGENERACY

  call get_free_unit( fileunit )
  open( &
       unit=fileunit, &
       file='dev_intersection_simple_surface/curves.dat', &
       action='write' )
  write ( fileunit, * ) interdat%nc
  do i = 1,interdat%nc
     write ( fileunit, * ) interdat%curve(i)%uvxyz(:,1)
     write ( fileunit, * ) interdat%curve(i)%uvxyz(:,2)
  end do
  close( fileunit )


  do isurf = 1,2
     write (strnum,'(I1)') isurf
     call export_region_tree( &
          region(isurf)%ptr, &
          'dev_intersection_simple_surface/tree_' // strnum // '.dat' )

     call free_region_tree( region(isurf)%ptr )
     call free_polynomial( region(isurf)%ptr%poly )
     deallocate( region(isurf)%ptr%poly )

     call free_polynomial( surf(isurf)%x )
     call free_polynomial( surf(isurf)%xu )
     call free_polynomial( surf(isurf)%xv )
     call free_polynomial( surf(isurf)%xuu )
     call free_polynomial( surf(isurf)%xuv )
     call free_polynomial( surf(isurf)%xvv )

     !deallocate( root(isurf)%uvbox )
  end do


  if ( allocated( ch2be_matrices ) ) then
     do i = 1,size(ch2be_matrices)
        if ( allocated(ch2be_matrices(i)%mat) ) deallocate( ch2be_matrices(i)%mat )
     end do
     deallocate( ch2be_matrices )
  end if

  if ( allocated(uvxyz) ) deallocate( uvxyz )



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










  subroutine intersect_border_surface( &
       root_s, &
       region, &
       icurv, &
       ivar, &
       ival, &
       uvxyz, &
       nuvxyz, &
       stat_degeneracy )
    use mod_util
    use mod_math
    use mod_chebyshev2
    use mod_polynomial
    use mod_diffgeom2
    use mod_regiontree
    use mod_tolerances
    implicit none
    LOGICAL, PARAMETER :: DEBUG = .false.
    integer, parameter                        :: npts_init = 10
    type(ptr_surface),          intent(in)    :: root_s(2)
    type(ptr_region),           intent(inout) :: region(2)
    integer,                    intent(in)    :: icurv, ivar, ival
    real(kind=fp), allocatable, intent(inout) :: uvxyz(:,:)
    integer,                    intent(inout) :: nuvxyz
    integer,                    intent(inout) :: stat_degeneracy
    type(type_polynomial)                     :: regc
    type(type_curve)                          :: root_c
    type(type_region)                         :: region_c
    type(type_region)                         :: region_s
    real(kind=fp), allocatable                :: tuvxyz(:,:)
    integer                                   :: ntuvxyz
    real(kind=fp)                             :: uv(2,2)
    integer                                   :: isurf, ipt, jpt

    !IF ( DEBUG ) THEN
       !PRINT *,''; PRINT *,''; PRINT *,''; PRINT *,''
       !PRINT *,'   ICURV, IVAR, IVAL =', ICURV, IVAR, IVAL
    !END IF

    isurf = 1 + mod(icurv,2)

    ! convert surface border to parametric curve
    ! <-------------------------------------------+
    call chgvar2( &                               !
        root_s(icurv)%ptr%x, &                    !
        regc, &                                   !-- METTRE DANS INTERSECT_SIMPLE_SURFACES
        region(icurv)%ptr%uvbox([1,3]), &         !
        region(icurv)%ptr%uvbox([2,4]) )          !
    ! <-------------------------------------------+
    call bivar2univar( &
         regc, &!root_s(icurv)%ptr%x, &
         root_c%x, &
         ivar, &
         ival )       

    call economize1( root_c%x, EPSmath )
    call compute_deriv1( root_c )
    call compute_deriv2( root_c )

    IF ( DEBUG ) THEN
       call write_polynomial( root_c%x, 'dev_intersection_simple_surface/root_c_x.cheb' )
       call write_polynomial( root_c%xt, 'dev_intersection_simple_surface/root_c_xt.cheb' )
       call write_polynomial( root_s(isurf)%ptr%x, 'dev_intersection_simple_surface/root_s_x.cheb' )
       call write_polynomial( root_s(isurf)%ptr%xu, 'dev_intersection_simple_surface/root_s_xu.cheb' )
       call write_polynomial( root_s(isurf)%ptr%xv, 'dev_intersection_simple_surface/root_s_xv.cheb' )
    END IF

    IF ( .false. ) then !MAXVAL( REGION(ISURF)%PTR%UVBOX([2,4]) - REGION(ISURF)%PTR%UVBOX([1,3]) ) < 1._FP ) THEN
       CALL WRITE_POLYNOMIAL( ROOT_S(ISURF)%PTR%X, 'dev_intersection_simple_surface/root_s.cheb' )

       CALL CHGVAR2( &
            ROOT_S(ISURF)%PTR%X, &
            REGC, &
            REGION(ISURF)%PTR%UVBOX([1,3]), &
            REGION(ISURF)%PTR%UVBOX([2,4]) )
       PRINT *,'UVBOX =',REGION(ISURF)%PTR%UVBOX
       CALL WRITE_POLYNOMIAL( REGC, 'dev_intersection_simple_surface/region_s.cheb' )
       CALL WRITE_POLYNOMIAL( REGION(ISURF)%PTR%POLY, 'dev_intersection_simple_surface/region_s.bern' )

       CALL CHGVAR2( &
            ROOT_S(ISURF)%PTR%X, &
            REGC, &
            REGION(ISURF)%PTR%PARENT%UVBOX([1,3]), &
            REGION(ISURF)%PTR%PARENT%UVBOX([2,4]) )
       PRINT *,'PARENT''S UVBOX =',REGION(ISURF)%PTR%PARENT%UVBOX
       CALL WRITE_POLYNOMIAL( REGC, 'dev_intersection_simple_surface/parent_s.cheb' )
       CALL WRITE_POLYNOMIAL( REGION(ISURF)%PTR%PARENT%POLY, 'dev_intersection_simple_surface/parent_s.bern' )
       STOP 'DEBUG REGIONS'
    END IF

    ! initialize curve region tree
    call init_region( &
         region_c, &
         1, &
         [ -1._fp, 1._fp ] ) 

    ! initialize surface region tree
    call init_region( &
         region_s, &
         2, &
         region(isurf)%ptr%uvbox ) 

    ! compute curve Bezier control points
    allocate( region_c%poly )
    call cheb2bern_poly( &
         root_c%x, &
         region_c%poly )

    ! copy surface Bezier control points
    region_s%poly => region(isurf)%ptr%poly

    IF ( DEBUG ) THEN
       CALL WRITE_POLYNOMIAL( REGION_C%POLY, 'dev_intersection_simple_surface/root_c_bezier.bern' )
       CALL WRITE_POLYNOMIAL( REGION_S%POLY, 'dev_intersection_simple_surface/root_s_bezier.bern' )
    END IF


    



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

    IF ( DEBUG ) THEN
       open( unit=13, file='dev_intersection_simple_surface/tuv_xyz.dat', action='write' )
       if ( ntuvxyz < 1 ) then
          write (13,*) ''
       else
          do ipt = 1,ntuvxyz
             write (13,*) tuvxyz(:,ipt)
          end do
       end if
       close(13)

       call export_region_tree( &
            region_s, &
            'dev_intersection_simple_surface/tree_s.dat' )
       call export_region_tree( &
            region_c, &
            'dev_intersection_simple_surface/tree_c.dat' )
    END IF

    IF ( STAT_DEGENERACY > 10 ) THEN
       PRINT *,'STAT_DEGENERACY =',STAT_DEGENERACY
       RETURN
       !STOP '********************'
    END IF

    !PRINT  *,'';PRINT  *,'';PRINT  *,'';PRINT  *,'';PRINT  *,'';PRINT  *,'';
    !PRINT *,'     BEFORE :',NUVXYZ,' POINTS'
    !IF ( NUVXYZ > 0 ) CALL PRINT_MAT( TRANSPOSE(UVXYZ(:,1:NUVXYZ)) )
    outer : do ipt = 1,ntuvxyz
       do jpt = 1,nuvxyz
          if ( sum( ( tuvxyz(4:6,ipt) - uvxyz(5:7,jpt) )**2 ) < EPSxyzsqr ) cycle outer
       end do

       uv(ivar,icurv) = region(icurv)%ptr%uvbox(2*(ivar-1)+ival) !region(icurv)%ptr%uvbox(ival,ivar)
       uv(1+mod(ivar,2),icurv) = n1p12ab( &
            tuvxyz(1,ipt), &
            region(icurv)%ptr%uvbox(2*(1+mod(ivar,2))-1), &!region(icurv)%ptr%uvbox(1,1+mod(ivar,2)), &
            region(icurv)%ptr%uvbox(2*(1+mod(ivar,2))) )!region(icurv)%ptr%uvbox(2,1+mod(ivar,2)) )
       uv(:,isurf) = tuvxyz(2:3,ipt)
       !PRINT *,uv,tuvxyz(4:6,ipt)
       !PRINT *,'       XYZ  =',TUVXYZ(4:6,IPT)
       !PRINT *,'       T    =',TUVXYZ(1,IPT)
       !PRINT *,'       UVBOX=',REGION(ICURV)%PTR%UVBOX(2*IVAR+[-1,0])
       !PRINT *,'       UV   =',UV(:,ICURV)

       !CALL WRITE_POLYNOMIAL( REGION_C%POLY, 'dev_intersection_simple_surface/root_c_bezier.bern' )

       call append_vector( &
            [ uv(:,1), uv(:,2), tuvxyz(4:6,ipt) ], &
            7, &
            uvxyz, &
            nuvxyz )

    end do outer
    !PRINT *,'     AFTER  :',NUVXYZ,' POINTS'
    !IF ( NUVXYZ > 0 ) CALL PRINT_MAT( TRANSPOSE(UVXYZ(:,1:NUVXYZ)) )

    call free_polynomial( region_c%poly )
    deallocate( region_c%poly )

    nullify( region_s%poly )

    call free_region_tree( region_s )
    call free_region_tree( region_c )

    call free_polynomial( root_c%x )
    call free_polynomial( root_c%xt )
    call free_polynomial( root_c%xtt )

  end subroutine intersect_border_surface










































  subroutine append_vector( &
       vec, &
       dim, &
       array, &
       n )
    use mod_math
    implicit none
    integer,                    intent(in)    :: dim
    real(kind=fp),              intent(in)    :: vec(dim)
    real(kind=fp), allocatable, intent(inout) :: array(:,:)
    integer,                    intent(inout) :: n
    real(kind=fp), allocatable                :: tmp(:,:)

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










  subroutine intersect_curve_surface_elsewhere( &
       b_c, &
       b_s, &
       xyzinter, &
       separable, &
       randomize )
    use mod_math
    use mod_polynomial
    use mod_geometry
    use mod_separation
    implicit none
    type(type_polynomial), intent(in)  :: b_c
    type(type_polynomial), intent(in)  :: b_s
    real(kind=fp),         intent(in)  :: xyzinter(3)
    logical,               intent(out) :: separable
    logical, optional,     intent(in)  :: randomize
    real(kind=fp)                      :: sep_c(b_c%degr(1)+1,3)
    real(kind=fp)                      :: sep_s((b_s%degr(1)+1)*(b_s%degr(2)+1),3)
    real(kind=fp)                      :: vec(3), rot(3,3)
    integer                            :: nbcps, nc, ns

    call rearrange_for_separability_test( &
         b_c%coef(1:b_c%degr(1)+1,1:b_c%dim,1), &
         b_c%degr(1)+1, &
         xyzinter, &
         sep_c, &
         nc )

    nbcps = ( b_s%degr(1) + 1 ) * ( b_s%degr(2) + 1 )
    call rearrange_for_separability_test( &
         reshape( b_s%coef(1:b_s%degr(1)+1,1:b_s%degr(2)+1,1:b_s%dim), [nbcps,3] ), &
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
             !IF ( .false. ) THEN
             !   PRINT *,'';PRINT *,'';PRINT *,''
             !   PRINT *,'RANDOM ROTATION MATRIX ='
             !   CALL PRINT_MAT( ROT )
             !END IF
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

       !IF ( .FALSE. ) THEN!.NOT.SEPARABLE ) THEN
       !   PRINT *,'';PRINT *,'';PRINT *,''
       !   PRINT *,'XYZ_SEP (C) ='
       !   CALL PRINT_MAT( SEP_C(1:NC,:) )
       !   PRINT *,'XYZ_SEP (S) ='
       !   CALL PRINT_MAT( SEP_S(1:NS,:) )     
       !END IF
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
    integer,       intent(in)  :: nbcp
    real(kind=fp), intent(in)  :: bcp(nbcp,3)
    real(kind=fp), intent(in)  :: origin(3)
    real(kind=fp), intent(out) :: sep(nbcp,3)
    integer,       intent(out) :: nsep
    real(kind=fp)              :: xyzi(3)
    integer                    :: i

    nsep = 0
    do i = 1,nbcp
       xyzi = bcp(i,:) - origin
       if ( sum(xyzi**2) > EPSxyzsqr ) then
          nsep = nsep + 1
          sep(nsep,:) = xyzi
       end if
    end do

  end subroutine rearrange_for_separability_test












  subroutine newton_curve_surface( &
       curv, &
       surf, &
       tbox, &
       uvbox, &
       tuv, &
       stat, &
       xyz )
    use mod_math
    use mod_diffgeom2
    use mod_tolerances    
    ! stat = 0 : converged
    !        1 : not converged
    !        2 : degeneracy
    implicit none
    real(kind=fp), parameter          :: THRESHOLD = real(1.e-2, kind=fp)
    real(kind=fp), parameter          :: EPS = real(1.e-12, kind=fp)
    real(kind=fp), parameter          :: TOL = EPSxyz
    real(kind=fp), parameter          :: EPSsqr = EPS**2
    real(kind=fp), parameter          :: TOLsqr = TOL**2
    integer, parameter                :: nitmax = ceiling(-log10(EPS))
    integer, parameter                :: nitcheck = 5

    type(type_curve),   intent(in)    :: curv
    type(type_surface), intent(in)    :: surf
    real(kind=fp),      intent(in)    :: tbox(2)
    real(kind=fp),      intent(in)    :: uvbox(4)
    real(kind=fp),      intent(inout) :: tuv(3)
    integer,            intent(out)   :: stat
    real(kind=fp),      intent(out)   :: xyz(3)
    real(kind=fp), dimension(3)       :: lowerb, upperb, rng
    real(kind=fp), dimension(3)       :: xyz_c, xyz_s, r
    real(kind=fp), dimension(3)       :: tuvtmp, dtuv
    real(kind=fp)                     :: rescheck, res, restmp
    real(kind=fp)                     :: jac(3,3), lambda
    logical                           :: singular
    integer                           :: rank
    integer                           :: it

    ! Feasible t,u,v-domain
    lowerb = [ tbox(1), uvbox([1,3]) ]
    upperb = [ tbox(2), uvbox([2,4]) ]
    rng = upperb - lowerb
    lowerb = lowerb - EPsuv*rng
    upperb = upperb + EPsuv*rng
    
    stat = 1
    restmp = huge(1._fp)
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
             xyz = 0.5_fp * ( xyz_s + xyz_c )
          end if
          if ( restmp < EPSsqr ) then
             tuv = tuvtmp
             return
          end if
       end if

       ! Jacobian matrix
       call evald1( jac(:,1), curv, tuv(1) )
       jac(:,1) = -jac(:,1)
       call evald1( jac(:,2), surf, tuv(2:3), 1 )
       call evald1( jac(:,3), surf, tuv(2:3), 2 )

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

       if ( lambda < -EPSmath ) return ! negative damped factor

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
    use mod_diffgeom2
    use mod_obb
    use mod_regiontree
    use mod_tolerances
    implicit none
    type(type_curve),           intent(in)            :: root_c
    type(type_surface),         intent(in)            :: root_s
    type(type_region),          intent(inout), target :: region_c
    type(type_region),          intent(inout), target :: region_s
    real(kind=fp), allocatable, intent(inout)         :: coords(:,:) ! t,u,v,x,y,z
    integer,                    intent(inout)         :: npts
    integer,                    intent(inout)         :: stat_newpoint
    !logical, allocatable                              :: mask(:)
    integer, allocatable                              :: sharedpts(:)
    integer                                           :: n_sharedpts
    logical                                           :: separable, overlap
    real(kind=fp)                                     :: tuv_subdiv(3)
    logical                                           :: interior(2)
    real(kind=fp)                                     :: tuv(3), xyz(3)
    integer                                           :: stat_subdiv, nchild(2)
    type(type_region), pointer                        :: newregion_c
    type(type_region), pointer                        :: newregion_s
    integer                                           :: ipt, jpt, i, j, k, ichild, jchild

    if ( stat_degeneracy > 1 ) return

    !PRINT *,''
    !PRINT *,''
    !PRINT *,' TBOX =',REGION_C%UVBOX
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


    call inherit_points( &
         region_c, &
         coords(1,1:npts), &
         npts )

    call inherit_points( &
         region_s, &
         coords(2:3,1:npts), &
         npts )



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
    !PRINT *,N_SHAREDPTS,' SHARED POINT(S)'
    !IF ( N_SHAREDPTS > 0 ) THEN
    !   DO IPT = 1,N_SHAREDPTS
    !      PRINT *,COORDS(:,SHAREDPTS(IPT))
    !   END DO
    !END IF
    !END IF
    IF (.FALSE.) THEN
       if ( n_sharedpts == 0 ) then
          ! check endpoint/corner pairs for possible intersection point
          do k = 1,region_s%poly%degr(2)+1,region_s%poly%degr(2)
             do j = 1,region_s%poly%degr(1)+1,region_s%poly%degr(1)
                do i = 1,region_c%poly%degr(1)+1,region_c%poly%degr(1)
                   !PRINT *, I, J, K, NORM2( region_c%poly%coef(i,:,1) - region_s%poly%coef(j,k,:) )
                   if ( sum( (region_c%poly%coef(i,:,1) - region_s%poly%coef(j,k,:))**2 ) < EPSxyzsqr ) then
                      xyz = 0.5_fp * ( region_c%poly%coef(i,:,1) + region_s%poly%coef(j,k,:) )
                      tuv(1) = real( i-1, kind=fp ) / real( region_c%poly%degr(1), kind=fp )
                      tuv(2) = real( j-1, kind=fp ) / real( region_s%poly%degr(1), kind=fp )
                      tuv(3) = real( k-1, kind=fp ) / real( region_s%poly%degr(2), kind=fp )
                      tuv = -1._fp + 2._fp * tuv

                      call append_vector( &
                           [tuv,xyz], &
                           6, &
                           coords, &
                           npts ) 

                      call append( &
                           sharedpts, &
                           npts, &
                           newlength=n_sharedpts )

                      call add_point_bottom_up( region_c, npts )
                      call add_point_bottom_up( region_s, npts )
                   end if
                end do
             end do
          end do
       end if
    END IF


    if ( n_sharedpts == 1 ) then
       ipt = sharedpts(1)
       interior(1) = is_in_interval_strict( coords(1,ipt), region_c%uvbox(1), region_c%uvbox(2) )
       interior(2) = ( &
            is_in_interval_strict( coords(2,ipt), region_s%uvbox(1), region_s%uvbox(2) ) .and. &
            is_in_interval_strict( coords(3,ipt), region_s%uvbox(3), region_s%uvbox(4) ) )

       if ( any(interior) ) then
          separable = .false. 
       else
          ! check if the curve and surface regions can intersect at other (not yet discovered) points
          call intersect_curve_surface_elsewhere( &
               region_c%poly, &
               region_s%poly, &
               coords(4:6,sharedpts(1)), &
               separable, &
               randomize=.true. )
       end if
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
          interior(1) = is_in_interval_strict( coords(1,ipt), region_c%uvbox(1), region_c%uvbox(2) )
          interior(2) = ( &
               is_in_interval_strict( coords(2,ipt), region_s%uvbox(1), region_s%uvbox(2) ) .and. &
               is_in_interval_strict( coords(3,ipt), region_s%uvbox(3), region_s%uvbox(4) ) )
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
               region_c%poly%coef(1:region_c%poly%degr(1)+1,1:3,1), &
               region_c%poly%degr(1), &
               region_c%xyzbox )
       end if

       if ( .not.associated(region_s%xyzbox) ) then
          allocate( region_s%xyzbox )
          call bernOBB2( &
               region_s%poly%coef(1:region_s%poly%degr(1)+1,1:region_s%poly%degr(2)+1,1:3), &
               region_s%poly%degr, &
               region_s%xyzbox )
       end if

       call overlap_OBBs( &
            region_c%xyzbox, &
            region_s%xyzbox, &
            overlap )

       if ( .not.overlap ) return


       !IF ( &
       !     REGION_C%UVBOX(2) < REGION_C%UVBOX(1) + 1.D-4 .OR. &
       !     REGION_S%UVBOX(2) < REGION_S%UVBOX(1) + 1.D-4 .OR. &
       !     REGION_S%UVBOX(4) < REGION_S%UVBOX(3) + 1.D-4 ) THEN
       !IF (.FALSE.) THEN
       !   PRINT *,' TBOX =',REGION_C%UVBOX
       !   PRINT *,'UVBOX =',REGION_S%UVBOX
       !   CALL WRITE_OBB( REGION_C%XYZBOX, 'dev_intersection_simple_surface/xyzbox_c.dat' )
       !   CALL WRITE_OBB( REGION_S%XYZBOX, 'dev_intersection_simple_surface/xyzbox_s.dat' )
       !   CALL WRITE_POLYNOMIAL( REGION_C%POLY, 'dev_intersection_simple_surface/region_c_bezier.bern' )
       !   CALL WRITE_POLYNOMIAL( REGION_S%POLY, 'dev_intersection_simple_surface/region_s_bezier.bern' )
       !   PRINT *,'--> VERIF OBBs'
       !   STAT_NEWPOINT = 44
       !   RETURN
       !END IF

    end if

    !PRINT *,'INTERIOR?',INTERIOR
    !PRINT *,'ALL(.NOT.INTERIOR) =',ALL(.NOT.INTERIOR)
    if ( all(.not.interior) ) then ! <------------------------------------------------------------------------------------------+
       !! Search for a new intersection point                                                                                   !
       tuv = 0.5_fp * [ &
            region_c%uvbox(1) + region_c%uvbox(2), &
            region_s%uvbox([1,3]) + region_s%uvbox([2,4]) ]
       !PRINT *,''; PRINT *,'';
       !PRINT *,'NEWTON <--- TUV0 =',TUV
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
       call newton_curve_surface( &
            root_c, &
            root_s, &
            region_c%uvbox, &
            region_s%uvbox, &
            tuv, &
            stat_newpoint, &
            xyz ) 
       !PRINT *,'----------------------'
       !PRINT *,'STAT_NEWPOINT =',STAT_NEWPOINT

       ! if we just found a degenerate point, return and report degeneracy
       if ( stat_newpoint > 1 ) return

       if ( stat_newpoint == 0 ) then ! <--------------------------------------------------------------------------------+
          ! check if the "new" point has not already been discovered                                                     !
          do ipt = 1,npts ! <-------------------------------------------------------------------+                        !
             if ( sum( ( xyz - coords(4:6,ipt) )**2 ) < EPSxyzsqr ) then ! <---------+          !                        ! 
                stat_newpoint = 1                                                    !          !                        !
                exit                                                                 !          !                        !
             end if ! <--------------------------------------------------------------+          !                        !
          end do ! <----------------------------------------------------------------------------+                        !
       end if ! <--------------------------------------------------------------------------------------------------------+

       if ( stat_newpoint == 0 ) then
          ! if we just found a new intersection point, add it to the lists, and subdivide at that point
          !PRINT *,'NEW POINT =',TUV, XYZ
          !PRINT *,'IN   TBOX =',REGION_C%TBOX
          !PRINT *,'IN  UVBOX =',REGION_S%UVBOX
          !PRINT *,'';PRINT *,'';PRINT *,'';
          call append_vector( &
               [tuv,xyz], &
               6, &
               coords, &
               npts )
          !PRINT *,'COORDS ='
          !CALL PRINT_MAT( TRANSPOSE(COORDS(:,1:NPTS)) )
          call add_point_bottom_up( region_c, npts )
          call add_point_bottom_up( region_s, npts )

          tuv_subdiv = tuv
       else
          ! else, subdivide at the parametric midpoint
          tuv_subdiv = 0.5_fp * [ &
               region_c%uvbox(1) + region_c%uvbox(2), &
               region_s%uvbox([1,3]) + region_s%uvbox([2,4]) ]
       end if
    else
       IF ( N_SHAREDPTS < 1 ) STOP 'INTERIOR BUT N_SHAREDPTS < 1'
       tuv_subdiv = coords(1:3,sharedpts(1))
       !PRINT *,'TUV_SUBDIV =',REAL( TUV_SUBDIV )
    end if



    !PRINT *,'TUV_SUBDIV =',TUV_SUBDIV
    !! Subdivide the curve region
    call subdiv_region( &
         region_c, &
         tuv_subdiv(1), &
         stat_subdiv )

    !PRINT *,'STAT_SUBDIV C =',STAT_SUBDIV
    tuv_subdiv(1) = 0.5_fp * ( ab2n1p1( tuv_subdiv(1), region_c%uvbox(1), region_c%uvbox(2) ) + 1._fp )
    !PRINT *,'SIZE(CHILD) =',SIZE(REGION_C%CHILD)

    if ( stat_subdiv == 1 ) then ! the subdivision point is at one of the curve's endpoints
       nchild(1) = 1
    else ! the subdivision point is interior to the curve
       nchild(1) = size(region_c%child)
       IF ( NCHILD(1) /= 2 ) THEN
          PRINT *,'ERROR : NCHILD(1) =',NCHILD(1)
          STOP
       END IF

       if ( stat_subdiv == 0 ) then ! the curve region has no children yet
          allocate( region_c%child(1)%poly, region_c%child(2)%poly )
          call subdiv_bezier1( &
               region_c%poly, &
               tuv_subdiv(1), &
               bl=region_c%child(1)%poly, &
               br=region_c%child(2)%poly )
       end if

    end if


    !! Subdivide the surface region
    call subdiv_region( &
         region_s, &
         tuv_subdiv(2:3), &
         stat_subdiv )


    !PRINT *,'STAT_SUBDIV S =',STAT_SUBDIV
    tuv_subdiv(2) = 0.5_fp * ( ab2n1p1( tuv_subdiv(2), region_s%uvbox(1), region_s%uvbox(2) ) + 1._fp )
    tuv_subdiv(3) = 0.5_fp * ( ab2n1p1( tuv_subdiv(3), region_s%uvbox(3), region_s%uvbox(4) ) + 1._fp )

    if ( stat_subdiv == 2 ) then ! the subdivision point is at one of the curve's corners
       nchild(2) = 1
    else ! the subdivision point is interior to the surface
       nchild(2) = size(region_s%child)
    end if

    if ( stat_subdiv == 0 ) then ! 4 children
       allocate( &
            region_s%child(1)%poly, region_s%child(2)%poly, &
            region_s%child(3)%poly, region_s%child(4)%poly )
       call subdiv_bezier2( &
            region_s%poly, &
            tuv_subdiv(2:3), &
            bsw=region_s%child(1)%poly, &
            bse=region_s%child(2)%poly, &
            bnw=region_s%child(3)%poly, &
            bne=region_s%child(4)%poly )

    elseif ( stat_subdiv == 1 ) then ! 2 children
       allocate( region_s%child(1)%poly, region_s%child(2)%poly )
       if ( region_s%child(2)%uvbox(1) <= region_s%uvbox(1) + EPSregion ) then
          call subdiv_bezier2_only_v( &
               region_s%poly, &
               v=tuv_subdiv(3), &
               bs=region_s%child(1)%poly, &
               bn=region_s%child(2)%poly )
       else
          call subdiv_bezier2_only_u( &
               region_s%poly, &
               u=tuv_subdiv(2), &
               bw=region_s%child(1)%poly, &
               be=region_s%child(2)%poly )
       end if

    end if



    !! Carry on the recursion with the children
    if ( all(nchild < 2) ) then  
       PRINT *,' TBOX =',REGION_C%UVBOX
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
       CALL WRITE_POLYNOMIAL( REGION_C%POLY, 'dev_intersection_simple_surface/region_c_bezier.bern' )
       CALL WRITE_POLYNOMIAL( REGION_S%POLY, 'dev_intersection_simple_surface/region_s_bezier.bern' )
       IF ( ASSOCIATED( REGION_C%XYZBOX ) ) THEN
          CALL WRITE_OBB( REGION_C%XYZBOX, 'dev_intersection_simple_surface/xyzbox_c.dat' )
       ELSE
          PRINT *,'OBB_C N/A'
       END IF
       IF ( ASSOCIATED( REGION_S%XYZBOX ) ) THEN
          CALL WRITE_OBB( REGION_S%XYZBOX, 'dev_intersection_simple_surface/xyzbox_s.dat' )
       ELSE
          PRINT *,'OBB_S N/A'
       END IF
       !IF ( ASSOCIATED( REGION_C%XYZBOX ) .AND. ASSOCIATED( REGION_S%XYZBOX ) ) THEN  
       !END IF
       RETURN
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













  subroutine inherit_points( &
       region, &
       coords, &
       npts )
    implicit none
    type(type_region), intent(inout) :: region
    integer,           intent(in)    :: npts
    real(kind=fp),     intent(in)    :: coords(region%dim,npts)
    logical, allocatable             :: mask(:)
    integer                          :: idim, ipt, jpt

    if ( .not.associated(region%parent) ) return
    if ( region%parent%npts < 1 ) return

    allocate( mask(region%parent%npts) )
    mask(:) = .true.
    outer : do jpt = 1,region%parent%npts
       ipt = region%parent%ipts(jpt)
       do idim = 1,region%dim
          if ( .not.is_in_interval( &
               coords(idim,ipt), &
               region%uvbox(2*idim-1), &
               region%uvbox(2*idim) ) ) then
             mask(jpt) = .false.
             cycle outer
          end if
       end do
    end do outer

    call append_n( &
         region%ipts, &
         region%npts, &
         pack( region%parent%ipts(1:region%parent%npts), mask ), &
         count(mask), &
         unique=.true. )
    deallocate( mask )

  end subroutine inherit_points

















  recursive subroutine intersect_simple_surfaces( &
       surfroot, &
       region, &
       param_vector, &
       interdat, &
       uvxyz, &
       nuvxyz, &
       stat_degeneracy )
    implicit none
    integer, parameter                        :: nuvxyz_init = 10
    type(ptr_surface),          intent(in)    :: surfroot(2)
    type(ptr_region),           intent(inout) :: region(2)
    real(kind=fp),              intent(in)    :: param_vector(3)
    integer,                    intent(inout) :: stat_degeneracy
    !INTEGER,                    intent(inout) :: interdat
    type(type_listcurves),      intent(inout) :: interdat
    real(kind=fp), allocatable, intent(inout) :: uvxyz(:,:)
    integer,                    intent(inout) :: nuvxyz
    integer                                   :: stat_point(2)
    real(kind=fp)                             :: uv_subdiv(2)
    integer                                   :: stat_subdiv
    type(ptr_region)                          :: newregion(2)
    logical, allocatable                      :: mask(:,:)
    real(kind=fp), allocatable                :: newuvxyz(:,:)
    integer                                   :: nnewuvxyz
    integer                                   :: icurv, ivar, ival, ipt
    integer                                   :: ichild, jchild
    CHARACTER                                 :: STR1

    if ( stat_degeneracy > 1 ) return
    
    !PRINT *,'';PRINT *,'';PRINT *,'';
    !PRINT *,'UVBOXES ='
    !DO ISURF = 1,2
    !   PRINT *,REGION(ISURF)%PTR%UVBOX
    !END DO

    !Inherit from parent the points of interdat contained in the current regions
    !do isurf = 1,2
    !call inherit_points( &
    !     region(isurf)%ptr, &
    !     coords(1,1:npts), &
    !     npts )
    !end do

    ! intersect the 4 borders of each surface with the other surface
    !nuvxyz = 0
    !allocate( uvxyz(7,nuvxyz_init) )
    if ( .not.allocated(uvxyz) ) allocate( uvxyz(7,nuvxyz_init) )
    !PRINT *,'BEFORE, ',NUVXYZ,' POINTS'
    do icurv = 1,2 ! <--------------------------------------------------------+
       do ivar = 1,2 ! <---------------------------------------------------+  !
          do ival = 1,2 ! <---------------------------------------------+  !  !
             !                                                          !  !  !
             call intersect_border_surface( &                           !  !  !
                  surfroot, &                                           !  !  !
                  region, &                                             !  !  !
                  icurv, &                                              !  !  !
                  ivar, &                                               !  !  !
                  ival, &                                               !  !  !
                  uvxyz, &                                              !  !  !
                  nuvxyz, &                                             !  !  !
                  stat_degeneracy )                                     !  !  !
             !                                                          !  !  !
             if ( stat_degeneracy > 1 ) then ! <-------------------+    !  !  !
                ! propagation erreur vers routine appelante ...    !    !  !  !
                PRINT *,'STAT_DEGENERACY =',STAT_DEGENERACY
                return                                             !    !  !  !
             end if ! <--------------------------------------------+    !  !  !
             !                                                          !  !  !
          end do ! <----------------------------------------------------+  !  !
       end do ! <----------------------------------------------------------+  !
    end do ! <----------------------------------------------------------------+
    !PRINT *,'AFTER, ',NUVXYZ,' POINTS'


    if ( nuvxyz > 0 ) then ! <---------------------------------------------------------------+
       !PRINT *,NUVXYZ,' POINTS:'
       !CALL PRINT_MAT( TRANSPOSE( UVXYZ(:,1:NUVXYZ) ) )
       if ( nuvxyz > 2 ) then ! <------------------------------------------------------+     !
          ! subdivide both regions at their center point                               !     !
          do isurf = 1,2 ! <-------------------------------------+                     !     !
             !                                                   !                     !     !
             uv_subdiv = 0.5_fp * ( &                            !                     !     !
                  region(isurf)%ptr%uvbox([1,3]) + &             !                     !     !
                  region(isurf)%ptr%uvbox([2,4]) )               !                     !     !
             !                                                   !                     !     !
             call subdiv_region( &                               !                     !     !
                  region(isurf)%ptr, &                           !                     !     !
                  uv_subdiv, &                                   !                     !     !
                  stat_subdiv )                                  !                     !     !
             !                                                   !                     !     !
             if ( stat_subdiv > 0 ) then ! <-----+               !                     !     !
                stat_degeneracy = 33             !               !                     !     !
                PRINT *,'STAT_DEGENERACY =',STAT_DEGENERACY
                return                           !               !                     !     !
             end if ! <--------------------------+               !                     !     !
             !                                                   !                     !     !
             if ( stat_subdiv == 0 ) then ! ----------------+    !                     !     !
                allocate( &                                 !    !                     !     !
                     region(isurf)%ptr%child(1)%poly, &     !    !                     !     !
                     region(isurf)%ptr%child(2)%poly, &     !    !                     !     !
                     region(isurf)%ptr%child(3)%poly, &     !    !                     !     !
                     region(isurf)%ptr%child(4)%poly )      !    !                     !     !
                call subdiv_bezier2( &                      !    !                     !     ! 
                     region(isurf)%ptr%poly, &              !    !                     !     !
                     [0.5_fp, 0.5_fp], &                    !    !                     !     !
                     bsw=region(isurf)%ptr%child(1)%poly, & !    !                     !     !
                     bse=region(isurf)%ptr%child(2)%poly, & !    !                     !     !
                     bnw=region(isurf)%ptr%child(3)%poly, & !    !                     !     !
                     bne=region(isurf)%ptr%child(4)%poly )  !    !                     !     !
             end if ! <-------------------------------------+    !                     !     !
             !                                                   !                     !     !
          end do ! <---------------------------------------------+                     !     !
          !                                                                            !     !
          ! carry on the recursion with the 4*4 new pairs of regions                   !     !
          allocate( mask(nuvxyz,2) )                                                   !     !
          mask(:,:) = .true.                                                           !     !
          do jchild = 1,size(region(2)%ptr%child) ! <-----------------------------+    !     !
             newregion(2)%ptr => region(2)%ptr%child(jchild)                      !    !     ! 
             !                                                                    !    !     !     
             do ipt = 1,nuvxyz ! <---------------------------------------------+  !    !     !
                if ( .not.is_in_interval( &                                    !  !    !     !
                     uvxyz(3,ipt), &                                           !  !    !     !
                     newregion(2)%ptr%uvbox(1), &                              !  !    !     !
                     newregion(2)%ptr%uvbox(2) ) .or. &                        !  !    !     !
                     .not.is_in_interval( &                                    !  !    !     !
                     uvxyz(4,ipt), &                                           !  !    !     !
                     newregion(2)%ptr%uvbox(3), &                              !  !    !     !
                     newregion(2)%ptr%uvbox(4) ) ) mask(ipt,2) = .false.       !  !    !     !
             end do ! <--------------------------------------------------------+  !    !     !
             !PRINT *,'JCHILD =',JCHILD,', UVBOX=',NEWREGION(2)%PTR%UVBOX
             !PRINT *,'MASK =',MASK(:,2)
             !                                                                    !    !     !     
             do ichild = 1,size(region(1)%ptr%child) ! <-----------------------+  !    !     !    
                newregion(1)%ptr => region(1)%ptr%child(ichild)                !  !    !     !
                !                                                              !  !    !     !
                do ipt = 1,nuvxyz ! <---------------------------------------+  !  !    !     !
                   if ( .not.is_in_interval( &                              !  !  !    !     !
                        uvxyz(1,ipt), &                                     !  !  !    !     !
                        newregion(1)%ptr%uvbox(1), &                        !  !  !    !     !
                        newregion(1)%ptr%uvbox(2) ) .or. &                  !  !  !    !     !
                        .not.is_in_interval( &                              !  !  !    !     !
                        uvxyz(2,ipt), &                                     !  !  !    !     !
                        newregion(1)%ptr%uvbox(3), &                        !  !  !    !     !
                        newregion(1)%ptr%uvbox(4) ) ) mask(ipt,1) = .false. !  !  !    !     !
                end do ! <--------------------------------------------------+  !  !    !     !
                !PRINT *,'ICHILD =',ICHILD,', UVBOX=',NEWREGION(1)%PTR%UVBOX
                !PRINT *,'MASK =',MASK(:,1)
                !                                                              !  !    !     !
                nnewuvxyz = count( mask(:,1) .and. mask(:,2) )                 !  !    !     !
                !PRINT *,'NNEWUVXYZ =',NNEWUVXYZ
                if ( nnewuvxyz > 0 ) then ! <-------------------+              !  !    !     !
                   allocate( newuvxyz(7,nnewuvxyz) )            !              !  !    !     !
                   newuvxyz(1:7,1:nnewuvxyz) = reshape( &       !              !  !    !     !
                        pack( &                                 !              !  !    !     !
                        uvxyz(1:7,1:nuvxyz), &                  !              !  !    !     !
                        spread( mask(:,1).and.mask(:,2), &      !              !  !    !     !
                        dim=1, ncopies=7 ) ), [7,nnewuvxyz] )   !              !  !    !     !
                   !PRINT *,'NEWUVXYZ ='
                   !CALL PRINT_MAT( TRANSPOSE( NEWUVXYZ(:,1:NNEWUVXYZ) ) )
                end if ! <--------------------------------------+              !  !    !     !
                !                                                              !  !    !     !
                call intersect_simple_surfaces( &                              !  !    !     !
                     surfroot, &                                               !  !    !     !
                     newregion, &                                              !  !    !     !
                     param_vector, &                                           !  !    !     !
                     interdat, &                                               !  !    !     !
                     newuvxyz, &                                               !  !    !     !
                     nnewuvxyz, &                                              !  !    !     !
                     stat_degeneracy )                                         !  !    !     !
                if ( allocated(newuvxyz) ) deallocate( newuvxyz )              !  !    !     !
                !                                                              !  !    !     !  
             end do ! <--------------------------------------------------------+  !    !     !
          end do ! <--------------------------------------------------------------+    !     !
          deallocate( mask )                                                           !     ! 
          !                                                                            !     !
       else ! -------------------------------------------------------------------------+     !
          ! (nuvxyz <= 2)                                                              !     !
          ! classifier points : entering, exiting, isolated                            !     !
          call classify_border_surface_intersection_point( &                           !     !
               surfroot, &                                                             !     !
               region, &                                                               !     !
               reshape( uvxyz(1:4,1:nuvxyz), [2,2,nuvxyz] ), &                         !     !
               nuvxyz, &                                                               !     !
               stat_point(1:nuvxyz) )                                                  !     !
          !                                                                            !     !
          if ( nuvxyz == 2 ) then ! <----------------------------+                     !     !
             if ( product(stat_point) < 0 ) then ! <----------+  !                     !     !
                ! trace courbe...                             !  !                     !     !
                PRINT *,'1 CURVE :)'   
                CALL PRINT_MAT( TRANSPOSE(UVXYZ(:,1:2)) )
                interdat%nc = interdat%nc + 1
                interdat%curve(interdat%nc)%uvxyz = uvxyz(:,1:2)
             elseif ( all(stat_point == 0) ) then ! ----------+  !                     !     !
                ! 2 points isolés...                          !  !                     !     !
                PRINT *,'2 ISOLATED POINTS :)'
                CALL PRINT_MAT( TRANSPOSE(UVXYZ(:,1:2)) )
             else ! ------------------------------------------+  !                     !     !
                ! erreur                                      !  !                     !     !
                stat_degeneracy = 555                         !  !                     !     !
                PRINT *,'***********'
                PRINT *,'2 POINTS, INCORRECT CONFIGURATION ='
                CALL PRINT_MAT( TRANSPOSE(UVXYZ(:,1:2)) )
                PRINT *,'STAT_POINT =',STAT_POINT(1:2)
                DO ISURF = 1,2
                   WRITE (STR1,'(I1)') ISURF
                   CALL WRITE_POLYNOMIAL( REGION(ISURF)%PTR%POLY, 'dev_intersection_simple_surface/region_' // STR1 // '.bern' )
                END DO
                PRINT *,'***********'
                return                                        !  !                     !     !
             end if ! <---------------------------------------+  !                     !     !
          else ! ------------------------------------------------+                     !     !
             ! nuvxyz = 1                                        !                     !     !
             if ( stat_point(1) == 0 ) then ! <---------------+  !                     !     !
                ! 1 point isolé...                            !  !                     !     !
                PRINT *,'1 ISOLATED POINT :)'
                PRINT *, UVXYZ(:,1)
             else ! ------------------------------------------+  !                     !     !
                ! erreur                                      !  !                     !     !
                stat_degeneracy = 666                         !  !                     !     !
                PRINT *,'***********'
                PRINT *,'1 POINT, INCORRECT CONFIGURATION ='
                PRINT *, UVXYZ(:,1)
                PRINT *,'STAT_POINT =',STAT_POINT(1)
                PRINT *,'***********'
                return                                        !  !                     !     !
             end if ! <---------------------------------------+  !                     !     !
          end if ! <---------------------------------------------+                     !     ! 
          !                                                                            !     !
       end if ! <----------------------------------------------------------------------+     !
    end if ! <-------------------------------------------------------------------------------+

    if ( allocated(uvxyz) ) deallocate( uvxyz )

  end subroutine intersect_simple_surfaces







  subroutine classify_border_surface_intersection_point( &
      surf, &
      region, &
      uv, &
      npts, &
      stat_point ) 
    implicit none
    !stat_pointt = -1 : exiting
    !               0 : isolated
    !               1 : entering
    type(ptr_surface), intent(in)  :: surf(2)
    type(ptr_region),  intent(in)  :: region(2)
    integer,           intent(in)  :: npts
    real(kind=fp),     intent(in)  :: uv(2,2,npts)
    integer,           intent(out) :: stat_point(npts)
    real(kind=fp)                  :: uv_s(2,2,2), xyz_s(3,2)
    integer                        :: stat_contactpoint
    integer                        :: stat(2), sumstat
    integer                        :: isurf, ipt

    do ipt = 1,npts ! <---------------------------------------------+
       call tangent_to_intersection_curve( &                        !
            surf, &                                                 !
            uv(:,:,ipt), &                                          !
            uv_s, &                                                 !
            xyz_s, &                                                !
            stat_contactpoint )                                     !
       !                                                            !
       if ( stat_contactpoint > 0 ) then ! <---------------------+  !
          PRINT *,'classify_border_surface_intersection_point : stat_contactpoint = ',stat_contactpoint
          STOP
          if ( stat_contactpoint == 1 ) then ! <-------------+   !  !
             ! two tangent directions                        !   !  !
             ! (...)                                         !   !  !
          else ! --------------------------------------------+   !  !                                               
             ! undefined tangent direction                   !   !  !
             ! (...)                                         !   !  !
          end if ! <-----------------------------------------+   !  !
       end if ! <------------------------------------------------+  !
       !                                                            !
       do isurf = 1,2                                               !
          call classify_endpoint( &                                 !
               uv(:,isurf,ipt), &                                   !
               uv_s(:,1,isurf), &                                   !
               region(isurf)%ptr%uvbox, &                           !
               stat(isurf) )                                        !
       end do                                                       !
       !                                                            !
       sumstat = sum(stat)                                          !
       if ( sumstat == 0 ) then ! <---------+                       !
          stat_point(ipt) = 0               !                       !
       elseif ( sumstat < 0 ) then ! -------+                       !
          stat_point(ipt) = -1              !                       !
       else ! ------------------------------+                       !
          stat_point(ipt) = 1               !                       !
       end if ! <---------------------------+                       !
    end do ! <------------------------------------------------------+

  end subroutine classify_border_surface_intersection_point






  subroutine classify_endpoint( &
      uv, &
      tng, &
      uvbox, &
      stat )
    implicit none
    !stat = -1 : exiting
    !        0 : interior
    !        1 : entering
    real(kind=fp), intent(in)  :: uv(2)
    real(kind=fp), intent(in)  :: tng(2) ! tangential direction
    real(kind=fp), intent(in)  :: uvbox(4)
    integer,       intent(out) :: stat
    real(kind=fp)              :: uvloc
    integer                    :: ivar
    
    stat = 0
    do ivar = 1,2
       uvloc = ab2n1p1( uv(ivar), uvbox(2*ivar-1), uvbox(2*ivar) )
       if ( abs(uvloc) > 1._fp - EPSmath ) then
          if ( tng(ivar)*uvloc > 0._fp ) then
             stat = -1 ! exiting
             return
          else
             stat = 1 ! entering
          end if
       end if
    end do  

  end subroutine classify_endpoint











  subroutine tangent_to_intersection_curve( &
      surf, &
      uv, &
      uv_s, &
      xyz_s, &
      stat )
    implicit none
    ! ( see "Shape interrogation for computer aided design and manufacturing", &
    ! Patrikalakis et al. (2009), pp.166-175)
    type(ptr_surface), intent(in)  :: surf(2)
    real(kind=fp),     intent(in)  :: uv(2,2)
    real(kind=fp),     intent(out) :: uv_s(2,2,2)
    real(kind=fp),     intent(out) :: xyz_s(3,2)
    integer,           intent(out) :: stat
    real(kind=fp)                  :: xyz_uv(3,2,2)
    real(kind=fp)                  :: n(3,2), nsqr(2), aux(3), normxyz_s
    integer                        :: isurf, ivar, i

    ! compute normal vectors to both surfaces
    do isurf = 1,2 ! <--------------------------------------------------+
       do ivar = 1,2 ! <-------------------------+                      !
          call evald1( &                         !                      !
               xyz_uv(:,ivar,isurf), &           !                      !
               surf(isurf)%ptr, &                !                      !
               uv(:,isurf), &                    !                      !
               ivar )                            !                      !
       end do ! <--------------------------------+                      !
       n(:,isurf) = cross( xyz_uv(:,1,isurf), xyz_uv(:,2,isurf) )       !
    end do ! <----------------------------------------------------------+

    xyz_s(:,1) = cross( n(:,1), n(:,2) )
    normxyz_s = norm2( xyz_s(:,1) )
    if ( normxyz_s > EPSmath ) then ! <------------------------------------------------+
       ! the normals are not collinear, this is a transversal intersection point       !
       stat = 0                                                                        !
       xyz_s(:,1) = xyz_s(:,1) / norm2( xyz_s(:,1) )                                   !
    else ! ----------------------------------------------------------------------------+
       ! this is a tangential intersection point, the tangent direction(s) to          !
       ! the intersection banch(es) is (are) computed using 2nd derivatives            !
       call characterize_tangential_intersection_point( &                              !
            surf, &                                                                    !
            uv, &                                                                      !
            xyz_uv, &                                                                  !
            n, &                                                                       !
            xyz_s, &                                                                   !
            stat )                                                                     !
    end if ! <-------------------------------------------------------------------------+
    
    if ( stat < 2 ) then ! <---------------------------------------------------------------------+
       ! compute tangential direction(s) in the uv-space of each surface                         !
       nsqr = sum( n**2, 1 )                                                                     !
       do isurf = 1,2 ! loop over surfaces ! <---------------------------------------------+     !
          do i = 1,stat+1 ! loop over tangential directions ! <-------------------------+  !     ! 
             aux = cross( xyz_s(:,i), n(:,isurf) )                                      !  !     !
             do ivar = 1,2 ! loop over parameters u,v ! <----------------------------+  !  !     !
                uv_s(ivar,i,isurf) = real( (-1)**ivar, kind=fp) * &                  !  !  !     !
                     dot_product( xyz_uv(:,1+mod(ivar,2),isurf), aux ) / nsqr(isurf) !  !  !     !
             end do ! <--------------------------------------------------------------+  !  !     !
          end do ! <--------------------------------------------------------------------+  !     !
       end do ! <--------------------------------------------------------------------------+     !
    else ! --------------------------------------------------------------------------------------+
       ! the tangential direction is undefined                                                   !
       ! ...                                                                                     !
    end if! <------------------------------------------------------------------------------------+

  end subroutine tangent_to_intersection_curve



  


  
  subroutine characterize_tangential_intersection_point( &
       surf, &
       uv, &
       s1, &
       n, &
       xyz_s, &
       stat )
    ! ( see "Shape interrogation for computer aided design and manufacturing", &
    ! Patrikalakis et al. (2009), pp.172-175 )
    ! stat = 0 : single tangential intersection curve
    ! stat = 1 : two tangential intersection branches
    ! stat = 2 : isolated tangential contact point
    ! stat = 3 : higher-order contact point
    implicit none
    real(kind=fp), parameter       :: EPS = EPSmath
    type(ptr_surface), intent(in)  :: surf(2)
    real(kind=fp),     intent(in)  :: uv(2,2)
    real(kind=fp),     intent(in)  :: s1(3,2,2)
    real(kind=fp),     intent(in)  :: n(3,2)
    real(kind=fp),     intent(out) :: xyz_s(3,2)
    integer,           intent(out) :: stat
    real(kind=fp)                  :: s2(3)
    real(kind=fp)                  :: aux(3,2), A(2,2), LMN(3,2), lambda, B(3)
    real(kind=fp)                  :: discr, w(2)
    integer                        :: isurf, ivar

    ! compute A
    do ivar = 1,2
       aux(:,ivar) = cross( n(:,1), s1(:,1+mod(ivar,2),2) )
    end do
    aux(:,1) = -aux(:,1)
    A = matmul( transpose(aux), s1(:,:,1) )

    ! compute second fundamental forms
    do isurf = 1,2
       do ivar = 1,3
          call evald2( &
               s2, &
               surf(isurf)%ptr, &
               uv(:,isurf), &
               ivar )
          LMN(ivar,isurf) = dot_product( n(:,isurf), s2 )
       end do
    end do

    ! compute B = [ b11 ; b12 ; b22 ]
    B(1) = LMN(1,2)*A(1,1)**2 + &
         LMN(2,2)*2*A(1,1)*A(2,1) + &
         LMN(3,2)*A(2,1)**2

    B(2) = LMN(1,2)*A(1,1)*A(1,2) + &
         LMN(2,2)*(A(1,1)*A(2,2) + A(1,2)*A(2,1)) + &
         LMN(3,2)*A(2,1)*A(2,2)

    B(3) = LMN(1,2)*A(1,2)**2 + &
         LMN(2,2)*2*A(1,2)*A(2,2) + &
         LMN(3,2)*A(2,2)**2

    lambda = 1._fp / &
         ( sqrt( product(sum(n**2, 1)) ) * dot_product( n(:,1), n(:,2) ) )

    B = LMN(:,1) - lambda * B

    ! analyze B
    discr = B(2)**2 - B(1)*B(3)

    if ( all( abs(B) < EPS ) ) then
       ! higher order contact point (xyz_s cannot be determined)
       stat = 3
       return

    elseif ( abs( discr ) < EPS ) then
       ! tangential intersection curve (one single xyz_s)
       stat = 0
       if ( abs(B(1)) > abs(B(3)) ) then
          w(1) = -B(2) / B(1)
          xyz_s(:,1) = w(1)*s1(:,1,1) + s1(:,2,1)
       else
          w(1) = -B(2) / B(3)
          xyz_s(:,1) = s1(:,1,1) + w(1)*s1(:,2,1)
       end if

    elseif ( discr < 0._fp ) then
       ! isolated tangential contact point (xyz_s undefined)
       stat = 2
       return

    else
       ! branch point (two disctinct xyz_s)
       stat = 1
       w = real( [-1,1], kind=fp ) * sqrt( discr )
       if ( abs(B(1)) > abs(B(3)) ) then
          w = (w - B(2)) / B(1)
          xyz_s = outer_product( s1(:,1,1), w ) + &
               spread( s1(:,2,1), dim=2, ncopies=2 )
       else
          w = (w - B(2)) / B(3)
          xyz_s = spread( s1(:,1,1), dim=2, ncopies=2 ) + &
               outer_product( s1(:,2,1), w )
       end if

    end if

    ! unitize xyz_s
    do ivar = 1,stat+1
       xyz_s(:,ivar) = xyz_s(:,ivar) / norm2( xyz_s(:,ivar) )
    end do

  end subroutine characterize_tangential_intersection_point




end program dev_intersection_simple_surface
