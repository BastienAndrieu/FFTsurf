program dev_intersection_simple_surface

  use mod_util
  use mod_math
  use mod_regiontree
  use mod_polynomial
  use mod_diffgeom2
  use mod_separation  

  ! =================================================================================
  logical, parameter         :: ECONOMIZE = .true.

  integer                    :: narg, numtest, icurv, ivar, ival
  character(100)             :: arg
  character                  :: strnum
  character(2)               :: strnum2
  integer                    :: fileunit

  type(type_surface), target :: surf(2)
  type(type_region), target  :: root(2)
  type(ptr_surface)          :: surfroot(2)
  type(ptr_region)           :: region(2)

  real(kind=fp), allocatable :: uvxyz(:,:)
  integer                    :: nuvxyz
  integer                    :: stat_degeneracy
  integer                    :: i, isurf!, ipt
  integer*8                  :: tic, toc, count_rate
  ! =================================================================================


  ! Lecture arguments
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
     call read_polynomial( &
          surf(isurf)%x, &
                                !'/stck/bandrieu/Bureau/coeffstest/C' // strnum // '_test' // strnum2 // '.txt', &
          '/home/bastien/Bureau/coeffstest/C' // strnum // '_test' // strnum2 // '.txt', &
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
  PRINT *,''; PRINT *,''; PRINT *,''
  PRINT *,'ELAPSED =',REAL( TOC - TIC ) / REAL( COUNT_RATE )


  PRINT *,NUVXYZ,' INTERSECTION POINT(S)'
  call get_free_unit( fileunit )
  open( unit=fileunit, file='dev_intersection_simple_surface/uv_xyz.dat', action='write' )
  if ( nuvxyz < 1 ) then
     write (fileunit,*) ''
  else
     do ipt = 1,nuvxyz
        PRINT *,uvxyz(:,ipt)
        write (fileunit,*) uvxyz(:,ipt)
     end do
  end if
  close(fileunit)







  do isurf = 1,2
     call free_polynomial( region(isurf)%ptr%poly )
     deallocate( region(isurf)%ptr%poly )

     call free_polynomial( surf(isurf)%x )
     call free_polynomial( surf(isurf)%xu )
     call free_polynomial( surf(isurf)%xv )
     call free_polynomial( surf(isurf)%xuu )
     call free_polynomial( surf(isurf)%xuv )
     call free_polynomial( surf(isurf)%xvv )

     deallocate( root(isurf)%uvbox )
  end do


  if ( allocated( ch2be_matrices ) ) then
     do i = 1,size(ch2be_matrices)
        if ( allocated(ch2be_matrices(i)%mat) ) deallocate( ch2be_matrices(i)%mat )
     end do
     deallocate( ch2be_matrices )
  end if

  deallocate( uvxyz )





contains




  subroutine cheb2bern_poly( c, b )
    implicit none
    type(type_polynomial), intent(in)  :: c
    type(type_polynomial), intent(out) :: b
    real(kind=fp)                      :: au(c%degr(1)+1,c%degr(1)+1)
    real(kind=fp)                      :: av(c%degr(2)+1,c%degr(2)+1)
    integer                            :: i

    if ( c%base /= 1 ) STOP 'cheb2bern_poly : input polynomial not in Chebyshev basis'

    call reset_polynomial( poly=b, nvar=c%nvar, base=2, degr=c%degr, dim=c%dim )

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
    use mod_polynomial
    use mod_diffgeom2
    use mod_regiontree
    use mod_tolerances
    implicit none
    integer, parameter                        :: npts_init = 10
    type(ptr_surface),          intent(in)    :: root_s(2)
    type(ptr_region),           intent(inout) :: region(2)
    integer,                    intent(in)    :: icurv, ivar, ival
    real(kind=fp), allocatable, intent(inout) :: uvxyz(:,:)
    integer,                    intent(inout) :: nuvxyz
    integer,                    intent(inout) :: stat_degeneracy
    type(type_curve)                          :: root_c
    type(type_region)                         :: region_c
    type(type_region)                         :: region_s
    real(kind=fp), allocatable                :: tuvxyz(:,:)
    integer                                   :: ntuvxyz
    real(kind=fp)                             :: uv(2,2)
    integer                                   :: isurf, ipt, jpt

    PRINT *,''; PRINT *,''; PRINT *,''; PRINT *,''
    PRINT *,'ICURV, IVAR, IVAL =', ICURV, IVAR, IVAL

    isurf = 1 + mod(icurv,2)

    ! convert surface border to parametric curve
    call bivar2univar( &
         root_s(icurv)%ptr%x, &
         root_c%x, &
         ivar, &
         ival )

    call economize1( root_c%x, EPSmath )
    call compute_deriv1( root_c )
    call compute_deriv2( root_c )

    call write_polynomial( root_c%x, 'dev_intersection_simple_surface/root_c_x.cheb' )
    call write_polynomial( root_c%xt, 'dev_intersection_simple_surface/root_c_xt.cheb' )

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


    CALL WRITE_POLYNOMIAL( REGION_C%POLY, 'dev_intersection_simple_surface/root_c_bezier.bern' )
    CALL WRITE_POLYNOMIAL( REGION_S%POLY, 'dev_intersection_simple_surface/root_s_bezier.bern' )


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

       uv(ivar,icurv) = region(icurv)%ptr%uvbox(2*(ivar-1)+ival) !region(icurv)%ptr%uvbox(ival,ivar)
       uv(1+mod(ivar,2),icurv) = n1p12ab( &
            tuvxyz(1,ipt), &
            region(icurv)%ptr%uvbox(2*(1+mod(ivar,2))-1), &!region(icurv)%ptr%uvbox(1,1+mod(ivar,2)), &
            region(icurv)%ptr%uvbox(2*(1+mod(ivar,2))) )!region(icurv)%ptr%uvbox(2,1+mod(ivar,2)) )
       uv(:,isurf) = tuvxyz(2:3,ipt)
       !PRINT *,uv,tuvxyz(4:6,ipt)

       call append_vector( &
            [ uv(:,1), uv(:,2), tuvxyz(4:6,ipt) ], &
            7, &
            uvxyz, &
            nuvxyz )

    end do outer

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


































  ! METTRE LES DEUX REGIONS DANS UNE SEULE VARIABLE (TABLEAU DE PTR_REGION)
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
    logical, allocatable                              :: mask(:)
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

       if ( .not.overlap ) then
          !PRINT *,'DISJOINT OBBs'
          return
       end if


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





  
  

end program dev_intersection_simple_surface
