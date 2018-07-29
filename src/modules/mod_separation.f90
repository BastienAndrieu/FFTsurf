module mod_separation

  USE MOD_UTIL
  use mod_math
  use mod_tolerances

  implicit none

  logical,       parameter       :: OPTcartsp  = .true.  ! cartesian separating plane
  logical,       parameter       :: OPTsphbbsp = .true.  ! spherical bounding box separation
  logical,       parameter       :: OPTlpsp    = .true.  ! linear programming separation
  real(kind=fp), parameter       :: EPSlpsep = real( 1.e-6, kind=fp )  

  type(type_matrix), allocatable :: ch2be_matrices(:)

contains

  subroutine separating_plane( &
       xyz1, &
       xyz2, &
       n1, &
       n2, &
       vec, &
       separable )
    ! Searches for a plane that strictly sepatates two sets of points 
    ! of coords. 'xyz1' and 'xyz2' by using different methods which
    ! can be switched on/off by the opt_* global parameters.
    ! Returns a unit normal vector 'vec' of such a plane (if any).
    implicit none
    logical, parameter          :: VERBOSE = .false.
    integer,       intent(in)   :: n1, n2
    real(kind=fp), intent(in)   :: xyz1(n1,3)
    real(kind=fp), intent(in)   :: xyz2(n2,3)
    real(kind=fp), intent(out)  :: vec(3)
    logical,       intent(out)  :: separable

    if ( OPTcartsp ) then
       ! Least expensive, least effective method 
       call separating_cartesian_plane( &
            xyz1, &
            xyz2, &
            n1, &
            n2, &
            vec, &
            separable )

       if ( separable ) then
          IF ( VERBOSE ) PRINT *,'CARTESIAN PLANE'
          return
       end if

    end if


    if ( OPTsphbbsp ) then
       ! More expensive, more effective method
       call separate_spherical_bounding_boxes( &
            xyz1, &
            xyz2, &
            n1, &
            n2, &
            vec, &
            separable )

       if ( separable ) then
          IF ( VERBOSE ) PRINT *,'SPHERICAL BOUNDING BOX'
          return
       end if

    end if


    if ( OPTlpsp ) then
       ! Most expensive, most effective method
       call separation_linearprogramming( &
            xyz1, &
            xyz2, &
            n1, &
            n2, &
            vec, &
            separable )

       if ( separable ) then
          IF ( VERBOSE ) PRINT *,'LINEAR PROGRAMMING'
          return
       end if

    end if

  end subroutine separating_plane













  subroutine separating_cartesian_plane( &
       xyz1, &
       xyz2, &
       n1, &
       n2, &
       vec, &
       separable )
    ! Tests if one of the "cartesian" planes (x=0), (y=0), (z=0) 
    ! strictly sepatates two sets of points of coords. 'xyz1' and 'xyz2'. 
    ! Returns a unit normal vector 'vec' of such a plane (if any).
    implicit none
    integer,       intent(in)  :: n1, n2
    real(kind=fp), intent(in)  :: xyz1(n1,3)
    real(kind=fp), intent(in)  :: xyz2(n2,3)
    real(kind=fp), intent(out) :: vec(3)
    logical,       intent(out) :: separable
    real(kind=fp)              :: m(2)
    integer                    :: iaxe

    separable = .false.

    do iaxe = 1,3
       m(1) = max( minval( xyz1(:,iaxe) ), &
            minval( xyz2(:,iaxe) ) ) ! max(min)
       m(2) = min( maxval( xyz1(:,iaxe) ), &
            maxval( xyz2(:,iaxe) ) ) ! min(max)

       if ( m(2) < -EPSfp .and. m(1) > EPSfp ) then
          ! separating cartesian plane found 
          vec(:) = 0._fp
          vec(iaxe) = 1._fp
          separable = .true.
          return
       end if
    end do

  end subroutine separating_cartesian_plane



  subroutine separate_spherical_bounding_boxes( &
       xyz1, &
       xyz2, &
       n1, &
       n2, &
       vec, &
       separable, &
       mask_axes )
    ! Tests if a plane containing one of the three cartesian axes
    ! strictly sepatates two sets of points of coords. 'xyz1' and 'xyz2'. 
    ! Returns a unit normal vector 'vec' of such a plane (if any).
    implicit none
    LOGICAL, PARAMETER :: DEBUG = .false.
    integer,           intent(in)  :: n1, n2
    real(kind=fp),     intent(in)  :: xyz1(n1,3)
    real(kind=fp),     intent(in)  :: xyz2(n2,3)
    real(kind=fp),     intent(out) :: vec(3)
    logical,           intent(out) :: separable
    logical, optional, intent(in)  :: mask_axes(3)
    real(kind=fp)                  :: wedge(2,2), d, alpha, angv
    logical                        :: bounded
    integer                        :: iaxe, iset

    separable = .false.

    IF (DEBUG) THEN
       if ( present(mask_axes) ) print *,'MASK AXES =',mask_axes
    END IF

    loop_axes : do iaxe = 1,3
       
       if ( present( mask_axes ) ) then
          if ( .not.mask_axes(iaxe) ) cycle
       end if

       IF (DEBUG) THEN
          PRINT *,''
          PRINT *,'IAXE =',IAXE
       END IF

       call minimal_bounding_sector( &
            xyz1(1:n1,[1+mod(iaxe,3), 1+mod(iaxe+1,3)]), &
            n1, &
            wedge(:,1), &
            bounded )
       if ( .not.bounded ) cycle loop_axes

       call minimal_bounding_sector( &
            xyz2(1:n2,[1+mod(iaxe,3), 1+mod(iaxe+1,3)]), &
            n2, &
            wedge(:,2), &
            bounded )
       if ( .not.bounded ) cycle loop_axes

       do iset = 1,2
          if ( wedge(2,iset) >= 0.5_fp * CSTpi - EPSmath ) then
             IF (DEBUG) PRINT *,'WEDGE',ISET,' WIDER THAN PI'
             cycle loop_axes
          end if
       end do

       IF (DEBUG) THEN
          PRINT *,'BOUNDING SECTORS ='
          CALL PRINT_MAT( WEDGE )
       END IF

       d = diff_angle( wedge(1,1), wedge(1,2) )
       alpha = 0.5_fp * ( abs(d) - sum(wedge(2,:)) )
       IF (DEBUG) PRINT *,'ALPHA',ALPHA
       if ( alpha > EPSfp ) then
          if ( alpha + 2._fp*wedge(2,1) > CSTpi + EPSfp ) then
             angv = wedge(1,1) + 0.5_fp * CSTpi
          elseif ( alpha + 2._fp*wedge(2,2) > CSTpi + EPSfp ) then
             angv = wedge(1,2) + 0.5_fp * CSTpi
          else
             if ( d > EPSfp ) then
                iset = 1
             else
                iset = 2
             end if
             angv = mean_angle( wedge(1,[iset,1+mod(iset,2)]) + &
                  [-1._fp, 1._fp]*wedge(2,[iset,1+mod(iset,2)]) )
          end if

          vec(:) = 0._fp
          vec(1+mod(iaxe,3)) = -sin( angv )
          vec(1+mod(iaxe+1,3)) = cos( angv )
          separable = .true.
          
          IF (DEBUG) THEN
             PRINT *,'SEPARABLE AROUND AXIS #', IAXE

             CALL WRITE_MATRIX( &
                  xyz1(1:n1,[1+mod(iaxe,3), 1+mod(iaxe+1,3)]), &
                  N1, &
                  2, &
                  '/stck/bandrieu/Bureau/CYPRES/Intersections/separation/xy1.dat' )
             CALL WRITE_MATRIX( &
                  xyz2(1:n2,[1+mod(iaxe,3), 1+mod(iaxe+1,3)]), &
                  N2, &
                  2, &
                  '/stck/bandrieu/Bureau/CYPRES/Intersections/separation/xy2.dat' )

             CALL WRITE_MATRIX( &
                  WEDGE, &
                  2, &
                  2, &
                  '/stck/bandrieu/Bureau/CYPRES/Intersections/separation/wedges.dat' )
          END IF
          return

       end if

    end do loop_axes

  end subroutine separate_spherical_bounding_boxes






  subroutine minimal_bounding_sector( &
       xy, &
       n, &
       wedge, &
       bounded )
    use mod_geometry
    implicit none
    LOGICAL, PARAMETER :: DEBUG = .false.
    integer,       intent(in)   :: n
    real(kind=fp), intent(in)   :: xy(n,2)
    real(kind=fp), intent(out)  :: wedge(2)
    logical,       intent(out)  :: bounded
    integer                     :: hull(n), nhull
    real(kind=fp), dimension(2) :: xy1, xy2, rng
    real(kind=fp)               :: ang(n), mnmx(2), mnpmxn(2)
    integer                     :: i


    rng = maxval( xy, dim=1 ) - minval( xy, dim=1 )
    IF (DEBUG) PRINT *,'RANGES =',RNG

    where( rng < EPSfp ) rng = 1._fp

    call quickhull( &
         matmul( diag(1._fp/rng), transpose(xy) ), &
         n, &
         hull, &
         nhull )
    
    IF (DEBUG) PRINT *,'HULL =',HULL(1:NHULL)

    bounded = .false.
    ! Check whether the origin is inside the convex hull :
    ! A point is inside a convex polygon if it lies on the left of each of 
    ! the counter-clockwise directed edges.
    ! The vector ni = (y_i - y_{i+1} , x_{i+1} - x_i) is normal to the edge e_i = [p_i, p_{i+1}]
    ! and points towards the interior of the polygon.
    ! The point z = (x,y) is strictly on the right of e_i if ( z - p_i ).n_i < 0.
    ! If z = (0,0) this reduces to x_i*(y_i - y_{i+1} + y_i*( x_{i+1} - x_i) > 0.
    if ( nhull == 1 ) then
       bounded = ( sum( xy(hull(1),:)**2 ) > 2._fp*EPSfp )
    elseif ( nhull == 2 ) then
       bounded = ( abs(distance_from_line( [0._fp, 0._fp], xy(hull(1),:), xy(hull(2),:) )) > EPSfp .or. &
            dot_product(xy(hull(1),:), xy(hull(2),:)) > EPSfp )
    else
       check_boundedness : do i = 1,nhull
          xy1 = xy(hull(i),:)
          xy2 = xy(hull(1+mod(i,nhull)),:)

          if ( xy1(1) * ( xy1(2) - xy2(2) ) + &
               xy1(2) * ( xy2(1) - xy1(1) ) > EPSfp ) then
             bounded = .true.
             exit check_boundedness
          end if
       end do check_boundedness
    end if

    IF (DEBUG) PRINT *,'BOUNDED?',BOUNDED

    if ( .not.bounded ) return

    ang(1:nhull) = atan2( xy(hull(1:nhull),2), xy(hull(1:nhull),1) )
    
    mnmx = [ minval(ang(1:nhull)), maxval(ang(1:nhull)) ]
    if ( mnmx(1) < -EPSfp .and. mnmx(2) > EPSfp ) then
       mnpmxn(1) = minval( ang(1:nhull), MASK=(ang(1:nhull) > 0._fp) )
       mnpmxn(2) = maxval( ang(1:nhull), MASK=(ang(1:nhull) < 0._fp) )
       if ( diff_angle( mnpmxn(2), mnpmxn(1) ) > 0._fp ) mnmx = mnpmxn
    end if

    wedge(1) = mean_angle( mnmx )
    wedge(2) = diff_angle( mnmx(2), wedge(1) )
    IF (DEBUG) PRINT *,'WEDGE =',WEDGE

  end subroutine minimal_bounding_sector











  subroutine separation_linearprogramming( &
       xyz1, &
       xyz2, &
       n1, &
       n2, &
       vec, &
       separable )
    use mod_linprog
    ! Searches for a plane that strictly sepatates two sets of points 
    ! of coords. 'xyz1' and 'xyz2' by solving a linear programming problem.
    ! Returns a unit normal vector 'vec' of such a plane (if any).
    implicit none
    integer,       intent(in)  :: n1, n2
    real(kind=fp), intent(in)  :: xyz1(n1,3)
    real(kind=fp), intent(in)  :: xyz2(n2,3)
    real(kind=fp), intent(out) :: vec(3)
    logical,       intent(out) :: separable
    integer                    :: stat
    integer                    :: perm(n1+n2)
    real(kind=fp)              :: LPmat( 6 + n1 + n2, 4 )
    integer                    :: i, l

    ! assemble the matrix of linear constraints
    LPmat(:,:) = 0._fp

    ! secondary constraints
    do i = 1,3
       LPmat(i,i) = 1._fp
    end do
    do i = 1,3
       LPmat(3+i,i) = -1._fp
    end do
    l = 6! l = 3
    LPmat(1:l,4) = 1._fp

    ! primary constraints
    LPmat( (l+1)   :(l+n1)   , 1:3 ) =  xyz1
    LPmat( (l+n1+1):(l+n1+n2), 1:3 ) = -xyz2

    LPmat(l+1:,4) = -EPSlpsep

    ! permute randomly the LP constraints
    call randperm( perm, n1+n2 )
    LPmat((l+1):(l+n1+n2),:) = LPmat(l+perm,:)

    !PRINT *,'LP_MAT ='
    !CALL PRINT_MAT( LPMAT(1:(l+n1+n2),1:4) )
    !PRINT *,''
    !CALL WRITE_MATRIX( LPMAT(1:(l+n1+n2),1:4), l+n1+n2, 4, 'linearprogramming/lpmat.dat' )

    ! solve the linear programming problem
    vec = real( [1, 1, 1], kind=fp ) / sqrt( 3._fp )
    call lp_solve( &
         vec, &
         stat, &
         LPmat(1:(l+n1+n2),1:4), &
         [0._fp, 0._fp, 0._fp], &
         3, &
         l+n1+n2 )

    separable = ( stat == 0 )
    !IF ( SEPARABLE ) THEN
    !   PRINT *,'VEC = ',VEC
    !   PRINT *,'LPRES =',MINVAL( MATMUL( LPMAT(:,1:3), VEC ) - LPMAT(:,4) )
    !END IF

    if ( stat < 0 ) PRINT *,'separation_linearprogramming : problème réalisable non borné'

  end subroutine separation_linearprogramming







  subroutine bounding_cone( &
       xyz, &
       axe, &
       angle )
    implicit none
    real(kind=fp), intent(in)  :: xyz(:,:)
    real(kind=fp), intent(out) :: axe(3)
    real(kind=fp), intent(out) :: angle
    real(kind=fp)              :: cosi
    integer                    :: i

    angle = 0._fp
    axe = sum( xyz, 1 ) / real( size(xyz,1), kind=fp )
    axe = axe / norm2( axe )

    do i = 1,size(xyz,1)
       cosi = dot_product( axe, xyz(i,:) ) / norm2( xyz(i,:) )
       angle = max( angle, acos( min( 1._fp, max( -1._fp, cosi ) ) ) )
    end do

  end subroutine bounding_cone




end module mod_separation
