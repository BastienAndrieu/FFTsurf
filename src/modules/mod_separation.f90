module mod_separation

  use mod_math

  implicit none

  logical,           parameter   :: OPTcartsp  = .true.  ! cartesian separating plane
  logical,           parameter   :: OPTsphbbsp = .true.  ! spherical bounding box separation
  logical,           parameter   :: OPTlpsp    = .true.  ! linear programming separation
  real(kind=MATHpr), parameter   :: EPSlpsep = real( 1.e-6, kind=MATHpr )  

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
    logical, parameter              :: VERBOSE = .false.
    integer,           intent(in)   :: n1, n2
    real(kind=MATHpr), intent(in)   :: xyz1(n1,3)
    real(kind=MATHpr), intent(in)   :: xyz2(n2,3)
    real(kind=MATHpr), intent(out)  :: vec(3)
    logical,           intent(out)  :: separable

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
    integer,           intent(in)  :: n1, n2
    real(kind=MATHpr), intent(in)  :: xyz1(n1,3)
    real(kind=MATHpr), intent(in)  :: xyz2(n2,3)
    real(kind=MATHpr), intent(out) :: vec(3)
    logical,           intent(out) :: separable
    real(kind=MATHpr)              :: m(2)
    integer                        :: iaxe

    separable = .false.
    !PRINT *,'SIZE(XYZ1)=',SIZE(XYZ1,1),SIZE(XYZ1,2)
    !PRINT *,'SIZE(XYZ2)=',SIZE(XYZ2,1),SIZE(XYZ2,2)

    do iaxe = 1,3
       m(1) = max( minval( xyz1(:,iaxe) ), &
            minval( xyz2(:,iaxe) ) ) ! max(min)
       m(2) = min( maxval( xyz1(:,iaxe) ), &
            maxval( xyz2(:,iaxe) ) ) ! min(max)

       if ( ( m(2) < m(1) ) .and. ( m(2)*m(1) < 0._MATHpr ) ) then
          ! separating cartesian plane found 
          vec(:) = 0._MATHpr
          vec(iaxe) = 1._MATHpr
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
    integer,           intent(in)  :: n1, n2
    real(kind=MATHpr), intent(in)  :: xyz1(n1,3)
    real(kind=MATHpr), intent(in)  :: xyz2(n2,3)
    real(kind=MATHpr), intent(out) :: vec(3)
    logical,           intent(out) :: separable
    logical, optional, intent(in)  :: mask_axes(3)
    real(kind=MATHpr)              :: wedge(2,2), d, alpha, angv
    logical                        :: bounded
    integer                        :: iaxe, iset

    separable = .false.

    loop_axes : do iaxe = 1,3
       if ( present( mask_axes ) ) then
          if ( .not.mask_axes(iaxe) ) cycle
       end if

       IF (.false.) THEN
          call minimal_2d_wedge( &
               xyz1(:,[1+mod(iaxe,3), 1+mod(iaxe+1,3)]), &
               n1, &
               wedge(:,1) )
          call minimal_2d_wedge( &
               xyz2(:,[1+mod(iaxe,3), 1+mod(iaxe+1,3)]), &
               n2, &
               wedge(:,2) )
          do iset = 1,2
             if ( wedge(2,iset) >= 0.5_MATHpr * MATHpi ) then
                !PRINT *,'WEDGE',ISET,' WIDER THAN PI'
                cycle loop_axes
             end if
          end do
       ELSE
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
       END IF


       d = diff_angle( wedge(1,1), wedge(1,2) )
       alpha = 0.5_MATHpr * ( abs(d) - sum(wedge(2,:)) )
       !PRINT *,'ALPHA',REAL(ALPHA)
       if ( alpha > MATHeps ) then
          if ( alpha + 2._MATHpr*wedge(2,1) > MATHpi ) then
             angv = wedge(1,1) + 0.5_MATHpr * MATHpi
          elseif ( alpha + 2._MATHpr*wedge(2,2) > MATHpi ) then
             angv = wedge(1,2) + 0.5_MATHpr * MATHpi
          else
             if ( d > 0._MATHpr ) then
                iset = 1
             else
                iset = 2
             end if
             angv = mean_angle( wedge(1,[iset,1+mod(iset,2)]) + &
                  [-1._MATHpr, 1._MATHpr]*wedge(2,[iset,1+mod(iset,2)]) )
          end if

          vec(:) = 0._MATHpr
          vec(1+mod(iaxe,3)) = -sin( angv )
          vec(1+mod(iaxe+1,3)) = cos( angv )
          separable = .true.
          return

       end if

    end do loop_axes

  end subroutine separate_spherical_bounding_boxes






  subroutine minimal_2d_wedge( &
       xy, &
       n, &
       wedge )
    ! Computes a minimal wedge bounding a two-dimensional set of points 'xy'.
    ! The wedge is the set of points between two rays originating from (0,0)
    ! wedge(1) = angle between the x-axis and the bissector of the two rays;
    ! wedge(2) = half angle spanned by the wedge.
    ! (all angles are measured in radians)
    implicit none
    integer,           intent(in)  :: n
    real(kind=MATHpr), intent(in)  :: xy(n,2)
    real(kind=MATHpr), intent(out) :: wedge(2)
    real(kind=MATHpr)              :: ang(size(xy,1)), meanang, e(2), d
    integer                        :: i

    ang = atan2( xy(:,2), xy(:,1) )
    meanang = mean_angle( ang )

    e(:) = 0._MATHpr
    do i = 1,n
       d = diff_angle( ang(i), meanang )
       e(1) = min( e(1), d )
       e(2) = max( e(2), d )
    end do

    wedge(2) = 0.5_MATHpr * ( e(2) - e(1) )
    wedge(1) = mean_angle( meanang + e )
    if ( e(2) - e(1) > MATHpi ) then
       if ( wedge(1) > 0._MATHpr ) then
          wedge(1) = wedge(1) + MATHpi
       else
          wedge(1) = wedge(1) - MATHpi
       end if
    end if

  end subroutine minimal_2d_wedge





  subroutine minimal_bounding_sector( &
       xy, &
       n, &
       wedge, &
       bounded )
    use mod_geometry
    implicit none
    integer,           intent(in)   :: n
    real(kind=MATHpr), intent(in)   :: xy(n,2)
    real(kind=MATHpr), intent(out)  :: wedge(2)
    logical,           intent(out)  :: bounded
    integer                         :: hull(n), nhull
    real(kind=MATHpr), dimension(2) :: xy1, xy2
    real(kind=MATHpr)               :: ang(n), mnmx(2), mnpmxn(2)
    integer                         :: i

    call convex_hull_2d( &
         transpose(xy), &
         n, &
         hull, &
         nhull )
    
    bounded = .false.
    check_boundedness : do i = 1,nhull
       xy1 = xy(hull(i),:)
       xy2 = xy(hull(1+mod(i,nhull)),:)
       if ( xy1(1) * ( xy2(2) - xy1(2) ) + &
            xy1(2) * ( xy1(1) - xy2(1) ) < -MATHeps ) then
          bounded = .true.
          exit check_boundedness
       end if
    end do check_boundedness

    if ( .not.bounded ) return

    ang(1:nhull) = atan2( xy(hull(1:nhull),2), xy(hull(1:nhull),1) )
    
    !if ( maxval( xy(hull(1:nhull),1) ) < 0._MATHpr .and. &
    !     maxval( xy(hull(1:nhull),2) ) * minval( xy(hull(1:nhull),2) ) < 0._MATHpr ) then
    !   mnmx(1) = minval( ang(1:nhull), MASK=(ang(1:nhull) > 0) )
    !   mnmx(2) = maxval( ang(1:nhull), MASK=(ang(1:nhull) < 0) )
    !else
    !   mnmx = [ minval(ang(1:nhull)), maxval(ang(1:nhull)) ]
    !end if
    mnmx = [ minval(ang(1:nhull)), maxval(ang(1:nhull)) ]
    if ( mnmx(1) < 0._MATHpr .and. mnmx(2) > 0._MATHpr ) then
       mnpmxn(1) = minval( ang(1:nhull), MASK=(ang(1:nhull) > 0) )
       mnpmxn(2) = maxval( ang(1:nhull), MASK=(ang(1:nhull) < 0) )
       if ( diff_angle( mnpmxn(2), mnpmxn(1) ) > 0._MATHpr ) mnmx = mnpmxn
    end if

    wedge(1) = mean_angle( mnmx )
    wedge(2) = diff_angle( mnmx(2), wedge(1) )

  end subroutine minimal_bounding_sector











  subroutine separation_linearprogramming( &
       xyz1, &
       xyz2, &
       n1, &
       n2, &
       vec, &
       separable )
    use mod_linearprogramming
    ! Searches for a plane that strictly sepatates two sets of points 
    ! of coords. 'xyz1' and 'xyz2' by solving a linear programming problem.
    ! Returns a unit normal vector 'vec' of such a plane (if any).
    implicit none
    integer,           intent(in)  :: n1, n2
    real(kind=MATHpr), intent(in)  :: xyz1(n1,3)
    real(kind=MATHpr), intent(in)  :: xyz2(n2,3)
    real(kind=MATHpr), intent(out) :: vec(3)
    logical,           intent(out) :: separable
    integer                        :: stat
    real(kind=MATHpr)              :: LPmat( 3 + n1 + n2, 4 )
    integer                        :: i

    ! assemble the matrix of linear constraints
    LPmat(:,:) = 0._MATHpr

    ! secondary constraints
    do i = 1,3
       LPmat(i,i) = 1._MATHpr
    end do
    LPmat(1:3,4) = -1._MATHpr

    ! primary constraints
    LPmat(4:3+n1,1:3) = -xyz1
    LPmat(4+n1:3+n1+n2,1:3) = xyz2

    LPmat(4:,4) = EPSlpsep

    ! solve the linear programming problem
    call lpsolve( &
         vec, &
         stat, &
         LPmat, &
         real( [0, 0, 0], kind=MATHpr ) )

    separable = ( stat == 0 )

    if ( stat < 0 ) then
       PRINT *,'separation_linearprogramming : problème réalisable non borné'
       !STOP !?!
    end if

  end subroutine separation_linearprogramming







  subroutine bounding_cone( &
       xyz, &
       axe, &
       angle )
    implicit none
    real(kind=MATHpr), intent(in)  :: xyz(:,:)
    real(kind=MATHpr), intent(out) :: axe(3)
    real(kind=MATHpr), intent(out) :: angle
    real(kind=MATHpr)              :: cosi
    integer                        :: i

    angle = 0._MATHpr
    axe = sum( xyz, 1 ) / real( size(xyz,1), kind=MATHpr )
    axe = axe / norm2( axe )

    do i = 1,size(xyz,1)
       cosi = dot_product( axe, xyz(i,:) ) / norm2( xyz(i,:) )
       angle = max( angle, acos( min( 1._MATHpr, max( -1._MATHpr, cosi ) ) ) )
    end do

  end subroutine bounding_cone




end module mod_separation
