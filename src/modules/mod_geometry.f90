module mod_geometry

  use mod_math

contains

  ! ---------------------------------------------------------
  function angle_between_vectors_2d( v1, v2 ) result( a )
    ! Angle orienté (radians) entre deux vecteurs du plan
    implicit none
    real(kind=fp), dimension(2), intent(in) :: v1, v2
    real(kind=fp)                           :: a

    a = atan2( v1(2)*v2(1) - v1(1)*v2(2), v1(1)*v2(1) + v1(2)*v2(2) )

  end function angle_between_vectors_2d
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  function angle_between_vectors_3d( v1, v2 ) result( a )
    ! Angle orienté (radians) entre deux vecteurs du plan
    implicit none
    real(kind=fp), dimension(3), intent(in) :: v1, v2
    real(kind=fp)                           :: a

    a = atan2(norm2(cross(v1,v2)), dot_product(v1,v2))

  end function angle_between_vectors_3d
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  function rotation_matrix2d( angle ) result( R )
    ! Matrice de rotation en 2d
    implicit none
    real(kind=fp), intent(in) :: angle
    real(kind=fp)             :: R(2,2)
    real(kind=fp)             :: c, s

    c = cos( angle )
    s = sin( angle )
    R(:,1) = [ c , -s ]
    R(:,2) = [ s, c ]

  end function rotation_matrix2d
  ! ---------------------------------------------------------



  function triangle_normal( xyz, unitize ) result( nor )
    implicit none
    real(kind=fp),     intent(in) :: xyz(3,3)
    logical, optional, intent(in) :: unitize
    real(kind=fp)                 :: nor(3)

    nor = cross(xyz(:,2) - xyz(:,1), xyz(:,3) - xyz(:,1))
    if ( present(unitize) ) then
       if ( unitize ) nor = nor / norm2(nor)
    end if
    
  end function triangle_normal

  

  ! ---------------------------------------------------------
  subroutine convex_hull_2d( &
       xy, &
       n, &
       hull, &
       nhull )
    ! Enveloppe convexe d'un ensemble de points du plan.
    implicit none
    integer,       intent(in)  :: n       ! nombre de points
    real(kind=fp), intent(in)  :: xy(2,n) ! coordonnées xy des points
    integer,       intent(out) :: hull(n) ! liste des points sur l'enveloppe convexe
    integer,       intent(out) :: nhull   ! nombre de points sur l'enveloppe convexe
    real(kind=fp)              :: dxy(2), dxy_prev(2), angle, anglenext, dist, distnext
    integer                    :: i, first, next

    !PRINT *,'';PRINT *,'';PRINT *,'';

    ! Cas particuliers n = 0, 1 ou 2
    if ( n < 1 ) then
       ! aucun point => l'enveloppe est vide
       nhull = 0
       return
    elseif ( n == 1 ) then
       ! un seul point => l'enveloppe est ce point
       nhull = 1
       hull(1) = 1
       return
    elseif ( n == 2 ) then
       ! deux points ...
       if ( norm2( xy(:,2) - xy(:,1) ) < MATHeps ) then
          ! ... coincidents => un seul point
          nhull = 1
          hull(1) = 1
       else
          ! ... disctincts => l'enveloppe est le segment entre ces points
          nhull = 2
          hull(1:2) = [1,2]
       end if
       return
    end if


    ! on commence par le point d'abscisse minimale
    ! ( en cas d'égalité, on garde le point d'ordonnée minimale )
    first = 1
    do i = 2,n
       if ( xy(1,i) < xy(1,first) .or. &
            abs(xy(1,i) - xy(1,first)) < MATHeps .and.  xy(2,i) < xy(2,first)  ) then
          first = i
       end if
    end do

    ! on ajoute ce premier point à l'enveloppe convexe
    nhull = 1
    hull(1) = first   

    ! Construction incrémentale de l'enveloppe convexe :
    ! connaissant $p_{i-1}$ et $p_{i}$, on cherche le point 
    ! $p_{i+1}$ qui maximise $\angle p_{i-1} p_{i} p_{i+1} $.
    dxy_prev = real( [0,-1],kind=fp ) ! point p_{0} fictif à p_{1} - (0,1)
    do
       !PRINT *,'HULL =',HULL(1:NHULL)
       !PRINT *,'DXY PREV =',DXY_PREV
       anglenext = MATHpi
       distnext = huge(1._fp)
       next = 0
       do i = 1,n
          if ( i == hull(nhull) ) cycle
          dxy = xy(:,i) - xy(:,hull(nhull))
          dist = sum( dxy**2 )
          if ( dist < MATHeps**2 ) cycle ! points coincidents

          angle = angle_between_vectors_2d( dxy_prev, dxy )
          !PRINT *,'I =',I
          !PRINT *,'    DIST =',DIST
          !PRINT *,'     DXY =',DXY
          !PRINT *,'   ANGLE =',ANGLE
          if ( next == 0 .or. angle > anglenext ) then
             anglenext = angle
             distnext = dist
             next = i
          elseif ( next /= 0 .and. angle > anglenext - MATHeps ) then
             ! en cas d'égalité, on retient le point le plus proche
             !PRINT *,'TIE, CHOSE CLOSEST POINT'
             if ( dist < distnext ) then
                anglenext = angle
                distnext = dist
                next = i
             end if
          end if
       end do
       ! si le point suivant est identique au premier, la construction est terminée
       if ( next == first ) return

       ! on ajoute le point suivant à l'enveloppe convexe
       nhull = nhull + 1
       IF ( NHULL > N ) THEN
          !PRINT *,'XY ='
          !CALL PRINT_MAT( TRANSPOSE(XY) )
          !PRINT *,'HULL='
          !PRINT *,HULL
          !PRINT *,'NEXT =', NEXT
          nhull = nhull - 1
          RETURN
       END IF
       hull(nhull) = next
       dxy_prev = xy(:,hull(nhull-1)) - xy(:,next)

    end do

  end subroutine convex_hull_2d
  ! ---------------------------------------------------------




  ! ---------------------------------------------------------
  subroutine minimum_area_OBB_2d( &
       xy, &
       n, &
       center, &
       ranges, &
       axes )
    ! Retourne la boîte orientée d'aire minimale qui englobe un ensemble de points du plan.
    ! La boîte est definie par son centre, ses axes, et ses demi-côtés (ranges)
    implicit none
    integer,       intent(in)  :: n         ! nombre de points
    real(kind=fp), intent(in)  :: xy(2,n)   ! coordonnées xy des points
    real(kind=fp), intent(out) :: center(2) ! coordonnées xy du centre de la boîte
    real(kind=fp), intent(out) :: ranges(2) ! demi-côtés de la boîte
    real(kind=fp), intent(out) :: axes(2,2) ! axes de la boîte
    integer                    :: hull(n), nhull
    real(kind=fp)              :: area_tmp, area
    real(kind=fp)              :: vec(2), rot(2,2)
    real(kind=fp), allocatable :: xy_rot(:,:)
    real(kind=fp)              :: mn(2), mx(2), ranges_tmp(2)
    integer                    :: i

    ! enveloppe convexe de l'ensemble des points
    call convex_hull_2d( xy, n, hull, nhull )

    ! cas particuliers nhull = 0, 1 ou 2
    if ( nhull < 1 ) then
       center(:) = 0._fp
       ranges(:) = 0._fp
       axes = identity_matrix(2) ! rotation_matrix2d( 0._fp )
       return

    elseif ( nhull == 1 ) then
       center = xy(:,hull(1))
       ranges(:) = 0._fp
       axes = identity_matrix(2) ! rotation_matrix2d( 0._fp )
       return

    elseif( nhull == 2 ) then
       center = 0.5_fp * sum( xy(:,hull(1:2)), 2 )
       vec = xy(:,hull(2)) - xy(:,hull(1))
       ranges(1) = 0.5_fp * norm2( vec )
       ranges(2) = 0._fp
       axes = rotation_matrix2d( -atan2( vec(2), vec(1) ) )
       return

    end if

    ! La boîte a un côté porté par une arête de l'enveloppe convexe
    allocate( xy_rot(2,nhull) )

    area = huge(1.d0)
    do i = 1,nhull
       ! i-ème arête de l'enveloppe convexe
       vec = xy(:,hull( 1+mod(i,nhull) )) - xy(:,hull(i))

       ! rotation pour aligner l'arête avec l'axe x
       rot = rotation_matrix2d( atan2( vec(2), vec(1) ) )
       xy_rot = matmul( rot, xy(:,hull(1:nhull)) )

       ! étendue de l'enveloppe convexe dans le repère après rotation
       mn = minval( xy_rot, 2 )
       mx = maxval( xy_rot, 2 )
       ranges_tmp = mx - mn
       area_tmp = ranges_tmp(1) * ranges_tmp(2)

       if ( area_tmp < area ) then
          area = area_tmp
          ! rotation vers le repère initial
          rot = transpose( rot )
          ! coordonnées du centre de la boîte dans le repère initial
          center = matmul( rot, 0.5_fp * (mn + mx) )
          ! on ordonne les axes pour avoir ranges(1) >= ranges(2)
          if ( ranges_tmp(2) > ranges_tmp(1) ) then
             ranges = 0.5_fp * ranges_tmp([2,1])
             axes(:,1) = rot(:,2)
             axes(:,2) = -rot(:,1)
          else
             ranges = 0.5_fp * ranges_tmp
             axes = rot
          end if
       end if

    end do

    deallocate( xy_rot )

  end subroutine minimum_area_OBB_2d
  ! ---------------------------------------------------------





  subroutine slerp( &
       p0, &
       p1, &
       t, &
       n, &
       p )
    ! Spherical linear interpolation between vectors p0 and p1
    ! ( see https://en.wikipedia.org/wiki/Slerp )
    ! p0 and p1 need not be unit vectors
    ! t should be a subset of [0,1]
    real(kind=fp), intent(in)  :: p0(3), p1(3)
    integer,       intent(in)  :: n
    real(kind=fp), intent(in)  :: t(n)
    real(kind=fp), intent(out) :: p(3,n)
    real(kind=fp)              :: r0, r1, a, s

    r0 = norm2(p0) 
    r1 = norm2(p1)

    a = acos( min( 1._fp, max( -1._fp, dot_product(p0, p1) / (r0 * r1) ) ) )

    s = sin(a)

    !p = matmul( p0, sin( (1._fp - t)*a ) ) + &
    !     matmul( p1, sin( t*a ) )
    p = outer_product( p0, sin( (1._fp - t)*a ) ) + &
         outer_product( p1, sin( t*a ) )
    if ( abs(s) > MATHeps ) p = p / s

  end subroutine slerp




  subroutine complete_orthonormal_basis( u, v, w )
    ! given a unit 3-vector u, returns 3-vectors v and w such
    ! that (u,v,w) is a direct orthonormal basis 
    implicit none
    real(kind=fp), intent(in)  :: u(3)
    real(kind=fp), intent(out) :: v(3), w(3)
    
    v(:) = 0._fp
    v(minloc(u,1)) = 1._fp
    v = v - dot_product( u, v ) * u
    v = v / norm2(v)

    w = cross(u, v)    

  end subroutine complete_orthonormal_basis



  ! ---------------------------------------------------------
  subroutine random_rotation_matrix3d( R )
    implicit none
    real(kind=fp), intent(out) :: R(3,3)

    call random_number( R(:,1) )
    R(:,1) = R(:,1) / norm2( R(:,1) )
    !R(:,2) = 0._fp
    !R(minloc(abs(R(:,1))),2) = 1._fp
    !R(:,2) = R(:,2) - dot_product( R(:,2), R(:,1) ) * R(:,1)
    !R(:,2) = R(:,2) / norm2( R(:,2) )
    !R(:,3) = cross( R(:,1), R(:,2) )
    call complete_orthonormal_basis( R(:,1), R(:,2), R(:,3) )

  end subroutine random_rotation_matrix3d
  ! ---------------------------------------------------------
















  subroutine quickhull( &
       xy, &
       nxy, &
       hull, &
       nhull )
    use mod_tolerances
    implicit none
    integer,       intent(in)  :: nxy
    real(kind=fp), intent(in)  :: xy(2,nxy)
    integer,       intent(out) :: hull(nxy)
    integer,       intent(out) :: nhull
    integer                    :: a, b, i
    integer                    :: part(nxy), npart

    IF ( .false. ) THEN
       PRINT *,'XY ='
       CALL PRINT_MAT( TRANSPOSE(XY) )
    END IF

    if ( nxy < 2 ) then
       nhull = nxy
       hull = 1
       return
    end if

    ! find leftmost (a) and rightmost (b) points
    a = 1
    b = 1
    do i = 2,nxy
       if ( xy(1,i) < xy(1,a) .or. &
            xy(1,i) < xy(1,a) + EPSfp .and. xy(2,i) < xy(2,a) ) then
          a = i
       elseif ( xy(1,i) > xy(1,b) .or. &
            xy(1,i) > xy(1,b) - epsilon(1._fp) .and. xy(2,i) > xy(2,b) ) then
          b = i
       end if
       !if ( xy(1,i) < xy(1,a) - EPSfp .or. &
       !     xy(1,i) < xy(1,a) + EPSfp .and. xy(2,i) < xy(2,a) - EPSfp ) then
       !   a = i
       !elseif ( xy(1,i) > xy(1,b) + EPSfp .or. &
       !     xy(1,i) > xy(1,b) - EPSfp .and. xy(2,i) > xy(2,b) + EPSfp ) then
       !   b = i
       !end if
    end do

    !PRINT *,'XYA =', XY(:,A)
    !PRINT *,'XYB =', XY(:,B)

    ! partition into two 
    nhull = 1
    hull(nhull) = a

    call points_right_to_line( &
         xy(:,a), &
         xy(:,b), &
         xy, &
         nxy, &
         part, &
         npart )
    call quickhull_partition( &
         xy, &
         nxy, &
         part, &
         npart, &
         a, &
         b, &
         hull, &
         nhull )

    nhull = nhull + 1
    hull(nhull) = b

    call points_right_to_line( &
         xy(:,b), &
         xy(:,a), &
         xy, &
         nxy, &
         part, &
         npart )
    call quickhull_partition( &
         xy, &
         nxy, &
         part, &
         npart, &
         b, &
         a, &
         hull, &
         nhull )

  end subroutine quickhull




  recursive subroutine quickhull_partition( &
       xy, &
       nxy, &
       part, &
       npart, &
       a, &
       b, &
       hull, &
       nhull )
    use mod_tolerances
    implicit none
    integer,       intent(in)    :: nxy
    real(kind=fp), intent(in)    :: xy(2,nxy)
    integer,       intent(in)    :: part(nxy)
    integer,       intent(in)    :: npart
    integer,       intent(in)    :: a, b
    integer,       intent(inout) :: hull(nxy)
    integer,       intent(inout) :: nhull
    integer                      :: newpart(nxy)
    integer                      :: nnewpart
    real(kind=fp)                :: abn(2), distp, distc
    integer                      :: c, p, i

    if ( npart < 1 ) return    

    ! find the point c with maximal distance from the line [ab]
    abn(1) = xy(2,a) - xy(2,b)
    abn(2) = xy(1,b) - xy(1,a) 
    c = 0
    distc = 0._fp
    do i = 1,npart
       p = part(i)
       if ( p == a .or. p == b ) cycle
       distp = &
            ( xy(1,a) - xy(1,p) ) * abn(1) + &
            ( xy(2,a) - xy(2,p) ) * abn(2)
       if ( distp > distc ) then!+ EPSfp ) then
          c = p
          distc = distp
       end if
    end do

    ! recursion on new partitions:
    ! 1) points right to [ac]
    call points_right_to_line( &
         xy(:,a), &
         xy(:,c), &
         xy, &
         nxy, &
         newpart, &
         nnewpart )
    call quickhull_partition( &
         xy, &
         nxy, &
         newpart, &
         nnewpart, &
         a, &
         c, &
         hull, &
         nhull )

    nhull = nhull + 1
    hull(nhull) = c

    ! 2) points right to [cb]
    call points_right_to_line( &
         xy(:,c), &
         xy(:,b), &
         xy, &
         nxy, &
         newpart, &
         nnewpart )
    call quickhull_partition( &
         xy, &
         nxy, &
         newpart, &
         nnewpart, &
         c, &
         b, &
         hull, &
         nhull )

  end subroutine quickhull_partition






  subroutine points_right_to_line( &
       xya, &
       xyb, &
       xy, &
       nxy, &
       list, &
       nlist )
    use mod_tolerances
    implicit none
    real(kind=fp), intent(in)  :: xya(2), xyb(2)
    integer,       intent(in)  :: nxy
    real(kind=fp), intent(in)  :: xy(2,nxy)
    integer,       intent(out) :: list(nxy)
    integer,       intent(out) :: nlist
    real(kind=fp)              :: abn(2), dist
    integer                    :: i

    abn(1) = xya(2) - xyb(2)
    abn(2) = xyb(1) - xya(1) 

    nlist = 0
    do i = 1,nxy
       dist = &
            ( xy(1,i) - xya(1) ) * abn(1) + &
            ( xy(2,i) - xya(2) ) * abn(2)
       if ( dist < -EPSfp ) then
          nlist = nlist + 1
          list(nlist) = i
       end if
    end do

  end subroutine points_right_to_line






  function distance_from_line( &
       p, &
       a, &
       b )
    implicit none
    real(kind=fp), dimension(2), intent(in) :: p, a, b
    real(kind=fp)                           :: distance_from_line
    real(kind=fp)                           :: lab
    
    distance_from_line = (a(2) - b(2))*(p(1) - a(1)) + (b(1) - a(1))*(p(2) - a(2))
    lab = sum((b - a)**2)
    if ( lab > EPSfp ) distance_from_line = distance_from_line / sqrt(lab)

  end function distance_from_line




  subroutine point_in_polygon( &
       px, &
       py, &
       np, &
       x, &
       y, &
       inside )
    implicit none
    integer,                      intent(in)  :: np
    real(kind=fp), dimension(np), intent(in)  :: px, py
    real(kind=fp),                intent(in)  :: x, y
    logical,                      intent(out) :: inside
    integer                                   :: i, j

    inside = .false.
    do i = 1,np
       j = 1 + mod(i,np)
       if ( (py(i) > y) .neqv. (py(j) > y) .and. &
            (py(i) > py(j)) .neqv. ((x - px(i))*(py(i) - py(j)) > (py(i) - y)*(px(j) - px(i))) ) then
          inside = .not.inside
       end if
    end do
    
  end subroutine point_in_polygon




  subroutine circumcircle( &
       p1, &
       p2, &
       p3, &
       ctr, &
       radsqr )
    implicit none
    real(kind=fp), intent(in), dimension(2) :: p1, p2, p3
    real(kind=fp), intent(out)              :: ctr(2)
    real(kind=fp), intent(out)              :: radsqr
    real(kind=fp), dimension(2)             :: u, v, w
    real(kind=fp)                           :: usqr, vsqr, denom
    
    u = p2 - p1
    v = p3 - p1
    usqr = sum(u**2)
    vsqr = sum(v**2)
    w = usqr*v - vsqr*u
    
    denom = usqr*vsqr - dot_product(u,v)**2
    if ( denom < EPSfp ) then
       ctr = p1
       radsqr = 0._fp
    else
       ctr = 0.5_fp * (dot_product(v,w)*u - dot_product(u,w)*w) / denom
       radsqr = sum(ctr**2)
    end if
    
  end subroutine circumcircle



  
  ! ==================================================================
  ! cf. "Computational Geometry, Algorithms, and Applications", M. de Berg et al. (2000)
  subroutine smallest_enclosing_ball( &
       n, &
       dim, &
       p, &
       ctr, &
       radsqr )
    implicit none
    integer,       intent(in)  :: n
    integer,       intent(in)  :: dim
    real(kind=fp), intent(in)  :: p(dim,n)
    real(kind=fp), intent(out) :: ctr(dim)
    real(kind=fp), intent(out) :: radsqr
    integer                    :: i

    if ( n < 2 ) then
       ctr = p(1:dim,1)
       radsqr = 0._fp
       return
    end if

    ctr = 0.5_fp*sum(p(1:dim,1:2), 2)
    radsqr = 0.25_fp*sum( (p(1:dim,1) - p(1:dim,2))**2 )

    do i = 3,n
       if ( sum( (p(1:dim,i) - ctr)**2 ) > radsqr ) then
          call smallest_enclosing_ball_with_point( &
               i-1, &
               dim, &
               p(1:dim,1:i-1), &
               p(1:dim,i), &
               ctr, &
               radsqr )
       end if
    end do

  end subroutine smallest_enclosing_ball


  subroutine smallest_enclosing_ball_with_point( &
       n, &
       dim, &
       p, &
       q, &
       ctr, &
       radsqr )
    implicit none
    integer,       intent(in)  :: n
    integer,       intent(in)  :: dim
    real(kind=fp), intent(in)  :: p(dim,n)
    real(kind=fp), intent(in)  :: q(dim)
    real(kind=fp), intent(out) :: ctr(dim)
    real(kind=fp), intent(out) :: radsqr
    integer                    :: i

    ctr = 0.5_fp*(p(1:dim,1) + q)
    radsqr = 0.25_fp*sum( (p(1:dim,1) - q)**2 )

    do i = 2,n
       if ( sum( (p(1:dim,i) - ctr)**2 ) > radsqr ) then
          call smallest_enclosing_ball_with_two_points( &
               n, &
               dim, &
               p(1:dim,1:i-1), &
               p(1:dim,i), &
               q, &
               ctr, &
               radsqr )
       end if
    end do

  end subroutine smallest_enclosing_ball_with_point


  subroutine smallest_enclosing_ball_with_two_points( &
       n, &
       dim, &
       p, &
       q1, &
       q2, &
       ctr, &
       radsqr )
    implicit none
    integer,       intent(in)  :: n
    integer,       intent(in)  :: dim
    real(kind=fp), intent(in)  :: p(dim,n)
    real(kind=fp), intent(in)  :: q1(dim)
    real(kind=fp), intent(in)  :: q2(dim)
    real(kind=fp), intent(out) :: ctr(dim)
    real(kind=fp), intent(out) :: radsqr
    integer                    :: i

    ctr = 0.5_fp*(q1 + q2)
    radsqr = 0.25_fp*sum( (q1 - q2)**2 )

    do i = 1,n
       if ( sum( (p(1:dim,i) - ctr)**2 ) > radsqr ) then
          call circumsphere( &
               dim, &
               p(1:dim,i), &
               q1, &
               q2, &
               ctr, &
               radsqr )
       end if
    end do

  end subroutine smallest_enclosing_ball_with_two_points
  ! ==================================================================





  subroutine circumsphere( &
       dim, &
       p1, &
       p2, &
       p3, &
       ctr, &
       radsqr )
    implicit none
    integer,       intent(in)                 :: dim
    real(kind=fp), intent(in), dimension(dim) :: p1, p2, p3
    real(kind=fp), intent(out)                :: ctr(dim)
    real(kind=fp), intent(out)                :: radsqr
    real(kind=fp), dimension(dim)             :: u, v, w
    real(kind=fp)                             :: usqr, vsqr, denom

    u = p2 - p1
    v = p3 - p1
    usqr = sum(u**2)
    vsqr = sum(v**2)
    w = usqr*v - vsqr*u

    denom = usqr*vsqr - dot_product(u,v)**2
    if ( denom < EPSfp ) then
       ctr = (p1 + p2 + p3)/3._fp
       radsqr = 0._fp
    else
       ctr = 0.5_fp * (dot_product(v,w)*u - dot_product(u,w)*w) / denom
       radsqr = sum(ctr**2)
       ctr = ctr + p3
    end if

  end subroutine circumsphere

end module mod_geometry
