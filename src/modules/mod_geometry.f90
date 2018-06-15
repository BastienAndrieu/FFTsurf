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
  ! ---------------------------------------------------------*

  
  
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
       axes = rotation_matrix2d( 0._fp )
       return

    elseif ( nhull == 1 ) then
       center = xy(:,hull(1))
       ranges(:) = 0._fp
       axes = rotation_matrix2d( 0._fp )
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
          if ( ranges_tmp(2) > 1 ) then
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


  
  
  ! ---------------------------------------------------------
  subroutine random_rotation_matrix3d( R )
    implicit none
    real(kind=fp), intent(out) :: R(3,3)
    
    call random_number( R(:,1) )
    R(:,1) = R(:,1) / norm2( R(:,1) )
    R(:,2) = 0._fp
    R(minloc(abs(R(:,1))),2) = 1._fp
    R(:,2) = R(:,2) - dot_product( R(:,2), R(:,1) ) * R(:,1)
    R(:,3) = cross( R(:,1), R(:,2) )

  end subroutine random_rotation_matrix3d
  ! ---------------------------------------------------------



end module mod_geometry
