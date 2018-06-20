module mod_obb
  ! Oriented Bounding Boxes (OBBs)
  use mod_math

  implicit none

  real(kind=fp), parameter :: MRGobb = real( 0.0, kind=fp )
  real(kind=fp), parameter :: EPSobb = real( 1.0e-14, kind=fp )

  type type_obb
     real(kind=fp) :: ctr(3)   ! center
     real(kind=fp) :: rng(3)   ! (half-)ranges
     real(kind=fp) :: axe(3,3) ! unit axes
  end type type_obb

contains

  ! ==================================================================
  subroutine overlap_OBBs( &
       box1, &
       box2, &
       overlap )
    ! Overlap test between two OBBs using the Separating Axis Theorem ( adapted from
    ! "Dynamic Collision Detection using Oriented Bounding Boxes", D. Eberly (2002)
    ! https://www.geometrictools.com/Documentation/DynamicCollisionDetection.pdf )
    implicit none
    type(type_obb), intent(in)  :: box1, box2
    logical,        intent(out) :: overlap
    real(kind=fp)               :: c(3,3), M(3), Rk, Rl, R
    integer                     :: i, j, k, l, n, o, p, q, s
    
    overlap = .true.

    M = box2%ctr - box1%ctr
    c = matmul( transpose(box1%axe), box2%axe )

    !PRINT *,'M =',M
    !PRINT *,'C ='
    !CALL PRINT_MAT( C )

    ! test all 15 possible separating axes
    do i = 1,15
       j = 1 + mod( i-1, 3 )

       if ( i == 4 .or. i == 7 ) c = transpose(c)

       if ( i < 7 ) then
          k = (i-j)/3 + 1
          if ( k == 1 ) then
             Rk = box1%rng(j)
             Rl = dot_product( box2%rng, abs(c(j,:)) )
             R = abs( dot_product( box1%axe(:,j), M ) )
          else
             Rk = box2%rng(j)
             Rl = dot_product( box1%rng, abs(c(j,:)) )
             R = abs( dot_product( box2%axe(:,j), M ) )
          end if
       else
          k = (i-j)/3 - 1
          l = max( 3-k, 1 )
          n = min( 5-k, 3 )
          o = max( 3-j, 1 )
          p = min( 5-j, 3 )
          q = 1 + mod(k, 3)
          s = 1 + mod(q, 3)

          Rk = box1%rng(l) * abs(c(n,j)) + box1%rng(n) * abs(c(l,j))
          Rl = box2%rng(o) * abs(c(k,p)) + box2%rng(p) * abs(c(k,o))
          R = abs( &
               c(q,j) * dot_product( box1%axe(:,s), M ) - &
               c(s,j) * dot_product( box1%axe(:,q), M ) &
               )
       end if
       
       !PRINT *,'RK =',RK
       !PRINT *,'RL =',RL
       !PRINT *,'R  =',R

       !PRINT *,'I=',I,', R-(RK+RL) =',R - RK - RL
       if ( R > Rk + Rl ) then
          !PRINT *,'R - (RK + RL) =',R - RK - RL
          !PRINT *,'DISJOINT OBBs, I =',I
          overlap = .false.
          return
       end if          
       
    end do

  end subroutine overlap_OBBs
  ! ==================================================================


  ! ==================================================================
  subroutine OBB_vertices( &
       box, &
       xyz )
    ! Returns the xyz coordinates of all eight vertices of an OBB
    implicit none
    type(type_obb), intent(in)  :: box
    real(kind=fp),  intent(out) :: xyz(8,3)
    integer                     :: i, j, k

    do i = 1,2
       do j = 1,2
          do k = 1,2
             xyz( 4*(i-1) + 2*(j-1) + k, : ) = box%ctr + &
                  matmul( real( (-1)**[i,j,k], kind=fp ) * box%rng, &
                  transpose(box%axe) )
          end do
       end do
    end do

  end subroutine OBB_vertices
  ! ==================================================================

  
  ! ==================================================================
  function is_inside_OBB( &
       xyz, &
       box, &
       tol_opt )
    ! Tests wether an xyz point lies inside the volume of an OBB
    implicit none
    real(kind=fp),  intent(in)           :: xyz(3)
    type(type_obb), intent(in)           :: box
    real(kind=fp),  intent(in), optional :: tol_opt
    logical                              :: is_inside_OBB
    real(kind=fp)                        :: tol
    integer                              :: iaxe

    if ( present(tol_opt) ) then
       tol = tol_opt
    else
       tol = 0._fp
    end if

    is_inside_OBB = .true.
    do iaxe = 1,3
       if ( abs( dot_product(xyz - box%ctr, box%axe(:,iaxe)) ) > &
            box%rng(iaxe) + tol ) then
          is_inside_OBB = .false.
          return
       end if
    end do

  end function is_inside_OBB
  ! ==================================================================


  ! ==================================================================
  subroutine chebOBB1( &
       c, &
       degr, &
       box )
    ! Computes an OBB for a curve represented by a Chebyshev series
    implicit none
    integer,        intent(in)  :: degr
    real(kind=fp),  intent(in)  :: c(degr+1,3)
    type(type_obb), intent(out) :: box
    real(kind=fp)               :: mag
    integer                     :: i

    if ( degr < 0 ) STOP 'chebOBB1 : degr < 0'
    
    ! center
    box%ctr = c(1,1:3)

    if ( degr < 1 ) then
       ! degenerate curve (single point)
       GO TO 99 ! axis-aligned bounding box
    end if

    ! axes
    mag = norm2( c(2,1:3) )
    if ( mag < EPSobb ) then
       GO TO 99 ! axis-aligned bounding box
    else
       ! oriented bounding box
       box%axe(:,1) = c(2,:) / mag

       i = minloc(abs(box%axe(:,1)),1)
       box%axe(:,2) = 0._fp
       box%axe(i,2) = 1._fp
       box%axe(:,2) = box%axe(:,2) - box%axe(i,1) * box%axe(:,1)
       box%axe(:,2) = box%axe(:,2) / norm2( box%axe(:,2) )

       box%axe(:,3) = cross( box%axe(:,1), box%axe(:,2) )

       ! ranges
       box%rng = sum( abs( matmul(c(2:,:), box%axe) ), 1 )    
       box%rng = box%rng + MRGobb

       return
    end if

    99  IF (.FALSE.) PRINT *,'chebOBB1 : axis-aligned bounding box'
    ! axes
    box%axe(:,:) = 0._fp
    do i = 1,3
       box%axe(i,i) = 1._fp
    end do

    ! ranges
    box%rng = sum( abs( c(2:,:) ), 1 )    
    box%rng = box%rng + MRGobb
    box%rng = max( EPSobb, box%rng )

  end subroutine chebOBB1
  ! ==================================================================
  

  ! ==================================================================
  subroutine chebOBB2( &
       c, &
       degr, &
       box )
    ! Computes an OBB for a surface represented by a Chebyshev Series
    implicit none
    integer,        intent(in)  :: degr(2)
    real(kind=fp),  intent(in)  :: c(degr(1)+1,degr(2)+1,3)
    type(type_obb), intent(out) :: box
    real(kind=fp)               :: vec(3,2)
    real(kind=fp)               :: mag(2)
    integer                     :: i
    
    if ( minval(degr) < 0 ) STOP 'chebOBB2 : degr < 0'

    ! center
    box%ctr = c(1,1,1:3)

    ! axes
    if ( all(degr == 0) ) then
       ! degenerate surface (single point)
       GO TO 99 ! axis-aligned bounding box

    elseif (degr(1) == 0) then
       ! degenerate surface (single point) s(u,v) = f(v)
       call chebOBB1( c(1,:,:), degr(2), box )
       return

    elseif (degr(2) == 0) then
       ! degenerate surface (single point) s(u,v) = f(u)
       call chebOBB1( c(:,1,:), degr(1), box )
       return

    else
       vec(:,1) = c(2,1,:)
       vec(:,2) = c(1,2,:)
       mag = norm2( vec, 1 )

       i = maxloc( mag, 1 )
       if ( mag(i) < EPSobb ) GO TO 99 ! axis-aligned bounding box

       box%axe(:,1) = vec(:,i) / mag(i)

       i = 1 + mod(i,2)

       box%axe(:,3) = cross( box%axe(:,1), vec(:,i) )
       mag(i) = norm2( box%axe(:,3) )

       if ( mag(i) < EPSobb ) then
          i = minloc( abs(box%axe(:,1)), 1 )
          box%axe(:,2) = 0._fp
          box%axe(i,2) = 1._fp
          box%axe(:,2) = box%axe(:,2) - &
               dot_product(box%axe(:,1),box%axe(:,2)) * box%axe(:,1)
          box%axe(:,2) = box%axe(:,2) / norm2( box%axe(:,2) )
          box%axe(:,3) = cross( box%axe(:,1), box%axe(:,2) )
       else
          box%axe(:,3) = box%axe(:,3) / mag(i)
          box%axe(:,2) = cross( box%axe(:,3), box%axe(:,1) )
       end if
    end if

    ! ranges
    box%rng = sum( abs( matmul( &
         reshape( c, [(degr(1)+1)*(degr(2)+1),3] ), &
         box%axe ) ), 1 ) - &
         abs( matmul( c(1,1,:), box%axe ) )
    box%rng = box%rng + MRGobb
    
    return

99  IF (.FALSE.) PRINT *,'chebOBB2 : axis-aligned bounding box'
    ! axes
    !box%axe(:,:) = 0._fp
    !do i = 1,3
    !   box%axe(i,i) = 1._fp
    !end do
    box%axe = identity_matrix( 3 )

    ! ranges
    box%rng = sum( abs( reshape( c,[(degr(1)+1)*(degr(2)+1),3] ) ), 1 ) - &
         abs( c(1,1,:) )
    box%rng = box%rng + MRGobb
    box%rng = max( EPSobb, box%rng )

  end subroutine chebOBB2
  ! ==================================================================




  
  ! ==================================================================
  subroutine bernOBB2( &
       b, &
       degr, &
       box )
    ! Computes an OBB for a surface represented by a Bernstein Series (Bezier surface)
    ! "Efficient bounding of displaced Bézier patches", Munkberg et al. (2010)
    implicit none
    integer,        intent(in)  :: degr(2)
    real(kind=fp),  intent(in)  :: b(degr(1)+1,degr(2)+1,3)
    type(type_obb), intent(out) :: box
    real(kind=fp)               :: vec(3,2)
    real(kind=fp)               :: mag(2)
    real(kind=fp)               :: X(size(b,1)*size(b,2),3)
    real(kind=fp), dimension(3) :: mx, mn
    integer                     :: i

    if ( minval(degr) < 0 ) STOP 'bernOBB2 : degr < 0'

    X = reshape( b, [size(b,1)*size(b,2),3] )

    ! axes
    vec(:,1) = b(degr(1)+1,1,:) - b(1,1,:) + b(degr(1)+1,degr(2)+1,:) - b(1,degr(2)+1,:)
    vec(:,2) = b(1,degr(2)+1,:) - b(1,1,:) + b(degr(1)+1,degr(2)+1,:) - b(degr(1)+1,1,:)
    mag = norm2( vec, 1 )
    i = maxloc( mag, 1 )

    if ( mag(i) < EPSobb ) then
       box%axe = identity_matrix( 3 )
    else
       box%axe(:,1) = vec(:,i) / mag(i)
       i = 1 + mod(i,2)

       box%axe(:,3) = cross( box%axe(:,1), vec(:,i) )
       mag(i) = norm2( box%axe(:,3) )
       
       if ( mag(i) < EPSobb ) then
          i = minloc( abs(box%axe(:,1)), 1 )
          box%axe(:,2) = 0._fp
          box%axe(i,2) = 1._fp
          box%axe(:,2) = box%axe(:,2) - &
               dot_product(box%axe(:,1),box%axe(:,2)) * box%axe(:,1)
          box%axe(:,2) = box%axe(:,2) / norm2( box%axe(:,2) )
          box%axe(:,3) = cross( box%axe(:,1), box%axe(:,2) )
       else
          box%axe(:,3) = box%axe(:,3) / mag(i)
          box%axe(:,2) = cross( box%axe(:,3), box%axe(:,1) )
       end if

       X = matmul( X, box%axe )
    end if

    ! ranges
    mn = minval( X, DIM=1 )
    mx = maxval( X, DIM=1 )
    box%rng = 0.5_fp * ( mx - mn )
    box%rng = box%rng + MRGobb
    box%rng = max( EPSobb, box%rng )

    ! center
    box%ctr = 0.5_fp * matmul( box%axe, mx + mn )

  end subroutine bernOBB2
    ! ==================================================================
  

  
  ! ==================================================================
  subroutine bernOBB1( &
       b, &
       degr, &
       box )
    implicit none
    integer,        intent(in)  :: degr
    real(kind=fp),  intent(in)  :: b(degr+1,3)
    type(type_obb), intent(out) :: box
    real(kind=fp)               :: vec(3), mag
    integer                     :: i
    real(kind=fp)               :: X(degr+1,3)
    real(kind=fp), dimension(3) :: mx, mn

    if ( degr < 0 ) STOP 'bernOBB1 : degr < 0'

    vec = b(degr+1,:) - b(1,:)
    mag = norm2( vec )

    if ( mag < EPSobb ) then
       if ( degr > 1 ) then
          vec = b(2,:) - b(1,:) + b(degr+1,:) - b(degr,:)
          mag = norm2( vec )
       end if
    end if

    if ( mag < EPSobb ) then
       box%axe = identity_matrix( 3 )
       X = b
    else
       box%axe(:,1) = vec / mag

       X(1:degr,:) = b(2:degr+1,:) - spread( b(1,:), dim=1, ncopies=degr )
       X(1:degr,:) = X(1:degr,:) - &
            matmul( X(1:degr,:), outer_product( box%axe(:,1), box%axe(:,1) ) )
       i = maxloc( sum( X(1:degr,:)**2, dim=2 ), 1 )
       mag = norm2( X(i,:) )
       if ( mag < EPSobb ) then
          i = minloc( abs(box%axe(:,1)), 1 )
          box%axe(:,2) = 0._fp
          box%axe(i,2) = 1._fp
          box%axe(:,2) = box%axe(:,2) - &
               dot_product(box%axe(:,1),box%axe(:,2)) * box%axe(:,1)
          box%axe(:,2) = box%axe(:,2) / norm2( box%axe(:,2) )
       else
          box%axe(:,2) = X(i,:) / mag
       end if
       box%axe(:,3) = cross( box%axe(:,1), box%axe(:,2) )
       
       X = matmul( b, box%axe )
    end if

    ! ranges
    mn = minval( X, DIM=1 )
    mx = maxval( X, DIM=1 )
    box%rng = 0.5_fp * ( mx - mn )
    box%rng = box%rng + MRGobb
    !box%rng = max( EPSobb, box%rng )

    ! center
    box%ctr = 0.5_fp * matmul( box%axe, mx + mn )

  end subroutine bernOBB1
  ! ==================================================================
  


subroutine print_obb( box )
  implicit none
  type(type_obb), intent(in) :: box
  integer                    :: i
  print *,'ctr =',box%ctr
  print *,'rng =',box%rng
  print *,'axe ='
  do i = 1,3
     print *,box%axe(:,i)
  end do
  
end subroutine print_obb


subroutine write_obb( box, filename )
  use mod_util
  implicit none
  type(type_obb), intent(in) :: box
  character(*),   intent(in) :: filename
  integer                    :: i, fileunit

  call get_free_unit( fileunit )
  if ( fileunit == 0 ) STOP "write_obb : could not find free unit"
  open( unit = fileunit, file = filename, action = "write" )
  write (fileunit,*) box%ctr
  write (fileunit,*) box%rng
  do i = 1,3
     write (fileunit,*) box%axe(:,i)
  end do
  close( fileunit )

end subroutine write_obb


end module mod_obb
