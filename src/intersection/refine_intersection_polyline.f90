subroutine refine_intersection_polyline( &
     curve, &
     tol2, &
     tol3 )
  use mod_math
  use mod_types_intersection
  use mod_diffgeom
  use mod_geometry
  USE MOD_UTIL
  implicit none
  LOGICAL, PARAMETER :: DEBUG = .false.
  type(type_intersection_curve), intent(inout) :: curve
  real(kind=fp),                 intent(in)    :: tol2
  real(kind=fp),                 intent(in)    :: tol3
  real(kind=fp)                                :: fractol2sqr, fractol3sqr
  type(type_intersection_polyline), pointer    :: polyline => null()
  real(kind=fp), dimension(4)                  :: lowerb, upperb
  real(kind=fp)                                :: uv(2,2), xyz(3), ctr(3), radsqr
  integer                                      :: stat
  logical                                      :: crit2, crit3
  integer                                      :: ivert, isurf
  INTEGER :: FID

  polyline => curve%polyline
  lowerb = reshape(curve%uvbox(1,1:2,1:2), [4])
  upperb = reshape(curve%uvbox(2,1:2,1:2), [4])

  fractol3sqr = 4._fp*(2._fp - tol3)*tol3**2
  fractol2sqr = 4._fp*(2._fp - tol2)*tol2**2
  
  IF ( DEBUG ) THEN
     CALL GET_FREE_UNIT(FID)
     OPEN(UNIT=FID, FILE='../debug/debug_rip_xyz0.dat', ACTION='WRITE')
     DO IVERT = 1,POLYLINE%NP
        WRITE (FID,*) POLYLINE%XYZ(1:3,IVERT)
     END DO
     CLOSE(FID)
     OPEN(UNIT=FID, FILE='../debug/debug_rip_uv0.dat', ACTION='WRITE')
     DO IVERT = 1,POLYLINE%NP
        WRITE (FID,*) POLYLINE%UV(1:2,1:2,IVERT)
     END DO
     CLOSE(FID)
  END IF
  
  
  ivert = 1
  outer_loop : do while ( ivert < polyline%np )
     IF ( DEBUG ) PRINT *,'IVERT =', IVERT, ', NP =', POLYLINE%NP
     !! get exact midpoint
     ! initial guess
     uv = 0.5_fp*(polyline%uv(1:2,1:2,ivert) + polyline%uv(1:2,1:2,ivert+1))
     ! relax to exact intersection
     call simultaneous_point_inversions( &
          curve%surf, &
          lowerb, &
          upperb, &
          stat, &
          uv, &
          xyz )
     if ( stat > 0 ) then
        PRINT *,'refine_intersection_polyline: failed to relax midpoint on segment #', ivert
        ivert = ivert + 1
        cycle outer_loop
     else
        IF ( DEBUG ) THEN
           PRINT *,'XYZ VERTICES ='
           CALL PRINT_MAT(TRANSPOSE(POLYLINE%XYZ(1:3,[IVERT,IVERT+1])))
           PRINT *,'XYZ MIDPOINT =',XYZ
           PRINT *,'UV VERTICES ='
           PRINT *,POLYLINE%UV(1:2,1:2,IVERT)
           PRINT *,POLYLINE%UV(1:2,1:2,IVERT+1)
           PRINT *,'UV MIDPOINT =',UV
        END IF
     end if
     !
     !! check chordal error
     ! in xyz space
     call circumcircle( &
       polyline%xyz(1:3,ivert), &
       xyz, &
       polyline%xyz(1:3,ivert+1), &
       ctr, &
       radsqr )
     !
     if ( radsqr < EPSfp ) radsqr = huge(1._fp)
     IF ( DEBUG ) PRINT *,'RAD3 =',SQRT(RADSQR)
     IF ( DEBUG ) PRINT *,SQRT(sum((polyline%xyz(1:3,ivert+1) - polyline%xyz(1:3,ivert))**2) / radsqr)
     crit3 = ( sum((polyline%xyz(1:3,ivert+1) - polyline%xyz(1:3,ivert))**2) < radsqr*fractol3sqr ) 
     !
     if ( crit3 ) then
        ! in uv space
        do isurf = 1,2
           call circumcircle( &
                [polyline%uv(1:2,isurf,ivert), 0._fp], &
                [uv(1:2,isurf), 0._fp], &
                [polyline%uv(1:2,isurf,ivert+1), 0._fp], &
                ctr, &
                radsqr )
           if ( radsqr < EPSfp ) radsqr = huge(1._fp)
           IF ( DEBUG ) PRINT *,'RAD2 =',SQRT(RADSQR)
           IF ( DEBUG ) PRINT *,SQRT(sum((polyline%uv(1:2,isurf,ivert+1) - polyline%uv(1:2,isurf,ivert))**2) / radsqr)
           crit2 = ( sum((polyline%uv(1:2,isurf,ivert+1) - polyline%uv(1:2,isurf,ivert))**2) < radsqr*fractol2sqr )
           if ( .not.crit2 ) exit
        end do
     else
        crit2 = .false.
     end if
     !
     IF ( DEBUG ) THEN
        PRINT *,'CRIT2, CRIT3 =', CRIT2, CRIT3
        PRINT *,''
     END IF
     if ( crit3 .and. crit2 ) then
        ivert = ivert + 1
        cycle outer_loop
     end if
     !
     !! insert midpoint
     call insert_polyline_point( &
          uv, &
          xyz, &
          stat, &
          polyline, &
          ivert )
     if ( stat == 0 ) then
        ! shift downstream polyline split points
        where ( curve%isplit(2,:) > ivert ) curve%isplit(2,:) = curve%isplit(2,:) + 1
     else
        PRINT *,'refine_intersection_polyline: failed to insert midpoint on segment #', ivert
        ivert = ivert + 1
        cycle outer_loop
     end if
     !
  end do outer_loop


  IF ( DEBUG ) THEN
     CALL GET_FREE_UNIT(FID)
     OPEN(UNIT=FID, FILE='../debug/debug_rip_xyz1.dat', ACTION='WRITE')
     DO IVERT = 1,POLYLINE%NP
        WRITE (FID,*) POLYLINE%XYZ(1:3,IVERT)
     END DO
     CLOSE(FID)
     OPEN(UNIT=FID, FILE='../debug/debug_rip_uv1.dat', ACTION='WRITE')
     DO IVERT = 1,POLYLINE%NP
        WRITE (FID,*) POLYLINE%UV(1:2,1:2,IVERT)
     END DO
     CLOSE(FID)

     PAUSE
  END IF
  
end subroutine refine_intersection_polyline
