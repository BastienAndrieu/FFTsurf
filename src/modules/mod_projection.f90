module mod_projection

  use mod_math

  implicit none

contains


  subroutine projection_hyperedge( &
       brep, &
       hyperedge, &
       iedge, &
       uv, &
       xyz, &
       duv, &
       dxyz, &
       iedgetmp, &
       uvtmp, &
       debug )
    use mod_linalg
    use mod_types_brep
    use mod_brep
    use mod_halfedge
    use mod_types_intersection
    use mod_intersection
    use mod_tolerances
    use mod_hypergraph
    implicit none
    logical, intent(in) :: debug
    type(type_brep),      intent(in), target  :: brep
    type(type_hyperedge), intent(in)          :: hyperedge
    integer,              intent(in)          :: iedge
    real(kind=fp),        intent(in)          :: uv(2,2)
    real(kind=fp),        intent(in)          :: xyz(3)
    real(kind=fp),        intent(in)          :: duv(2,2)
    real(kind=fp),        intent(in)          :: dxyz(3)
    integer,              intent(out)         :: iedgetmp
    real(kind=fp),        intent(out)         :: uvtmp(2,2)
    real(kind=fp)                             :: duvtmp(2,2)
    real(kind=fp), dimension(3)               :: xyzprev, xyztmp, xyztarget, dxyztmp, p
    type(type_intersection_curve), pointer    :: curve => null()
    type(type_intersection_polyline), pointer :: polyline => null()
    real(kind=fp)                             :: w, wfirst, wlast, ds, duv_ds(2,2,2), dxyz_ds(3,2)
    integer                                   :: jedge, ihedg(2), ifirst, ilast, sens, np, k, stat
    real(kind=fp), dimension(4)               :: lowerb, upperb

    IF ( DEBUG ) PRINT *,'PROJECTION HYPEREDGE'
    upperb(:) = 2._fp
    lowerb = -upperb
    
    do jedge = 1,hyperedge%ne
       if ( hyperedge%halfedges(1,jedge) == iedge ) then
          ihedg = hyperedge%halfedges(:,jedge)
          iedgetmp = ihedg(1)
          curve => brep%edges(iedgetmp)%curve
          polyline => curve%polyline
          call get_polyline_endpoints( &
               brep, &
               ihedg, &
               ifirst, &
               ilast, &
               sens, &
               np )
          exit
       end if
    end do

    uvtmp = uv
    duvtmp = duv
    xyztmp = xyz
    dxyztmp = dxyz

    IF ( DEBUG ) THEN
       PRINT *,'IHEDGE =', IHEDG
       PRINT *,' FACES =', get_face(brep, [iedgetmp,1]), get_face(brep, [iedgetmp,2])
       PRINT *,'    UV =', UV
       PRINT *,'   XYZ =', XYZ
       PRINT *,'  DXYZ =', DXYZ
    END IF
    
    do k = 1,hyperedge%ne
       xyztarget = xyztmp + dxyztmp
       
       p = real((-1)**sens, kind=fp) * curve%param_vector
       wfirst = dot_product(p, polyline%xyz(:,ifirst))
       wlast  = dot_product(p, polyline%xyz(:,ilast))
       w = dot_product(p, xyztarget)
       
       if ( is_in_closed_interval(w, wfirst, wlast, tolerance=EPSuv) ) then
          ! relax temporary point onto true intersection underlying the hyperedge
          uvtmp = uvtmp + duvtmp ! initial iterate
          xyzprev = xyztmp
          call newton_intersection_polyline( &
               curve%surf, &
               lowerb, &
               upperb, &
               xyzprev, &
               sum(dxyztmp**2), &
               stat, &
               uvtmp, &
               xyztmp )
          if ( stat > 0 ) then
             PRINT *,'projection_hyperedge : failed to reproject onto transversal intersection curve, k =',k
             PAUSE
          else
             IF ( DEBUG ) THEN
                PRINT *,'CONVERGED :'
                PRINT *,'IHEDGE =', IHEDG
                PRINT *,' FACES =', get_face(brep, [iedgetmp,1]), get_face(brep, [iedgetmp,2])
                PRINT *,'    UV =', UVTMP
                PRINT *,'   XYZ =', XYZTMP
             END IF
             return
          end if
       else
          ! the target point is on another halfedge
          if ( w < wfirst ) then
             ! move on to previous halfedge on the hyperedge
             jedge = 1 + mod(jedge + hyperedge%ne - 2,hyperedge%ne) ! previous
          else
             ! move on to next halfedge on the hyperedge
             jedge = 1 + mod(jedge,hyperedge%ne) ! next
          end if

          ihedg = hyperedge%halfedges(:,jedge)
          iedgetmp = ihedg(1)
          curve => brep%edges(iedgetmp)%curve
          polyline => curve%polyline
          call get_polyline_endpoints( &
               brep, &
               ihedg, &
               ifirst, &
               ilast, &
               sens, &
               np )

          if ( w < wfirst ) then
             uvtmp = polyline%uv(:,:,ilast)
             xyztmp = polyline%xyz(:,ilast)
          else
             uvtmp = polyline%uv(:,:,ifirst)
             xyztmp = polyline%xyz(:,ifirst)
          end if

          dxyztmp = xyztarget - xyztmp
          call diffgeom_intersection( &
               curve%surf, &
               uvtmp, &
               duv_ds, &
               dxyz_ds, &
               stat )
          ds = dot_product(dxyztmp, dxyz_ds(:,1))
          duvtmp = ds * duv_ds(:,1,:)
          dxyztmp = ds * dxyz_ds(:,1)

       end if
    end do
    STOP 'projection_hyperedge : failed to reproject'
    
  end subroutine projection_hyperedge









  

  subroutine projection_hyperface( &
       brep, &
       iface, &
       uv, &
       xyz, &
       duv, &
       dxyz, &
       ifacetmp, &
       uvtmp, &
       debug )
    use mod_linalg
    use mod_diffgeom
    use mod_types_brep
    use mod_halfedge
    use mod_types_intersection
    use mod_intersection
    use mod_tolerances
    use mod_geometry
    use mod_polynomial
    implicit none
    logical, intent(in) :: debug
    type(type_brep), intent(in)               :: brep
    integer,         intent(in)               :: iface
    real(kind=fp),   intent(in)               :: uv(2)
    real(kind=fp),   intent(in)               :: xyz(3)
    real(kind=fp),   intent(in)               :: duv(2)
    real(kind=fp),   intent(in)               :: dxyz(3)
    integer,         intent(out)              :: ifacetmp
    real(kind=fp),   intent(out)              :: uvtmp(2)
    real(kind=fp)                             :: duvtmp(2)
    real(kind=fp)                             :: xyztmp(3)
    real(kind=fp)                             :: dxyztmp(3)
    type(type_intersection_polyline), pointer :: polyline => null()
    type(type_matrix)                         :: uvpoly(2)
    integer                                   :: i1, i2, np, sens
    integer                                   :: ninter
    integer, allocatable                      :: isegments(:,:)
    real(kind=fp), allocatable                :: lambda(:,:)
    real(kind=fp)                             :: det, w, jac(3,2)
    logical                                   :: changeface, singular
    integer                                   :: iwire, ihedg(2), istart, it, ivar, iinter
    INTEGER :: I
    real(kind=fp)                             :: uvtng(2,2), duv_ds(2,2,2), dxyz_ds(3,2)
    real(kind=fp), dimension(4)               :: lowerb, upperb
    integer                                   :: stat

    upperb(:) = 1._fp + EPSuv
    lowerb = -upperb
    
    IF ( DEBUG ) PRINT *,'PROJECTION HYPERFACE'
    
    xyztmp = xyz
    uvtmp = uv
    ifacetmp = iface
    duvtmp = duv
    dxyztmp = dxyz

    IF ( DEBUG ) THEN
       PRINT *,'DXYZTMP  =',dxyztmp
       PRINT *,'DUVTMP   =',duvtmp
    END IF

    allocate(uvpoly(1)%mat(2,2))

    it = 0
    outer_loop : do ! <------------------------------------------------------------------------------+
       IF ( DEBUG ) PRINT *,'FA',ifacetmp
       IF ( DEBUG ) PAUSE
       changeface = .false.                                                                          !
       uvpoly(1)%mat(:,1) = uvtmp                                                                    !
       uvpoly(1)%mat(:,2) = uvtmp + duvtmp                                                           !
       !                                                                                             !
       brep_wires : do iwire = 0,brep%faces(ifacetmp)%ninner ! <----------------------------------+  !
          ! get index of first halfedge on the wire                                               !  !
          if ( iwire == 0 ) then ! <---------------------+                                        !  !
             ihedg = brep%faces(ifacetmp)%outer          !                                        !  !
          else ! ----------------------------------------+                                        !  !
             ihedg = brep%faces(ifacetmp)%inner(:,iwire) !                                        !  !
          end if ! <-------------------------------------+                                        !  !
          ! start brep vertex                                                                     !  !
          istart = get_orig(brep, ihedg)                                                          !  !
          ! cycle the wire's halfedges                                                            !  !
          brep_halfedges : do ! <--------------------------------------------------------------+  !  !
             IF ( DEBUG ) PRINT *,'HE',ihedg
             ! get current halfedge's uv polyline                                              !  !  !
             i1 = brep%edges(ihedg(1))%curve%isplit(2,brep%edges(ihedg(1))%isplit)             !  !  !
             i2 = brep%edges(ihedg(1))%curve%isplit(2,brep%edges(ihedg(1))%isplit + 1)         !  !  !
             np = i2 - i1 + 1                                                                  !  !  !
             sens = 1 + mod(ihedg(2),2)                                                        !  !  !
             polyline => brep%edges(ihedg(1))%curve%polyline                                   !  !  !
             allocate(uvpoly(2)%mat(2,np))                                                     !  !  !
             uvpoly(2)%mat(1:2,1:np) = polyline%uv(1:2,sens,i1:i2)                             !  !  !
             !                                                                                 !  !  !
             ! compute intersection between the polyline and the uv displacement               !  !  !
             ninter = 0                                                                        !  !  !
             call intersect_2Dpolylines( &                                                     !  !  !
                  uvpoly, &                                                                    !  !  !
                  [1,1], &                                                                     !  !  !
                  [2,np], &                                                                    !  !  !
                  ninter, &                                                                    !  !  !
                  isegments, &                                                                 !  !  !
                  lambda )                                                                     !  !  !
             IF ( DEBUG ) THEN
                OPEN(UNIT=13, FILE='polylines/xy1.dat', ACTION='WRITE')
                WRITE (13,*) 2
                DO I = 1,2
                   WRITE (13,*) UVPOLY(1)%MAT(:,I)
                END DO
                CLOSE(13)

                OPEN(UNIT=13, FILE='polylines/xy2.dat', ACTION='WRITE')
                WRITE (13,*) NP
                DO I = 1,NP
                   WRITE (13,*) UVPOLY(2)%MAT(:,I)
                END DO
                CLOSE(13)

                OPEN(UNIT=13, FILE='polylines/result.dat', ACTION='WRITE')
                WRITE (13,*) NINTER
                DO I = 1,NINTER
                   WRITE (13,*) ISEGMENTS(:,I)
                END DO
                DO I = 1,NINTER
                   WRITE (13,*) LAMBDA(:,I)
                END DO
                CLOSE(13)

                PRINT *,'--> CHECK POLYINES'
                PAUSE
             END IF
             !                                                                                 !  !  !
             if ( ninter > 0 ) then ! <-----------------------------------------------------+  !  !  !
                ! pick closest point to target uv point                                     !  !  !  !
                iinter = maxloc(lambda(1,1:ninter), 1)                                      !  !  !  !
                i1 = i1 + isegments(2,iinter) - 1                                           !  !  !  !
                i2 = i1 + 1                                                                 !  !  !  !
                ! on which side of the polyline lies the target point?                      !  !  !  !
                det = real((-1)**sens, kind=fp) * ( &
                     (uvpoly(1)%mat(2,2) - uvpoly(2)%mat(2,isegments(2,iinter)))*&
                     (uvpoly(2)%mat(1,isegments(2,iinter)+1) - uvpoly(2)%mat(1,isegments(2,iinter))) - &
                     (uvpoly(1)%mat(1,2) - uvpoly(2)%mat(1,isegments(2,iinter)))*&
                     (uvpoly(2)%mat(2,isegments(2,iinter)+1) - uvpoly(2)%mat(2,isegments(2,iinter))) )
                IF ( DEBUG ) THEN
                   PRINT *,'IT       =',IT
                   PRINT *,'XYZ      =',xyztmp
                   PRINT *,'UV       =',uvtmp
                   PRINT *,'IFACE    =',ifacetmp
                   PRINT *,'LAMBDA   =',LAMBDA(:,IINTER)
                   PRINT *,'DET      =',DET

                   PRINT *,'POLYLINE SEGMENT : (SENS =', SENS,')'
                   PRINT *,'XYZ'
                   PRINT *,polyline%xyz(1:3,i1)
                   PRINT *,polyline%xyz(1:3,i2)
                   PRINT *,'UV'
                   PRINT *,polyline%uv(1:2,sens,i1)
                   PRINT *,polyline%uv(1:2,sens,i2)
                END IF
                
                if ( det < -EPSfp ) then
                   changeface = .true.
                   ! else, the target point is inside another brep face                        !  !  !  !
                   ifacetmp = get_face(brep, get_twin(ihedg))                                  !  !  !  !
                   ! get coordinates of the point on the edge's polyline                       !  !  !  !
                   IF ( DEBUG ) THEN
                      PRINT *,'CHANGE BREP FACE'
                      PRINT *,'POLYLINE SEGMENT :'
                      PRINT *,'UV'
                      PRINT *,polyline%uv(1:2,ihedg(2),i1)
                      PRINT *,polyline%uv(1:2,ihedg(2),i2)
                   END IF
                   w = lambda(2,iinter)                                                        !  !  !  !
                   dxyztmp = dxyztmp + xyztmp                                                  !  !  !  !
                   xyztmp = (1._fp - w) * polyline%xyz(1:3,i1) + w * polyline%xyz(1:3,i2)      !  !  !  !
                   !uvtmp = (1._fp - w) * polyline%uv(1:2,ihedg(2),i1) + &                      !  !  !  !
                   !     w * polyline%uv(1:2,ihedg(2),i2)                                       !  !  !  !

                   ! relax on tangential intersection
                   uvtng = (1._fp - w) * polyline%uv(1:2,1:2,i1) + &                           !  !  !  !
                        w * polyline%uv(1:2,1:2,i2)
                   call simultaneous_point_inversions( &
                        brep%edges(ihedg(1))%curve%surf, &
                        lowerb, &
                        upperb, &
                        stat, &
                        uvtng, &
                        xyztmp )
                   IF ( STAT > 0 ) THEN
                      PRINT *,'FAILED TO RELAX ON TANGENTIAL INTERSECTION'
                      CALL WRITE_POLYNOMIAL(brep%edges(ihedg(1))%curve%surf(1)%ptr%x, 'Jouke/debug_relaxtng_surf1.cheb')
                      CALL WRITE_POLYNOMIAL(brep%edges(ihedg(1))%curve%surf(2)%ptr%x, 'Jouke/debug_relaxtng_surf2.cheb')
                      PRINT *,'UVTNG =',uvtng
                      PAUSE 
                   END IF
                   uvtmp = uvtng(1:2,ihedg(2))
                   IF ( maxval(abs(uvtmp)) < 1._fp - EPSuv ) THEN
                      call diffgeom_intersection( &
                           brep%edges(ihedg(1))%curve%surf, &
                           uvtng, &
                           duv_ds, &
                           dxyz_ds, &
                           stat )
                      IF ( STAT > 1 ) THEN
                         PRINT *,'STAT DIFFGEOM_INTERSECTION =',STAT
                         CALL WRITE_POLYNOMIAL(brep%edges(ihedg(1))%curve%surf(1)%ptr%x, 'Jouke/debug_diffgeominter_surf1.cheb')
                         CALL WRITE_POLYNOMIAL(brep%edges(ihedg(1))%curve%surf(2)%ptr%x, 'Jouke/debug_diffgeominter_surf2.cheb')
                         PRINT *,'UV =',UVTNG
                         STOP
                      END IF
                      PRINT *,'DUV_DS =',DUV_DS(:,1,ihedg(2))
                   END IF
                   
                   IF ( DEBUG ) THEN
                      PRINT *,'XYZTMP   =',xyztmp
                      PRINT *,'UVTMP    =',uvtmp
                      PRINT *,'IFACETMP =',ifacetmp
                   END IF
                   ! update xyz displacement                                                   !  !  !  !
                   dxyztmp = dxyztmp - xyztmp                                                  !  !  !  !
                   ! update uv displacement                                                    !  !  !  !
                   do ivar = 1,2 ! <-----------------------+                                   !  !  !  !
                      call evald1( &                       !                                   !  !  !  !
                           jac(:,ivar), &                  !                                   !  !  !  !
                           brep%faces(ifacetmp)%surface, & !                                   !  !  !  !
                           uvtmp, &                        !                                   !  !  !  !
                           ivar )                          !                                   !  !  !  !
                   end do ! <------------------------------+                                   !  !  !  !
                   call solve_NxN( &                                                           !  !  !  !
                        duvtmp, &                                                              !  !  !  !
                        matmul(transpose(jac), jac), &                                         !  !  !  !
                        matmul(transpose(jac), dxyztmp), &                                     !  !  !  !
                        singular )                                                             !  !  !  !
                   dxyztmp = matmul(jac, duvtmp)
                   IF ( DEBUG ) THEN
                      PRINT *,'DXYZTMP  =',dxyztmp
                      PRINT *,'DUVTMP   =',duvtmp
                   END IF
                end if
                !                                                                           !  !  !  !
             end if ! <---------------------------------------------------------------------+  !  !  !
             deallocate(uvpoly(2)%mat)                                                         !  !  !
             if ( changeface ) exit brep_wires
             ! move on to the next halfedge on the wire                                        !  !  !
             ihedg = get_next(brep, ihedg)                                                     !  !  !
             if ( get_orig(brep, ihedg) == istart ) exit brep_halfedges                        !  !  !
             !                                                                                 !  !  !
          end do brep_halfedges ! <------------------------------------------------------------+  !  !
          !                                                                                       !  !
       end do brep_wires ! <----------------------------------------------------------------------+  !
       !                                                                                             !
       if ( .not.changeface ) then
          IF ( DEBUG ) PRINT *,'CONVERGED'
          exit outer_loop                                                                            !
       end if
       it = it + 1                                                                                   !
       IF ( IT > BREP%NF ) THEN
          PRINT *,'*** FAILED TO PROJECT ON HYPERFACE ***'
          DUVTMP(:) = 0._FP
          EXIT outer_loop
       END IF
    end do outer_loop ! <----------------------------------------------------------------------------+

    uvtmp = uvtmp + duvtmp
    deallocate(uvpoly(1)%mat)

  end subroutine projection_hyperface

end module mod_projection
