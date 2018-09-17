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
       uvtmp )
    use mod_linalg
    use mod_types_brep
    use mod_brep
    use mod_halfedge
    use mod_types_intersection
    use mod_intersection
    use mod_tolerances
    use mod_hypergraph
    implicit none
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
    STOP 'projection_hyperedge : failed to reproject onto hyperedge'
    
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
    real(kind=fp)                             :: w, jac(3,2)
    logical                                   :: changeface, singular
    integer                                   :: iwire, ihedg(2), istart, it, ivar, iinter

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
             deallocate(uvpoly(2)%mat)                                                         !  !  !
             !                                                                                 !  !  !
             if ( ninter > 0 ) then ! <-----------------------------------------------------+  !  !  !
                ! pick closest point to temporary uv point                                  !  !  !  !
                iinter = maxloc(lambda(1,1:ninter), 1)                                      !  !  !  !
                IF ( DEBUG ) THEN
                   PRINT *,'IT       =',IT
                   PRINT *,'XYZ      =',xyztmp
                   PRINT *,'UV       =',uvtmp
                   PRINT *,'IFACE    =',ifacetmp
                   PRINT *,'LAMBDA   =',LAMBDA(:,IINTER)
                END IF
                if ( it > 0 .and. lambda(1,iinter) < EPSuv ) exit brep_wires                !  !  !  !
                changeface = .true.                                                         !  !  !  !
                ! new supporting brep face                                                  !  !  !  !
                ifacetmp = get_face(brep, get_twin(ihedg))                                  !  !  !  !
                i1 = i1 + isegments(2,iinter) - 1                                           !  !  !  !
                i2 = i1 + 1                                                                 !  !  !  !
                w = lambda(2,iinter)                                                        !  !  !  !
                ! get coordinates of the point on the edge's polyline                       !  !  !  !
                dxyztmp = dxyztmp + xyztmp                                                  !  !  !  !
                xyztmp = (1._fp - w) * polyline%xyz(1:3,i1) + w * polyline%xyz(1:3,i2)      !  !  !  !
                uvtmp = (1._fp - w) * polyline%uv(1:2,ihedg(2),i1) + &                      !  !  !  !
                     w * polyline%uv(1:2,ihedg(2),i2)                                       !  !  !  !
                ! update xyz displacement                                                   !  !  !  !
                IF ( DEBUG ) THEN
                   PRINT *,'XYZTMP   =',xyztmp
                   PRINT *,'UVTMP    =',uvtmp
                   PRINT *,'IFACETMP =',ifacetmp
                END IF
                dxyztmp = dxyztmp - xyztmp                                                  !  !  !  !
                ! update uv displacement                                                    !  !  !  !
                do ivar = 1,2 ! <-----------------------+                                   !  !  !  !
                   call evald1( &                       !                                   !  !  !  !
                        jac(:,ivar), &                  !                                   !  !  !  !
                        brep%faces(ifacetmp)%surface, & !                                   !  !  !  !
                        uvtmp, &                        !                                   !  !  !  !
                        ivar )                          !                                   !  !  !  !
                end do ! <------------------------------+                                   !  !  !  !
                IF ( DEBUG ) CALL PRINT_MAT(TRANSPOSE(jac))
                call solve_NxN( &                                                           !  !  !  !
                     duvtmp, &                                                              !  !  !  !
                     matmul(transpose(jac), jac), &                                         !  !  !  !
                     matmul(transpose(jac), dxyztmp), &                                     !  !  !  !
                     singular )                                                             !  !  !  !
                IF ( DEBUG ) THEN
                   PRINT *,'DXYZTMP  =',dxyztmp
                   PRINT *,'DUVTMP   =',duvtmp
                END IF
                !                                                                           !  !  !  !
                exit brep_wires                                                             !  !  !  !
             end if ! <---------------------------------------------------------------------+  !  !  !
             ! move on to the next halfedge on the wire                                        !  !  !
             ihedg = get_next(brep, ihedg)                                                     !  !  !
             if ( get_orig(brep, ihedg) == istart ) exit brep_halfedges                        !  !  !
             !                                                                                 !  !  !
          end do brep_halfedges ! <------------------------------------------------------------+  !  !
          !                                                                                       !  !
       end do brep_wires ! <----------------------------------------------------------------------+  !
       !                                                                                             !
       if ( .not.changeface ) exit outer_loop                                                        !
       it = it + 1                                                                                   !
    end do outer_loop ! <----------------------------------------------------------------------------+

    uvtmp = uvtmp + duvtmp
    deallocate(uvpoly(1)%mat)

  end subroutine projection_hyperface

end module mod_projection
