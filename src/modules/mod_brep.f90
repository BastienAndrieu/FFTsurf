module mod_brep

  use mod_math
  use mod_types_brep
  use mod_halfedge

  implicit none

  integer, parameter :: BREP_xtra_nv = 20
  integer, parameter :: BREP_xtra_ne = 20
  integer, parameter :: BREP_xtra_nf = 20
  
contains


  subroutine make_brep( &
       surf, &
       nsurf, &
       interdata, &
       brep, &
       tolchord, &
       hmin, &
       hmax )
    use mod_diffgeom
    use mod_types_intersection
    use mod_types_brep
    use mod_intersection
    implicit none
    integer,                          intent(inout) :: nsurf
    type(type_surface),  allocatable, intent(inout) :: surf(:)
    type(type_intersection_data),     intent(inout) :: interdata
    type(type_brep),                  intent(inout) :: brep
    real(kind=fp),                    intent(in)    :: tolchord, hmin, hmax
    logical                                         :: mask(nsurf)

    mask(1:nsurf) = .true.
    ! compute transversal surface-surface intersections
    call intersect_all_surfaces( &
         surf, &
         nsurf, &
         interdata, &
         mask, &
         tolchord, &
         hmin, &
         hmax )

    IF ( .TRUE. ) THEN
       call write_intersection_data( &
            interdata, &
            '../debug/intersection_points.dat', &
            '../debug/intersection_curves.dat' )
    END IF
    

    ! make BREP
    brep%nv = 0
    brep%ne = 0
    brep%nf = 0
    call make_brep_from_intersection_data( &
         surf, &
         nsurf, &
         interdata, &
         brep )

  end subroutine make_brep

  

  
  subroutine write_brep_files( &
       brep, &
       fileve, &
       filehe, &
       filefa )
    use mod_util
    implicit none
    type(type_BREP), intent(in) :: brep
    character(*),    intent(in) :: fileve, filehe, filefa
    integer                     :: fid, ive, ied, ihe, ifa, iin

    call get_free_unit(fid)

    open(unit=fid, file=fileve, action='write')
    do ive = 1,brep%nv
       write (fid,*) brep%verts(ive)%point%xyz, brep%verts(ive)%halfedge
    end do
    close(fid)

    open(unit=fid, file=filehe, action='write')
    do ied = 1,brep%ne
       do ihe = 1,2
          write (fid,*) &
               brep%edges(ied)%halfedges(ihe)%face, &
               brep%edges(ied)%halfedges(ihe)%orig, &
               brep%edges(ied)%halfedges(ihe)%prev, &
               brep%edges(ied)%halfedges(ihe)%next
       end do
    end do
    close(fid)

    open(unit=fid, file=filefa, action='write')
    do ifa = 1,brep%nf
       write (fid,*) brep%faces(ifa)%outer
       write (fid,*) brep%faces(ifa)%ninner
       do iin = 1,brep%faces(ifa)%ninner
          write (fid,*) brep%faces(ifa)%inner(:,iin)
       end do
    end do
    close(fid)
    
  end subroutine write_brep_files
  


  

































  


  subroutine generate_face_mesh( &       
       brep, &
       iface, &
       filec, &
       filebpts, &
       filebedg, &
       fileinfo, &
       filetri, &
       fileuv, &
       filexyz )
    use mod_util
    use mod_polynomial
    implicit none
    type(type_BREP), intent(in) :: brep
    integer,         intent(in) :: iface
    character(*),    intent(in) :: filec, filebpts, filebedg, fileinfo, filetri, fileuv, filexyz
    integer                     :: npts, nedg
    real(kind=fp), allocatable  :: uv(:,:)
    integer, allocatable        :: edg(:,:)
    integer                     :: iwire, fid, ipt, ied
    
    ! write polynomial parametric surface 
    call write_polynomial(brep%faces(iface)%surface%x, filec)
    
    ! get boundary points & edges
    allocate(uv(2,1000), edg(2,1000))
    npts = 0
    nedg = 0
    call get_wire_uv_polyline( &
         brep, &
         brep%faces(iface)%outer, &
         npts, &
         uv, &
         nedg, &
         edg )
    do iwire = 1,brep%faces(iface)%ninner
       call get_wire_uv_polyline( &
         brep, &
         brep%faces(iface)%inner(:,iwire), &
         npts, &
         uv, &
         nedg, &
         edg )
    end do

    call get_free_unit(fid)
    ! write boundary points
    open(unit=fid, file=filebpts, action='write')
    !write (fid,*) npts, 2
    do ipt = 1,npts
       write (fid,*) uv(:,ipt)
    end do
    close(fid)

    ! write boundary edges
    open(unit=fid, file=filebedg, action='write')
    !write (fid,*) nedg
    do ied = 1,nedg
       write (fid,*) edg(:,ied)
    end do
    close(fid)

    ! run surface mesher
    call system('/stck/bandrieu/Bureau/MeshGen/./meshgen.out &
         & '//trim(filec)//'&
         & '//trim(filebpts)//'&
         & '//trim(filebedg)//'&
         & '//trim(fileinfo)//'&
         & '//trim(filetri)//'&
         & '//trim(fileuv)//'&
         & '//trim(filexyz))
    
  end subroutine generate_face_mesh


  
  subroutine get_wire_uv_polyline( &
       brep, &
       ifirstedge, &
       npts, &
       uv, &
       nedg, &
       edg )
    implicit none
    type(type_BREP),             intent(in)    :: brep
    integer,                     intent(in)    :: ifirstedge(2)
    integer,                     intent(inout) :: npts
    real(kind=fp), allocatable,  intent(inout) :: uv(:,:)
    integer,                     intent(inout) :: nedg
    integer, allocatable,        intent(inout) :: edg(:,:)
    real(kind=fp), allocatable                 :: uvtmp(:,:), uvedge(:,:)
    integer, allocatable                       :: edgtmp(:,:)
    integer                                    :: ihedg(2), nedgprev, npedge, i

    nedgprev = nedg
    ihedg = ifirstedge
    do ! <--------------------------------------------------------------+
       call get_uv_polyline( &                                          !
            brep, &                                                     !
            ihedg, &                                                    !
            npedge, &                                                   !
            uvedge )                                                    !
       !                                                                !
       if ( npts + npedge - 1 > size(uv,2) ) then ! <------+            !
          call move_alloc(from=uv, to=uvtmp)               !            !
          allocate(uv(2,npts + npedge + 100))              !            !
          uv(1:2,1:npts) = uvtmp(1:2,1:npts)               !            !
          deallocate(uvtmp)                                !            !
       end if ! <------------------------------------------+            !
       uv(1:2,npts+1:npts+npedge-1) = uvedge(1:2,1:npedge-1)            !
       !                                                                !
       if ( npts + npedge - 1 > size(edg,2) ) then ! <-----+            !
          call move_alloc(from=edg, to=edgtmp)             !            !
          allocate(edg(2,nedg + npedge + 100))             !            !
          edg(1:2,1:nedg) = edgtmp(1:2,1:nedg)             !            !
          deallocate(edgtmp)                               !            !
       end if ! <------------------------------------------+            !
       edg(1,nedg+1:nedg+npedge-1) = npts + [(i, i=1,npedge-1)]         !
       edg(2,nedg+1:nedg+npedge-1) = edg(1,nedg+1:nedg+npedge-1) + 1    !
       !                                                                !
       npts = npts + npedge - 1                                         !
       nedg = nedg + npedge - 1                                         !
       !                                                                !
       ihedg = get_next(brep, ihedg)                                    !
       if ( all(ihedg - ifirstedge == 0) ) then ! <----------+          !
          edg(2,nedg) = edg(1,nedgprev+1)                    !          !
          return                                             !          !
       end if ! <--------------------------------------------+          !
    end do ! <----------------------------------------------------------+
    
  end subroutine get_wire_uv_polyline
  



  
  subroutine get_uv_polyline( &
       brep, &
       ihedg, &
       np, &
       uv )
    implicit none
    type(type_BREP),             intent(in)    :: brep
    integer,                     intent(in)    :: ihedg(2)
    integer,                     intent(out)   :: np
    real(kind=fp), allocatable,  intent(inout) :: uv(:,:)
    integer                                    :: sens, ifirst, ilast

    call get_polyline_endpoints( &
         brep, &
         ihedg, &
         ifirst, &
         ilast, &
         sens, &
         np )

    if ( allocated(uv) .and. size(uv,2) < np ) deallocate(uv)
    if ( .not.allocated(uv) ) allocate(uv(2,np))
    uv(1:2,1:np) = brep%edges(ihedg(1))%curve%polyline%uv(:,sens,ifirst:ilast:(-1)**sens)
    
  end subroutine get_uv_polyline


  
  

  subroutine make_brep_from_intersection_data( &
       surf, &
       nsurf, &
       interdata, &
       brep )
    USE MOD_UTIL
    use mod_graph
    use mod_intersection
    use mod_tolerances
    implicit none
    LOGICAL, PARAMETER :: DEBUG = .false.
    integer,                      intent(in)    :: nsurf
    type(type_surface), target,   intent(in)    :: surf(nsurf)
    type(type_intersection_data), intent(inout) :: interdata
    type(type_BREP),              intent(inout) :: brep
    integer                                     :: stat
    integer, dimension(interdata%nc)            :: arc2curve, arc2split, arc_sens
    integer                                     :: arc2nod(2,interdata%nc)
    real(kind=fp)                               :: arc_angles(2,interdata%nc)
    integer                                     :: point2nod(interdata%np), nod2point(interdata%np)
    real(kind=fp)                               :: nod_uv(2,interdata%np)
    integer                                     :: narc, nnod
    integer                                     :: nwires
    integer, allocatable                        :: wire2arc(:,:), lenwire(:)
    logical                                     :: inside
    integer, allocatable                        :: nsuper(:), super(:,:)
    integer                                     :: nfaces
    integer, allocatable                        :: outer(:), inner(:,:), ninner(:)
    integer                                     :: isurf, iwire, jwire, iface, jface, iedge, ihedg, ivert
    INTEGER :: FID, IARC, INOD
    CHARACTER(3) :: STRNUM

    IF ( DEBUG ) PRINT *,'>>>> make_brep_from_intersection_data'
    
    do isurf = 1,nsurf
       IF ( DEBUG ) THEN
          PRINT *,''
          PRINT *,'ISURF =',ISURF
       END IF
       ! build an embedded graph of intersection curves & points
       call build_intersections_graph( &
            surf(isurf), &
            interdata, &
            arc2curve, &
            arc2split, &
            arc_sens, &
            arc_angles, &
            arc2nod, &
            narc, &
            nod2point, &
            point2nod, &
            nod_uv, &
            nnod, &
            stat )
       
       IF ( DEBUG ) THEN
          WRITE (STRNUM, '(I3.3)') ISURF
          CALL GET_FREE_UNIT(FID)
          OPEN(UNIT=FID, FILE='/d/bandrieu/GitHub/FFTsurf/test/Jouke/graph/graph_'//strnum//'.dat', ACTION='WRITE')
          WRITE (FID,*) NARC
          DO IARC = 1,NARC
             WRITE (FID,*) ARC2NOD(:,IARC)
          END DO
          DO IARC = 1,NARC
             WRITE (FID,*) ARC_ANGLES(:,IARC)
          END DO
          WRITE (FID,*) NNOD
          DO INOD = 1,NNOD
             WRITE (FID,*) NOD_UV(:,INOD)
          END DO
          CLOSE(FID)
       END IF

       IF ( STAT > 0 ) then
          PRINT *,'ISURF =',ISURF
          PRINT *,'STAT =',STAT
          STOP ':('
       end IF

       ! make wires
       nwires = min(nnod, narc)
       allocate(lenwire(nwires), wire2arc(nwires,nwires))
       call make_wires( &
            arc2nod(1:2,1:narc), &
            arc_angles(1:2,1:narc), &
            narc, &
            nnod, &
            nwires, &
            wire2arc, &
            lenwire )

       IF ( DEBUG ) THEN
          WRITE (STRNUM, '(I3.3)') ISURF
          CALL GET_FREE_UNIT(FID)
          OPEN(UNIT=FID, FILE='/d/bandrieu/GitHub/FFTsurf/test/Jouke/graph/wires_'//strnum//'.dat', ACTION='WRITE')
          WRITE (FID,*) NWIRES
          DO IWIRE = 1,NWIRES
             WRITE (FID,*) WIRE2ARC(1:LENWIRE(IWIRE),IWIRE)
             WRITE (FID,*) ARC2NOD(1,WIRE2ARC(1:LENWIRE(IWIRE),IWIRE))!WIRENOD(1:LENWIRE(IWIRE),IWIRE)
          END DO
          CLOSE(FID)
       END IF

       if ( nwires < 1 ) then
          deallocate(lenwire, wire2arc)
          cycle
       end if

       ! compute inclusion relationships between nested wires
       allocate(nsuper(nwires), super(nwires,nwires))
       nsuper(:) = 0
       do iwire = 1,nwires-1 ! <---------------------------------------------+
          do jwire = iwire+1,nwires ! <----------------------------------+   !
             ! get a point on wire #jwire, and test whether it is        !   !
             ! inside the polygon of wire #iwire                         !   !
             call point_in_wire( &                                       !   !
                  nod_uv(:,arc2nod(1,wire2arc(1,jwire))), &              !   !
                  interdata, &                                           !   !
                  lenwire(iwire), &                                      !   !
                  arc2curve(wire2arc(1:lenwire(iwire),iwire)), &         !   !
                  arc2split(wire2arc(1:lenwire(iwire),iwire)), &         !   !
                  arc_sens(wire2arc(1:lenwire(iwire),iwire)), &          !   !
                  inside )                                               !   !
             if ( inside ) then ! <----------------------------------+   !   !
                ! wire #jwire is nested inside wire #iwire           !   !   !
                nsuper(jwire) = nsuper(jwire) + 1                    !   !   !
                super(nsuper(jwire),jwire) = iwire                   !   !   !
                cycle                                                !   !   !
             else ! -------------------------------------------------+   !   !
                ! get a point on wire #iwire, and test whether it is !   !   !
                ! inside the polygon of wire #jwire                  !   !   !
                call point_in_wire( &                                !   !   !
                     nod_uv(:,arc2nod(1,wire2arc(1,iwire))), &       !   !   !
                     interdata, &                                    !   !   !
                     lenwire(jwire), &                               !   !   !
                     arc2curve(wire2arc(1:lenwire(jwire),jwire)), &  !   !   !
                     arc2split(wire2arc(1:lenwire(jwire),jwire)), &  !   !   !
                     arc_sens(wire2arc(1:lenwire(jwire),jwire)), &   !   !   !
                     inside )                                        !   !   !
                if ( inside ) then ! <---------------------------+   !   !   !
                   ! wire #iwire is nested inside wire #jwire    !   !   !   !
                   nsuper(iwire) = nsuper(iwire) + 1             !   !   !   !
                   super(nsuper(iwire),iwire) = jwire            !   !   !   !
                end if ! <---------------------------------------+   !   !   !
             end if ! <----------------------------------------------+   !   !
          end do ! <-----------------------------------------------------+   !
       end do ! <------------------------------------------------------------+

       !! Make faces
       allocate(outer(nwires), inner(nwires,nwires), ninner(nwires))
       call make_faces( &                                           
            nsuper(1:nwires), &                                
            super(1:nwires,1:nwires), &                        
            nwires, &                                               
            outer, &                                                
            inner, &                                                
            ninner, &                                               
            nfaces )

       IF ( DEBUG ) THEN
          WRITE (STRNUM, '(I3.3)') ISURF
          CALL GET_FREE_UNIT(FID)
          OPEN(UNIT=FID, FILE='/d/bandrieu/GitHub/FFTsurf/test/Jouke/graph/faces_'//strnum//'.dat', ACTION='WRITE')
          WRITE (FID,*) NFACES
          DO IFACE = 1,NFACES
             WRITE (FID,*) WIRE2ARC(1:LENWIRE(OUTER(IFACE)),OUTER(IFACE))
             WRITE (FID,*) NINNER(IFACE)
             DO IWIRE = 1,NINNER(IFACE)
                WRITE (FID,*) WIRE2ARC(1:LENWIRE(INNER(IWIRE,IFACE)),INNER(IWIRE,IFACE))
             END DO
          END DO
          CLOSE(FID)
       END IF

       !! Insert faces in BREP
       IF ( DEBUG ) THEN
          PRINT *,'   NWIRES =',NWIRES
          PRINT *,'   NINNER =',NINNER
          PRINT *,'   NFACES =',NFACES
       END IF

       do jface = 1,nfaces ! <-------------------------------+
          ! + 1 BREPface                                     !
          iface = brep%nf + 1                                !
          !                                                  !
          if ( .not.allocated(brep%faces) .or. &             !
               iface > size(brep%faces) ) then ! <---+       !
             call reallocate_faces( &                !       !
                  brep, &                            !       !
                  iface + BREP_xtra_nf )             !       !
          end if ! <---------------------------------+       !
          !                                                  !
          brep%nf = brep%nf + 1                              !
          brep%faces(brep%nf)%surface => surf(isurf)         !
          !                                                  !
          ! add outer wire                                   !
          iwire = outer(jface)                               !
          call insert_wire_in_brep( &                        !
               interdata, &                                  !
               brep, &                                       !
               iface, &                                      !
               .true., &                                     !
               wire2arc(1:lenwire(iwire),iwire), &           !
               lenwire(iwire), &                             !
               arc2curve(1:narc), &                          !
               arc2split(1:narc), &                          !
               arc2nod(1:2,1:narc), &                        !
               arc_sens(1:narc), &                           !
               narc, &                                       !
               nod2point(1:nnod), &                          !
               nnod )                                        !
          !                                                  !
          ! add potential inner wires                        !
          do jwire = 1,ninner(jface) ! <-----------------+   !
             iwire = inner(jwire,jface)                  !   !
             call insert_wire_in_brep( &                 !   !
                  interdata, &                           !   !
                  brep, &                                !   !
                  iface, &                               !   !
                  .false., &                             !   !
                  wire2arc(1:lenwire(iwire),iwire), &    !   !
                  lenwire(iwire), &                      !   !
                  arc2curve(1:narc), &                   !   !
                  arc2split(1:narc), &                   !   !
                  arc2nod(1:2,1:narc), &                 !   !
                  arc_sens(1:narc), &                    !   !
                  narc, &                                !   !
                  nod2point(1:nnod), &                   !   !
                  nnod )                                 !   !
          end do ! <-------------------------------------+   !
          !                                                  !
       end do ! <--------------------------------------------+

       deallocate(lenwire, wire2arc, nsuper, super, outer, inner, ninner)
    end do
    

    ! finally, set halfedge record for BREPvertices
    PRINT *,'BREP%NE =',BREP%NE
    PRINT *,'BREP%NF =',BREP%NF
    PRINT *,'BREP%NV =',BREP%NV
    do ivert = 1,brep%nv
       brep%verts(ivert)%halfedge(:) = 0
    end do
    
    do iedge = 1,brep%ne ! <--------------------------------------------+
       do ihedg = 1,2 ! <-------------------------------------------+   !
          if ( brep%edges(iedge)%halfedges(ihedg)%face < 1 ) cycle  !   !
          ivert = brep%edges(iedge)%halfedges(ihedg)%orig           !   !
          if ( brep%verts(ivert)%halfedge(1) > 0 ) then ! <---+     !   !
             if ( .not.is_boundary_edge(brep, iedge) ) cycle  !     !   !
          end if ! <------------------------------------------+     !   !
          !                                                         !   !
          brep%verts(ivert)%halfedge = [iedge, ihedg]               !   !
       end do ! <---------------------------------------------------+   !
    end do ! <----------------------------------------------------------+

    IF ( DEBUG ) PRINT *,'make_brep_from_intersection_data >>>>'
    
  end subroutine make_brep_from_intersection_data
  




  
  

  subroutine build_intersections_graph( &
       surf, &
       interdata, &
       arc2curve, &
       arc2split, &
       arc_sens, &
       arc_angles, &
       arc2nod, &
       narc, &
       nod2point, &
       point2nod, &
       nod_uv, &
       nnod, &
       stat )
    use mod_diffgeom
    use mod_types_intersection
    use mod_tolerances
    implicit none
    type(type_surface), target,   intent(in)    :: surf
    type(type_intersection_data), target, intent(inout) :: interdata
    integer,                      intent(inout) :: arc2curve(:)
    integer,                      intent(inout) :: arc2split(:)
    integer,                      intent(inout) :: arc_sens(:)
    real(kind=fp),                intent(inout) :: arc_angles(:,:)
    integer,                      intent(inout) :: arc2nod(:,:)
    integer,                      intent(out)   :: narc
    integer,                      intent(inout) :: nod2point(:)
    integer,                      intent(out)   :: point2nod(interdata%np)
    real(kind=fp),                intent(inout) :: nod_uv(:,:)
    integer,                      intent(out)   :: nnod
    integer,                      intent(out)   :: stat
    type(type_intersection_curve), pointer      :: curve => null()
    real(kind=fp)                               :: duv_ds(2)
    integer                                     :: icurv, jsurf
    integer                                     :: isplit, jsplit, ksplit
    integer                                     :: ipoint, jpoint, jnod

    narc = 0
    nnod = 0
    point2nod(:) = 0

    stat = 0
    do icurv = 1,interdata%nc ! <-----------------------------------------------------+
       curve => interdata%curves(icurv)
       do jsurf = 1,2 ! <---------------------------------------------------------+   !
          if ( associated(curve%surf(jsurf)%ptr, surf) ) exit   !   !
       end do ! <-----------------------------------------------------------------+   !
       if ( jsurf > 2 ) cycle                                                         !
       !                                                                              !
       !PRINT *,'ICURV, JSURF =', icurv, jsurf
       !! Add a node for each split point                                             !
       ! set the first split point index:                                             !
       !   if jsurf == 1, traverse the split points in reverse order                  !
       !   "         " 2, "                          " direct order                   !
       isplit = 1 + (2 - jsurf)*(curve%nsplit - 1)                  !
       do jsplit = 1,curve%nsplit ! <---------------------+         !
          ! get intersection point index                                    !         !
          ipoint = curve%isplit(1,isplit)                 !         !
          !PRINT *,'JSPLIT, ISPLIT, IPOINT =', jsplit, isplit, ipoint
          if ( point2nod(ipoint) == 0 ) then ! <------------------------+   !         !
             ! add a new node                                           !   !         !
             nnod = nnod + 1                                            !   !         !
             ! assign it to the intersection point...                   !   !         !
             nod2point(nnod) = ipoint                                   !   !         !
             ! ...and vice versa                                        !   !         !
             point2nod(ipoint) = nnod                                   !   !         !
             ! get polyline point index                                 !   !         !
             ipoint = curve%isplit(2,isplit)          !   !         !
             ! get coordinates of the embedded node                     !   !         !
             IF ( ipoint < 1 .OR. &
                  ipoint > curve%polyline%NP ) THEN
                STAT = 1
                RETURN
             END IF
             nod_uv(1:2,nnod) = curve%polyline%uv&    !   !         !
                  (1:2,jsurf,ipoint)                                    !   !         !
          end if ! <----------------------------------------------------+   !         !
          ! move on to the next split point                                 !         !
          isplit = isplit + (-1)**jsurf                                     !         !
       end do ! <-----------------------------------------------------------+         !
       !                                                                              !
       !! Add an arc for each curve segment (between 2 split points)                  !
       isplit = 1 + (2 - jsurf)*(curve%nsplit - 2)                  !
       do jsplit = 1,curve%nsplit - 1 ! <-----------------+         !
          !PRINT *,'JSPLIT, ISPLIT =', jsplit, isplit
          ! add a new arc                                                   !         !
          narc = narc + 1                                                   !         !
          ! assign it to the intersection curve segment                     !         !
          arc2curve(narc) = icurv                                           !         !
          arc2split(narc) = isplit                                          !         !
          ! set orientation                                                 !         !
          arc_sens(narc) = jsurf                                            !         !
          ! set endpoints data                                              !         !
          ksplit = isplit                                                   !         !
          do jnod = 1,2 ! <---------------------------------------------+   !         !
             ! get intersection point index                             !   !         !
             ipoint = curve%isplit(1,ksplit+2-jsurf)                    !   !         !
             !PRINT *,'JNOD, IPOINT =', jnod, ipoint
             ! set endpoint node index                                  !   !         !
             arc2nod(jnod,narc) = point2nod(ipoint)                     !   !         !
             ! get polyline point index                                 !   !         !
             ipoint = curve%isplit(2,ksplit+2-jsurf)                    !   !         !
             jpoint = ipoint + (-1)**(jsurf + jnod + 1)            
             !PRINT *,'IPOLYPOINT, JPOLYPOINT, NP =', ipoint, jpoint, curve%polyline%np
             IF ( ipoint < 1 .OR. &
                  ipoint > curve%polyline%NP .OR. &
                  jpoint < 1 .OR. &
                  jpoint > curve%polyline%NP ) THEN
                STAT = 2
                RETURN
             END IF
             ! compute tangent to the embedded arc at current endpoint  !   !         !
             do
                duv_ds = &
                     curve%polyline%uv(1:2,jsurf,jpoint) - &
                     curve%polyline%uv(1:2,jsurf,ipoint)
                if ( sum(duv_ds**2) > EPSuvsqr ) then
                   exit
                else
                   jpoint = jpoint + (-1)**(jsurf + jnod + 1)
                   !PRINT *,'|DUV_DS| << 1, JPOLYPOINT =',jpoint
                   !PAUSE
                end if
             end do
             duv_ds = real((-1)**(jnod+1), kind=fp) * duv_ds
             !PRINT *,'IARC, INOD, DUV_DS =', narc, arc2nod(jnod,narc), duv_ds
             !duv_ds = &                                                 !   !         !
             !     curve%polyline%uv&                                    !   !         !
             !     (:,jsurf,ipoint + (2-jnod)*(-1)**jsurf) - &           !   !         !
             !     curve%polyline%uv&                                    !   !         !
             !     (:,jsurf,ipoint + (1-jnod)*(-1)**jsurf)               !   !         !
             ! convert to angle                                         !   !         !
             arc_angles(jnod,narc) = atan2(duv_ds(2), duv_ds(1))        !   !         !
             ksplit = ksplit + (-1)**jsurf                              !   !         !
          end do ! <----------------------------------------------------+   !         !
          ! move on to the next split point                                 !         !
          isplit = isplit + (-1)**jsurf                                     !         !
       end do ! <-----------------------------------------------------------+         !
       !                                                                              !
    end do ! <------------------------------------------------------------------------+

  end subroutine build_intersections_graph




  
  
  

  subroutine insert_wire_in_brep( &
       interdata, &
       brep, &
       iface, &
       is_outer_wire, &
       wire2arc, &
       lenwire, &
       arc2curve, &
       arc2split, &
       arc2nod, &
       arc_sens, &
       narc, &
       nod2point, &
       nnod )
    use mod_util
    use mod_types_intersection
    implicit none
    type(type_intersection_data), target, intent(inout) :: interdata
    type(type_BREP),                      intent(inout) :: brep
    integer,                              intent(in)    :: iface
    logical,                              intent(in)    :: is_outer_wire
    integer,                              intent(in)    :: lenwire
    integer,                              intent(in)    :: wire2arc(lenwire)
    integer,                              intent(in)    :: narc
    integer,                              intent(in)    :: arc2curve(narc)
    integer,                              intent(in)    :: arc2split(narc)
    integer,                              intent(in)    :: arc2nod(2,narc)
    integer,                              intent(in)    :: arc_sens(narc)
    integer,                              intent(in)    :: nnod
    integer,                              intent(in)    :: nod2point(nnod)
    integer                                             :: ihedg(2,lenwire)
    integer                                             :: iarc, jarc, icurve, isplit, ipoint, ivert

    do jarc = 1,lenwire ! <------------------------------------------------+
       iarc = wire2arc(jarc)                                               !
       !                                                                   !
       icurve = arc2curve(iarc)                                            !
       isplit = arc2split(iarc)                                            !
       !                                                                   !
       ihedg(1,jarc) = interdata%curves(icurve)%iedge(isplit)              !
       ihedg(2,jarc) = 1 + mod(arc_sens(iarc),2)                           !
       !                                                                   !
       if ( ihedg(1,jarc) == 0 ) then ! <------------------------------+   !
          ! +1 BREPedge                                                !   !
          ihedg(1,jarc) = brep%ne + 1                                  !   !
          interdata%curves(icurve)%iedge(isplit) = ihedg(1,jarc)       !   !
          !                                                            !   !
          if ( .not.allocated(brep%edges) .or. &                       !   !
               ihedg(1,jarc) > size(brep%edges) ) then ! <---+         !   !
             call reallocate_edges( &                        !         !   !
                  brep, &                                    !         !   !
                  ihedg(1,jarc) + BREP_xtra_ne )             !         !   !
          end if ! <-----------------------------------------+         !   !
          brep%ne = brep%ne + 1                                        !   !
          ! set pointer to intersection curve segment                  !   !
          brep%edges(ihedg(1,jarc))%curve => interdata%curves(icurve)  !   !
          brep%edges(ihedg(1,jarc))%isplit = isplit                    !   !
       end if ! <------------------------------------------------------+   !
       !                                                                   !
       brep%edges(ihedg(1,jarc))%halfedges(ihedg(2,jarc))%face = iface     !
       !                                                                   !
       ipoint = nod2point(arc2nod(1,iarc))                                 !
       ivert = interdata%points(ipoint)%ivert                              !
       if ( ivert == 0 ) then ! <-------------------------------+          !
          ! +1 BREPvertex                                       !          !
          ivert = brep%nv + 1                                   !          !
          if ( .not.allocated(brep%verts) .or. &                !          !
               ivert > size(brep%verts) ) then ! <---+          !          !
             call reallocate_vertices( &             !          !          !
                  brep, &                            !          !          !
                  ivert + BREP_xtra_nv )             !          !          !
          end if ! <---------------------------------+          !          !
          interdata%points(ipoint)%ivert = ivert                !          !
          brep%nv = brep%nv + 1                                 !          !
          ! set pointer to intersection point                   !          !
          brep%verts(ivert)%point => interdata%points(ipoint)   !          !
       end if ! <-----------------------------------------------+          !
       !                                                                   !
       brep%edges(ihedg(1,jarc))%halfedges(ihedg(2,jarc))%orig = ivert     !
    end do ! <-------------------------------------------------------------+

    ! set prev/next records for new halfedges
    do jarc = 1,lenwire ! <-----------------------------------------+
       brep%edges(ihedg(1,jarc))%halfedges(ihedg(2,jarc))%prev = &  !
            ihedg(:,1 + mod(jarc + lenwire - 2,lenwire))            !
       brep%edges(ihedg(1,jarc))%halfedges(ihedg(2,jarc))%next = &  !
            ihedg(:,1 + mod(jarc,lenwire))                          !
    end do ! <------------------------------------------------------+

    ! set outer/inner record for current face
    if ( is_outer_wire ) then ! <-----------+
       brep%faces(iface)%outer = ihedg(:,1) !
    else ! ---------------------------------+
       call insert_column_after( &          !
            brep%faces(iface)%inner, &      !
            2, &                            !
            brep%faces(iface)%ninner, &     !
            ihedg(:,1), &                   !
            brep%faces(iface)%ninner )      !
    end if ! <------------------------------+
    
  end subroutine insert_wire_in_brep



  
  
  subroutine point_in_wire( &
       uv, &
       interdata, &
       lenwire, &
       arc2curve, &
       arc2split, &
       arc_sens, &
       inside )
    use mod_types_intersection
    ! tests whether a point lies inside the region outlined by a wire of
    ! intersection curves
    implicit none
    real(kind=fp), parameter                  :: M = 10._fp
    real(kind=fp),                intent(in)  :: uv(2)
    type(type_intersection_data), intent(in)  :: interdata
    integer,                      intent(in)  :: lenwire
    integer,                      intent(in)  :: arc2curve(lenwire)
    integer,                      intent(in)  :: arc2split(lenwire)
    integer,                      intent(in)  :: arc_sens(lenwire)
    logical,                      intent(out) :: inside
    real(kind=fp), pointer                    :: poly(:,:) => null()
    real(kind=fp)                             :: a, c, s, ub, vb
    integer                                   :: np
    integer                                   :: iarc, sens, icurv, isplit, ifirst, ilast, i

    inside = .false.

    call random_number(a)
    a = 0.5_fp * CSTpi * a
    c = cos(a)
    s = sin(a)
    ub = uv(1) + M*c
    vb = uv(2) + M*s

    do iarc = 1,lenwire
       sens   = arc_sens(iarc)
       icurv  = arc2curve(iarc)
       isplit = arc2split(iarc) + 2 - sens

       ifirst = interdata%curves(icurv)%isplit(2,isplit)
       ilast  = interdata%curves(icurv)%isplit(2,isplit + (-1)**sens)
       np = (-1)**sens * (ilast - ifirst) + 1

       poly => interdata%curves(icurv)%polyline%uv(:,sens,ifirst:ilast:(-1)**sens)
       do i = 1,np-1
          if ( ( ((poly(1,i) - uv(1))*s > (poly(2,i) - uv(2))*c) .neqv. &
               ((poly(1,i+1) - uv(1))*s > (poly(2,i+1) - uv(2))*c) ) .and. &
               ( ((uv(1) - poly(1,i))*(poly(2,i) - poly(2,i+1)) > &
               (uv(2) - poly(2,i))*(poly(1,i) - poly(1,i+1))) .neqv. &
               ((ub - poly(1,i))*(poly(2,i) - poly(2,i+1)) > &
               (vb - poly(2,i))*(poly(1,i) - poly(1,i+1))) ) ) then
             inside = .not.inside
          end if
       end do

       nullify(poly)
    end do
    
  end subroutine point_in_wire





  ! =============== (HALF)EDGE QUERIES ===============
  function is_smooth( &
       brep, &
       iedge )
    implicit none
    type(type_BREP), intent(in) :: brep
    integer,         intent(in) :: iedge
    logical                     :: is_smooth

    is_smooth = brep%edges(iedge)%curve%smooth

  end function is_smooth
  
  function is_boundary_edge( &
       brep, &
       iedge )
    implicit none
    type(type_BREP), intent(in) :: brep
    integer,         intent(in) :: iedge
    logical                     :: is_boundary_edge

    is_boundary_edge = &
         ( brep%edges(iedge)%halfedges(1)%face == 0 ) .or. &
         ( brep%edges(iedge)%halfedges(2)%face == 0 )
    
  end function is_boundary_edge


  subroutine get_polyline_endpoints( &
       brep, &
       ihedg, &
       ifirst, &
       ilast, &
       sens, &
       np )
    implicit none
    type(type_BREP), intent(in)  :: brep
    integer,         intent(in)  :: ihedg(2)
    integer,         intent(out) :: ifirst
    integer,         intent(out) :: ilast
    integer,         intent(out) :: sens
    integer,         intent(out) :: np
    integer                      :: isplit

    isplit = brep%edges(ihedg(1))%isplit

    sens = 1 + mod(ihedg(2),2)
    ifirst = brep%edges(ihedg(1))%curve%isplit(2,isplit + 2 - sens)
    ilast  = brep%edges(ihedg(1))%curve%isplit(2,isplit + 2 - sens + (-1)**sens)
    np = (-1)**sens * (ilast - ifirst) + 1
    
  end subroutine get_polyline_endpoints


  function get_polyline_length( &
       brep, &
       ihedg )
    implicit none
    type(type_BREP), intent(in)  :: brep
    integer,         intent(in)  :: ihedg(2)
    real(kind=fp)                :: get_polyline_length
    integer                      :: ifirst, ilast, sens, np
    
    call get_polyline_endpoints( &
       brep, &
       ihedg, &
       ifirst, &
       ilast, &
       sens, &
       np )
    get_polyline_length = real((-1)**sens, kind=fp) * ( &
         brep%edges(ihedg(1))%curve%polyline%s(ilast) - &
         brep%edges(ihedg(1))%curve%polyline%s(ifirst) )
    
  end function get_polyline_length
  ! ==================================================

  


  subroutine get_v2f( &
       brep, &
       ivert, &
       faces, &
       nfaces )
    ! returns all faces incident to a vertex, traversed CCW
    use mod_util
    implicit none
    type(type_BREP),      intent(in)    :: brep
    integer,              intent(in)    :: ivert
    integer, allocatable, intent(inout) :: faces(:)
    integer,              intent(out)   :: nfaces
    integer                             :: ihedg(2), iface

    ihedg = brep%verts(ivert)%halfedge
    iface = get_face(brep, ihedg)

    nfaces = 0
    do
       call insert_after( &
            faces, &
            nfaces, &
            iface, &
            nfaces )

       ! traverse halfedges counter-clockwise
       ihedg = get_prev(brep, ihedg) ! previous halfedge
       ihedg = get_twin(ihedg) ! twin halfedge
       iface = get_face(brep, ihedg)

       if ( iface < 1 .or. iface == faces(1) ) return
    end do
    
  end subroutine get_v2f



  


  subroutine reallocate_vertices( &
       brep, &
       nv )
    implicit none
    integer,             intent(in)    :: nv
    type(type_BREP),     intent(inout) :: brep
    type(type_BREPvertex), allocatable :: tmp(:)

    if ( .not.allocated(brep%verts) ) then
       allocate(brep%verts(nv))
       return
    end if

    allocate(tmp(brep%nv))
    call transfer_vertices( &
       from=brep%verts, &
       to=tmp, &
       n=brep%nv )
    deallocate(brep%verts)
    
    allocate(brep%verts(nv))
    call transfer_vertices( &
       from=tmp, &
       to=brep%verts, &
       n=brep%nv )
    deallocate(tmp)

  end subroutine reallocate_vertices



  subroutine reallocate_edges( &
       brep, &
       ne )
    implicit none
    integer,           intent(in)    :: ne
    type(type_BREP),   intent(inout) :: brep
    type(type_BREPedge), allocatable :: tmp(:)

    if ( .not.allocated(brep%edges) ) then
       allocate(brep%edges(ne))
       return
    end if

    allocate(tmp(brep%ne))
    call transfer_edges( &
       from=brep%edges, &
       to=tmp, &
       n=brep%ne )
    deallocate(brep%edges)
    
    allocate(brep%edges(ne))
    call transfer_edges( &
       from=tmp, &
       to=brep%edges, &
       n=brep%ne )
    deallocate(tmp)
    
  end subroutine reallocate_edges


 subroutine reallocate_faces( &
       brep, &
       nf )
    implicit none
    integer,           intent(in)    :: nf
    type(type_BREP),   intent(inout) :: brep
    type(type_BREPface), allocatable :: tmp(:)

    if ( .not.allocated(brep%faces) ) then
       allocate(brep%faces(nf))
       return
    end if

    allocate(tmp(brep%nf))
    call transfer_faces( &
       from=brep%faces, &
       to=tmp, &
       n=brep%nf )
    deallocate(brep%faces)
    
    allocate(brep%faces(nf))
    call transfer_faces( &
       from=tmp, &
       to=brep%faces, &
       n=brep%nf )
    deallocate(tmp)
    
  end subroutine reallocate_faces
  


  
  subroutine transfer_vertices( &
       from, &
       to, &
       n )
    implicit none
    integer,               intent(in)    :: n
    type(type_BREPvertex), intent(inout) :: from(:)
    type(type_BREPvertex), intent(inout) :: to(:)
    integer                              :: i

    do i = 1,n
       to(i)%point => from(i)%point
       nullify(from(i)%point)
       to(i)%halfedge  =  from(i)%halfedge
    end do

  end subroutine transfer_vertices
  
  

  subroutine transfer_edges( &
       from, &
       to, &
       n )
    implicit none
    integer,             intent(in)    :: n
    type(type_BREPedge), intent(inout) :: from(:)
    type(type_BREPedge), intent(inout) :: to(:)
    integer                            :: i

    do i = 1,n
       to(i)%curve     => from(i)%curve
       nullify(from(i)%curve)
       to(i)%isplit    =  from(i)%isplit
       to(i)%halfedges =  from(i)%halfedges
       to(i)%hyperedge =  from(i)%hyperedge
    end do
    
  end subroutine transfer_edges



  subroutine transfer_faces( &
       from, &
       to, &
       n )
    implicit none
    integer,             intent(in)    :: n
    type(type_BREPface), intent(inout) :: from(:)
    type(type_BREPface), intent(inout) :: to(:)
    integer                            :: i

    do i = 1,n
       to(i)%surface   => from(i)%surface
       nullify(from(i)%surface)
       to(i)%outer     =  from(i)%outer
       call move_alloc(from(i)%inner, to(i)%inner)
       to(i)%ninner    =  from(i)%ninner
       to(i)%hyperface =  from(i)%hyperface
    end do
    
  end subroutine transfer_faces
  



  subroutine free_brep( &
       self )
    implicit none
    type(type_brep), intent(inout) :: self
    integer                        :: ivert, iedge, iface

    do ivert = 1,self%nv
       nullify(self%verts(ivert)%point)
    end do

    do iedge = 1,self%ne
       nullify(self%edges(iedge)%curve)
    end do

    do iface = 1,self%nf
       nullify(self%faces(iface)%surface)
    end do

    if ( allocated(self%verts) ) deallocate(self%verts)
    if ( allocated(self%edges) ) deallocate(self%edges)
    if ( allocated(self%faces) ) deallocate(self%faces)
    
    self%nv = 0
    self%ne = 0
    self%nf = 0
    
  end subroutine free_brep
  
end module mod_brep
