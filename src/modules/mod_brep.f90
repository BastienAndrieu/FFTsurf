module mod_brep

  use mod_math
  use mod_diffgeom
  use mod_types_intersection
  
  implicit none

  integer, parameter :: PARAM_xtra_nv = 20
  integer, parameter :: PARAM_xtra_ne = 20
  integer, parameter :: PARAM_xtra_nf = 20

  type type_BREPface
     type(type_surface), pointer            :: surface => null()
     integer                                :: hyperface = 0
  end type type_BREPface

  type type_BREPedge
     type(type_intersection_curve), pointer :: curve => null()
     integer                                :: isplit = 0
     type(type_DCELhalfedge), pointer       :: halfedge => null()
     integer                                :: hyperedge = 0
     real(kind=fp)                          :: tangents(3,2) ! tangent directions at endpoints in R^3
  end type type_BREPedge
  

  
  type type_DCELface     
     integer                                :: outer = 0 ! DCELhalfedge index
     integer, allocatable                   :: inner(:) ! list of DCELhalfedge indices
     integer                                :: ninner = 0
  end type type_DCELface

  type type_DCELhalfedge
     integer                                :: face = 0 ! DCELface index
     type(type_DCELhalfedge), pointer       :: twin => null()
     type(type_DCELhalfedge), pointer       :: prev => null()
     type(type_DCELhalfedge), pointer       :: next => null()
     !integer                                :: twin = 0 ! DCELhalfedge index
     !integer                                :: prev = 0 ! DCELhalfedge index
     !integer                                :: next = 0 ! DCELhalfedge index
     integer                                :: orig = 0 ! DCELvertex index
     integer                                :: edge = 0 ! BREPedge index
     !type(type_BREPedge), pointer           :: edge => null()
  end type type_DCELhalfedge
  
  type type_DCELvertex
     integer                                :: halfedge = 0 ! outgoing DCELhalfedge index
     ! (if the vertex is on a boundary, the halfedge must be on a boundary as well)
  end type type_DCELvertex

  type type_DCEL
     ! Doubly-Connected Edge List
     integer                                :: nv = 0, ne = 0, nf = 0
     type(type_DCELvertex), allocatable     :: vertices(:)
     type(type_DCELhalfedge), allocatable   :: halfedges(:)
     type(type_DCELface), allocatable       :: faces(:)
  end type type_DCEL
  

  type type_BREP
     ! Boundary Representation
     integer                                :: ne = 0
     type(type_DCEL)                        :: dcel
     type(type_BREPface), allocatable       :: faces(:)
     type(type_BREPedge), allocatable       :: edges(:)
  end type type_BREP


  
contains

  function is_smooth(self)
    ! returns the logical value of a BREPedge's "smooth" tag
    implicit none
    type(type_BREPedge), intent(in) :: self
    logical                         :: is_smooth
    
    is_smooth = self%curve%smooth
    
  end function is_smooth

  
  
  function get_BREPedge_vertices(self) result(verts)
    ! returns the DCELvertex indices of a BREPedge's endpoints
    implicit none
    type(type_BREPedge), intent(in) :: self
    integer                         :: verts(2)
    
    verts(1) = self%halfedge%orig
    verts(2) = self%halfedge%next%orig
    
  end function get_BREPedge_vertices

  


  subroutine brep_from_intersection_data( &
       surf, &
       nsurf, &
       interdata, &
       brep )
    USE MOD_UTIL
    use mod_graph
    use mod_intersection
    use mod_tolerances
    implicit none
    LOGICAL, PARAMETER :: DEBUG = .TRUE.
    integer,                      intent(in)    :: nsurf
    type(type_surface), target,   intent(in)    :: surf(nsurf)
    type(type_intersection_data), intent(inout) :: interdata
    type(type_BREP),              intent(inout) :: brep
    integer, dimension(interdata%nc)            :: arc2curve, arc2split, arc_sens
    integer                                     :: arc2nod(2,interdata%nc)
    real(kind=fp)                               :: arc_angles(2,interdata%nc)
    integer                                     :: point2nod(interdata%np), nod2point(interdata%np)
    real(kind=fp)                               :: nod_uv(2,interdata%np)
    integer                                     :: narc, nnod
    real(kind=fp)                               :: duv_ds(2,2,2), dxyz_ds(3,2)
    integer                                     :: stat
    integer                                     :: nloops
    integer, allocatable                        :: looparc(:,:), loopnod(:,:), lenloop(:)
    integer, allocatable                        :: nascendants(:), ascendants(:,:)
    logical                                     :: inside
    integer, allocatable                        :: outer(:), inner(:,:), ninner(:)
    integer                                     :: nfaces
    integer                                     :: isurf, jsurf, icurv, isplit, jsplit, ipoint, iloop, jloop, i, iface
    real(kind=fp), allocatable                  :: uv_contour(:,:)
    integer                                     :: np_contour
    INTEGER :: FID, IARC, INOD
    CHARACTER(3) :: STRNUM, STRNUM2

    do isurf = 1,nsurf
       IF ( DEBUG ) THEN
          PRINT *,''
          PRINT *,'SURF #',ISURF
       END IF
                 
       ! build graph of intersection curves/points for the current surface
       narc = 0
       nnod = 0
       point2nod(:) = 0
       do icurv = 1,interdata%nc ! <------------------------------------------------------------------------+
          if ( .not.allocated(interdata%curves(icurv)%iedge) ) then ! <----------------+                    !
             allocate(interdata%curves(icurv)%iedge(interdata%curves(icurv)%nsplit))   !                    !
             interdata%curves(icurv)%iedge(:) = 0                                      !                    !
          end if ! <-------------------------------------------------------------------+                    !
          do jsurf = 1,2 ! <----------------------------------------------------------------------------+   !
             if ( associated(interdata%curves(icurv)%surf(jsurf)%ptr, surf(isurf)) ) then ! <-------+   !   !
                !                                                                                   !   !   !
                !IF (ISURF == 3) PRINT *,'ICURV =',ICURV,', JSURF =',JSURF
                do jsplit = 1,interdata%curves(icurv)%nsplit ! <--------------------+               !   !   !
                   if ( jsurf == 1 ) then ! <-----------------------------------+   !               !   !   !
                      isplit = interdata%curves(icurv)%nsplit - jsplit + 1      !   !               !   !   !
                   else ! ------------------------------------------------------+   !               !   !   !
                      isplit = jsplit                                           !   !               !   !   !
                   end if ! <---------------------------------------------------+   !               !   !   !
                   ipoint = interdata%curves(icurv)%isplit(1,isplit)                !               !   !   !
                   !                                                                !               !   !   !
                   if ( point2nod(ipoint) == 0 ) then ! <-----------------------+   !               !   !   !
                      nnod = nnod + 1                                           !   !               !   !   !
                      point2nod(ipoint) = nnod                                  !   !               !   !   !
                      nod2point(nnod) = ipoint                                  !   !               !   !   !
                      nod_uv(1:2,nnod) = interdata%curves(icurv)%polyline%uv&   !   !               !   !   !
                           (1:2,jsurf,interdata%curves(icurv)%isplit(2,isplit)) !   !               !   !   !
                   end if ! <---------------------------------------------------+   !               !   !   !
                end do ! <----------------------------------------------------------+               !   !   !
                !                                                                                   !   !   !
                do jsplit = 1,interdata%curves(icurv)%nsplit-1 ! <------------------------------+   !   !   !
                   if ( jsurf == 1 ) then ! <-----------------------------------+               !   !   !   !
                      isplit = interdata%curves(icurv)%nsplit - jsplit + 1      !               !   !   !   !
                   else ! ------------------------------------------------------+               !   !   !   !
                      isplit = jsplit                                           !               !   !   !   !
                   end if ! <---------------------------------------------------+               !   !   !   !
                   narc = narc + 1                                                              !   !   !   !
                   arc2curve(narc) = icurv                                                      !   !   !   !
                   arc2split(narc) = isplit                                                     !   !   !   !
                   arc_sens(narc) = jsurf
                   IF (ISURF == 3) PRINT *,'arc_sens(narc) =',arc_sens(narc)
                   IF ( .FALSE. ) THEN
                      if ( jsplit == 1 ) then ! <----------------------------------------------+   !   !   !   !
                         arc2nod(1,narc) = point2nod(interdata%curves(icurv)%isplit(1,isplit)) !   !   !   !   !
                         ipoint = interdata%curves(icurv)%isplit(2,isplit)                     !   !   !   !   !
                         call diffgeom_intersection( &                                         !   !   !   !   !
                              interdata%curves(icurv)%surf, &                                  !   !   !   !   !
                              interdata%curves(icurv)%polyline%uv(:,:,ipoint), &               !   !   !   !   !
                              duv_ds, &                                                        !   !   !   !   !
                              dxyz_ds, &                                                       !   !   !   !   !
                              stat )                                                           !   !   !   !   !
                         duv_ds = real((-1)**jsurf, kind=fp) * duv_ds                          !   !   !   !   !
                         arc_angles(1,narc) = atan2(duv_ds(2,1,jsurf), duv_ds(1,1,jsurf))      !   !   !   !   !
                      else ! ------------------------------------------------------------------+   !   !   !   !
                         arc2nod(1,narc) = arc2nod(2,narc-1)                                   !   !   !   !   !
                         arc_angles(1,narc) = arc_angles(2,narc-1)                             !   !   !   !   !
                      end if ! <---------------------------------------------------------------+   !   !   !   !
                      isplit = isplit + (-1)**jsurf                                                !   !   !   !
                      arc2nod(2,narc) = point2nod(interdata%curves(icurv)%isplit(1,isplit))        !   !   !   !
                      ipoint = interdata%curves(icurv)%isplit(2,isplit)                            !   !   !   !
                      call diffgeom_intersection( &                                                !   !   !   !
                           interdata%curves(icurv)%surf, &                                         !   !   !   !
                           interdata%curves(icurv)%polyline%uv(:,:,ipoint), &                      !   !   !   !
                           duv_ds, &                                                               !   !   !   !
                           dxyz_ds, &                                                              !   !   !   !
                           stat )                                                                  !   !   !   !
                      duv_ds = real((-1)**jsurf, kind=fp) * duv_ds                                 !   !   !   !
                      arc_angles(2,narc) = atan2(duv_ds(2,1,jsurf), duv_ds(1,1,jsurf))             !   !   !   !
                   ELSE
                      do i = 1,2
                         arc2nod(i,narc) = point2nod(interdata%curves(icurv)%isplit(1,isplit))
                         ipoint = interdata%curves(icurv)%isplit(2,isplit)
                         IF ( .TRUE. ) THEN
                            duv_ds(:,1,jsurf) = &
                                 interdata%curves(icurv)%polyline%uv(:,jsurf,ipoint + (2-i)*(-1)**jsurf) - &
                                 interdata%curves(icurv)%polyline%uv(:,jsurf,ipoint + (1-i)*(-1)**jsurf)
                         ELSE
                            call diffgeom_intersection( &
                                 interdata%curves(icurv)%surf, &
                                 interdata%curves(icurv)%polyline%uv(:,:,ipoint), &
                                 duv_ds, &
                                 dxyz_ds, &
                                 stat )
                            IF ( STAT > 1 ) PRINT *,'STAT =',STAT
                            duv_ds = real((-1)**jsurf, kind=fp) * duv_ds
                         END IF
                         arc_angles(i,narc) = atan2(duv_ds(2,1,jsurf), duv_ds(1,1,jsurf))
                         isplit = isplit + (-1)**jsurf
                      end do
                   END IF
                end do ! <----------------------------------------------------------------------+   !   !   !
                !                                                                                   !   !   !
                exit                                                                                !   !   !
             end if ! <-----------------------------------------------------------------------------+   !   !
          end do ! <------------------------------------------------------------------------------------+   !
       end do ! <-------------------------------------------------------------------------------------------+

       IF ( DEBUG ) THEN
          PRINT *,'NARC =',NARC
          PRINT *,'NNOD =',NNOD
          WRITE (STRNUM,'(I3.3)') ISURF
          CALL GET_FREE_UNIT(FID)
          OPEN(UNIT=FID, FILE='Jouke/graph/graph_'//strnum//'.dat', ACTION='WRITE')
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
       !
       ! make loops
       nloops = min(nnod, narc)
       allocate(lenloop(nloops), looparc(nloops,nloops), loopnod(nloops,nloops))
       call make_loops( &
            arc2nod(1:2,1:narc), &
            arc_angles(1:2,1:narc), &
            narc, &
            nnod, &
            nloops, &
            looparc, &
            !loopnod, &
            lenloop )
       IF ( DEBUG ) THEN
          PRINT *,'NLOOPS =',NLOOPS
          WRITE (STRNUM,'(I3.3)') ISURF
          CALL GET_FREE_UNIT(FID)
          OPEN(UNIT=FID, FILE='Jouke/graph/loops_'//strnum//'.dat', ACTION='WRITE')
          WRITE (FID,*) NLOOPS
          DO ILOOP = 1,NLOOPS
             WRITE (FID,*) LOOPARC(1:LENLOOP(ILOOP),ILOOP)
             WRITE (FID,*) LOOPNOD(1:LENLOOP(ILOOP),ILOOP)
          END DO
          CLOSE(FID)
       END IF
       !
       if ( nloops > 0 ) then ! <------------------------------------------------------+
          ! get inclusion relationships between nested loops                           !
          allocate(nascendants(nloops), ascendants(nloops,nloops))                     !
          nascendants(:) = 0                                                           !
          do iloop = 1,nloops-1 ! <---------------------------------------------+      !
             do jloop = iloop+1,nloops ! <----------------------------------+   !      !
                ! get a point on loop #jloop, and test whether it is        !   !      !
                ! inside the polygon of loop #iloop                         !   !      !
                call point_in_loop( &                                       !   !      !
                     nod_uv(:,arc2nod(1,looparc(1,jloop))), &               !   !      !
                     interdata, &                                           !   !      !
                     lenloop(iloop), &                                      !   !      !
                     arc2curve(looparc(1:lenloop(iloop),iloop)), &          !   !      !
                     arc2split(looparc(1:lenloop(iloop),iloop)), &          !   !      !
                     arc_sens(looparc(1:lenloop(iloop),iloop)), &           !   !      !
                     inside )                                               !   !      !
                if ( inside ) then ! <----------------------------------+   !   !      !
                   ! loop #jloop is nested inside loop #iloop           !   !   !      !
                   nascendants(jloop) = nascendants(jloop) + 1          !   !   !      !
                   ascendants(nascendants(jloop),jloop) = iloop         !   !   !      !
                   cycle                                                !   !   !      !
                else ! -------------------------------------------------+   !   !      !
                   ! get a point on loop #iloop, and test whether it is !   !   !      !
                   ! inside the polygon of loop #jloop                  !   !   !      !
                   call point_in_loop( &                                !   !   !      !
                        nod_uv(:,arc2nod(1,looparc(1,iloop))), &        !   !   !      !
                        interdata, &                                    !   !   !      !
                        lenloop(jloop), &                               !   !   !      !
                        arc2curve(looparc(1:lenloop(jloop),jloop)), &   !   !   !      !
                        arc2split(looparc(1:lenloop(jloop),jloop)), &   !   !   !      !
                        arc_sens(looparc(1:lenloop(jloop),jloop)), &    !   !   !      !
                        inside )                                        !   !   !      !
                   if ( inside ) then ! <---------------------------+   !   !   !      !
                      ! loop #iloop is nested inside loop #jloop    !   !   !   !      !
                      nascendants(iloop) = nascendants(iloop) + 1   !   !   !   !      !
                      ascendants(nascendants(iloop),iloop) = jloop  !   !   !   !      !
                   end if ! <---------------------------------------+   !   !   !      !
                end if ! <----------------------------------------------+   !   !      !
             end do ! <-----------------------------------------------------+   !      !
          end do ! <------------------------------------------------------------+      !
          !                                                                            !
          ! make faces                                                                 !
          allocate(outer(nloops), inner(nloops,nloops), ninner(nloops))                !
          call make_faces( &                                                           !
               nascendants(1:nloops), &                                                !
               ascendants(1:nloops,1:nloops), &                                        !
               nloops, &                                                               !
               outer, &                                                                !
               inner, &                                                                !
               ninner, &                                                               !
               nfaces )                                                                !
          IF ( DEBUG ) THEN
             PRINT *,NFACES,' FACE(S) :'
             DO IFACE = 1,NFACES
                PRINT *,'   ',NINNER(IFACE),' HOLE(S)'
             END DO

             WRITE (STRNUM,'(I3.3)') ISURF
             CALL GET_FREE_UNIT(FID)
             OPEN(UNIT=FID, FILE='Jouke/graph/faces_'//strnum//'.dat', ACTION='WRITE')
             WRITE (FID,*) NFACES
             DO IFACE = 1,NFACES
                WRITE (FID,*) LOOPARC(1:LENLOOP(OUTER(IFACE)),OUTER(IFACE))
                WRITE (FID,*) NINNER(IFACE)
                DO ILOOP = 1,NINNER(IFACE)
                   WRITE (FID,*) LOOPARC(1:LENLOOP(INNER(ILOOP,IFACE)),INNER(ILOOP,IFACE))
                END DO
             END DO
             CLOSE(FID)
          END IF
          !
          do iface = 1,nfaces
             brep%dcel%nf = brep%dcel%nf + 1
             brep%faces(brep%dcel%nf)%surface => surf(isurf)

             call insert_loop_in_dcel( &
                  interdata, &
                  brep, &
                  brep%dcel%nf, &
                  .true., &
                  looparc(1:lenloop(outer(iface)),outer(iface)), &
                  lenloop(outer(iface)), &
                  arc2curve, &
                  arc2split, &
                  arc2nod, &
                  narc, &
                  nod2point, &
                  nnod )

             IF ( DEBUG ) THEN
                WRITE (STRNUM,'(I3.3)') brep%dcel%nf
                CALL GET_FREE_UNIT(FID)
                OPEN(UNIT=FID, FILE='Jouke/meshgen/numsurf_'//strnum//'.dat', ACTION='WRITE')
                WRITE (FID,*) ISURF
                CLOSE(FID)

                call get_contour_uv( &
                     interdata, &
                     lenloop(outer(iface)), &
                     arc2curve(looparc(1:lenloop(outer(iface)),outer(iface))), &
                     arc2split(looparc(1:lenloop(outer(iface)),outer(iface))), &
                     arc_sens(looparc(1:lenloop(outer(iface)),outer(iface))), &
                     uv_contour, &
                     np_contour )
                OPEN(UNIT=FID, FILE='Jouke/meshgen/contours/uv_'//strnum//'.dat', ACTION='WRITE')
                WRITE (FID,*) np_contour, 2
                do ipoint = 1,np_contour
                   write (fid,*) uv_contour(:,ipoint)
                end do
                CLOSE(FID)

                OPEN(UNIT=FID, FILE='Jouke/meshgen/contours/edges_'//strnum//'.dat', ACTION='WRITE')
                WRITE (FID,*) np_contour
                do ipoint = 1,np_contour
                   write (fid,*) ipoint, 1+mod(ipoint,np_contour)
                end do
                CLOSE(FID)

                if ( allocated(uv_contour) ) deallocate(uv_contour)

                WRITE (STRNUM2,'(I3.3)') isurf
                call system('/stck/bandrieu/Bureau/MeshGen/./meshgen.out &
                     &Jouke/propergol/C_'//strnum2//'.cheb &
                     &Jouke/meshgen/contours/uv_'//strnum//'.dat &
                     &Jouke/meshgen/contours/edges_'//strnum//'.dat &
                     &Jouke/meshgen/info.dat &
                     &Jouke/meshgen/tri_'//strnum//'.dat &
                     &Jouke/meshgen/uv_'//strnum//'.dat')
             END IF
             
          end do
          !                                                                            !
          deallocate(nascendants, ascendants, outer, inner, ninner)                    !
       end if ! <----------------------------------------------------------------------+
       !
       deallocate(lenloop, looparc, loopnod)
       !
    end do

    IF ( DEBUG ) PRINT *,BREP%DCEL%NF,' FACES'

    !IF ( DEBUG ) THEN
    !   OPEN(UNIT=FID, FILE='Jouke/meshgen/info.dat', ACTION='write')
    !   WRITE (FID,*) 5.D-4    ! LMIN
    !   WRITE (FID,*) 1.D-1    ! LMAX
    !   WRITE (FID,*) TOLchord ! TOL
    !   CLOSE(FID)
    !END IF
  end subroutine brep_from_intersection_data





  subroutine insert_loop_in_dcel( &
       interdata, &
       brep, &
       iface, &
       is_outer_loop, &
       looparc, &
       lenloop, &
       arc2curve, &
       arc2split, &
       arc2nod, &
       narc, &
       nod2point, &
       nnod )
    use mod_util
    implicit none
    type(type_intersection_data), intent(inout), target :: interdata
    type(type_BREP),              intent(inout), target :: brep
    integer,                      intent(in)            :: iface
    logical,                      intent(in)            :: is_outer_loop
    integer,                      intent(in)            :: lenloop
    integer,                      intent(in)            :: looparc(lenloop)
    integer,                      intent(in)            :: narc
    integer,                      intent(in)            :: arc2curve(narc)
    integer,                      intent(in)            :: arc2split(narc)
    integer,                      intent(in)            :: arc2nod(2,narc)
    integer,                      intent(in)            :: nnod
    integer,                      intent(in)            :: nod2point(nnod)
    integer                                             :: nheprev
    integer                                             :: ihalfedge, iarc, jarc, ipoint, iprev, inext, iedge

    nheprev = brep%dcel%ne
    if ( .not.allocated(brep%dcel%halfedges) .or. &
         nheprev + lenloop > size(brep%dcel%halfedges) ) then ! <----+
       call reallocate_dcel_halfedges( &                             !
            brep%dcel, &                                             !
            nheprev + lenloop + 2*PARAM_xtra_ne )                    !
    end if ! <-------------------------------------------------------+
    
    do jarc = 1,lenloop ! <---------------------------------------------------------+
       ihalfedge = nheprev + jarc                                                   !
       if ( jarc == 1 ) then ! <------------------------------------+               !
          if ( is_outer_loop ) then ! <-------------------------+   !               !
             brep%dcel%faces(iface)%outer = ihalfedge           !   !               !
          else ! -----------------------------------------------+   !               !
             call insert_after( &                               !   !               !
                  brep%dcel%faces(iface)%inner, &               !   !               !
                  brep%dcel%faces(iface)%ninner, &              !   !               !
                  ihalfedge, &                                  !   !               !
                  brep%dcel%faces(iface)%ninner )               !   !               !
          end if ! <--------------------------------------------+   !               !
       end if ! <---------------------------------------------------+               !
       iarc = looparc(jarc)                                                         !
       !                                                                            !
       brep%dcel%ne = brep%dcel%ne + 1                                              !
       brep%dcel%halfedges(ihalfedge)%face = brep%dcel%nf                           !
       iprev = nheprev + 1 + mod(jarc + lenloop - 2,lenloop)                        !
       inext = nheprev + 1 + mod(jarc,lenloop)                                      !
       brep%dcel%halfedges(ihalfedge)%prev => brep%dcel%halfedges(iprev)            !
       brep%dcel%halfedges(ihalfedge)%next => brep%dcel%halfedges(inext)            !
       !                                                                            !
       ipoint = nod2point(arc2nod(1,iarc))                                          !
       if ( interdata%points(ipoint)%ivert == 0 ) then ! <----+                     !
          ! +1 DECLvertex                                     !                     !
          brep%dcel%nv = brep%dcel%nv + 1                     !                     !
          interdata%points(ipoint)%ivert = brep%dcel%nv       !                     !
       end if ! <---------------------------------------------+                     !
       brep%dcel%halfedges(ihalfedge)%orig = interdata%points(ipoint)%ivert         !
       !                                                                            !
       iedge = interdata%curves(arc2curve(iarc))%iedge(arc2split(iarc))             !
       if ( iedge == 0 ) then ! <-----------------------------------------------+   !
          ! +1 BREPedge                                                         !   !
          brep%ne = brep%ne + 1                                                 !   !
          iedge = brep%ne                                                       !   !
          brep%edges(iedge)%curve => interdata%curves(arc2curve(iarc))          !   !
          brep%edges(iedge)%isplit = arc2split(iarc)                            !   !
          interdata%curves(arc2curve(iarc))%iedge(arc2split(iarc)) = iedge      !   !
          !                                                                     !   !
          brep%edges(brep%ne)%halfedge => brep%dcel%halfedges(ihalfedge)        !   !
          brep%dcel%halfedges(ihalfedge)%edge = iedge                           !   !
       else ! ------------------------------------------------------------------+   !
          brep%dcel%halfedges(ihalfedge)%edge = iedge                           !   !
          brep%dcel%halfedges(ihalfedge)%twin => brep%edges(iedge)%halfedge     !   !
          brep%edges(iedge)%halfedge%twin => brep%dcel%halfedges(ihalfedge)     !   !
       end if ! <---------------------------------------------------------------+   !
    end do ! <----------------------------------------------------------------------+
    
  end subroutine insert_loop_in_dcel





  subroutine reallocate_dcel_halfedges( &
       dcel, &
       nhe )
    implicit none
    integer,               intent(in)    :: nhe
    type(type_DCEL),       intent(inout) :: dcel
    type(type_DCELhalfedge), allocatable :: tmp(:)
    integer                              :: ie

    if ( .not.allocated(dcel%halfedges) ) then
       allocate(dcel%halfedges(nhe))
       return
    end if

    allocate(tmp(dcel%ne))
    do ie = 1,dcel%ne
       tmp(ie)%face =  dcel%halfedges(ie)%face
       tmp(ie)%twin => dcel%halfedges(ie)%twin
       tmp(ie)%prev => dcel%halfedges(ie)%prev
       tmp(ie)%next => dcel%halfedges(ie)%next
       tmp(ie)%orig =  dcel%halfedges(ie)%orig
       tmp(ie)%edge =  dcel%halfedges(ie)%edge
    end do
    deallocate(dcel%halfedges)
    
    allocate(dcel%halfedges(nhe))
    do ie = 1,dcel%ne
       dcel%halfedges(ie)%face =  tmp(ie)%face
       dcel%halfedges(ie)%twin => tmp(ie)%twin
       dcel%halfedges(ie)%prev => tmp(ie)%prev
       dcel%halfedges(ie)%next => tmp(ie)%next
       dcel%halfedges(ie)%orig =  tmp(ie)%orig
       dcel%halfedges(ie)%edge =  tmp(ie)%edge
    end do
    deallocate(tmp)
    
  end subroutine reallocate_dcel_halfedges

  
  

  subroutine get_contour_uv( &
       interdata, &
       lenloop, &
       arc2curve, &
       arc2split, &
       arc_sens, &
       uv, &
       np )
    implicit none
    integer, parameter                          :: PARAM_xtra_np = 20
    type(type_intersection_data), intent(in)    :: interdata
    integer,                      intent(in)    :: lenloop
    integer,                      intent(in)    :: arc2curve(lenloop)
    integer,                      intent(in)    :: arc2split(lenloop)
    integer,                      intent(in)    :: arc_sens(lenloop)
    real(kind=fp), allocatable,   intent(inout) :: uv(:,:)
    integer,                      intent(out)   :: np
    real(kind=fp), allocatable                  :: tmp(:,:)
    integer                                     :: npi
    integer                                     :: iarc, icurv, isplit, ifirst, ilast
    
    np = 0
    if ( .not.allocated(uv) ) allocate(uv(2,PARAM_xtra_np))
    
    do iarc = 1,lenloop
       icurv = arc2curve(iarc)
       isplit = arc2split(iarc)
       ifirst = interdata%curves(icurv)%isplit(2,isplit)
       ilast  = interdata%curves(icurv)%isplit(2,isplit + (-1)**arc_sens(iarc)) - (-1)**arc_sens(iarc)
       npi = (-1)**arc_sens(iarc) * (ilast - ifirst) + 1

       if ( np + npi > size(uv,2) ) then
          call move_alloc(from=uv, to=tmp)
          !allocate(uv(2,np + npi + PARAM_xtra_np), source=tmp)
          allocate(uv(2,np + npi + PARAM_xtra_np))
          uv(1:2,1:size(tmp,2)) = tmp(1:2,1:size(tmp,2))
          deallocate(tmp)
       end if

       uv(1:2,np+1:np+npi) = interdata%curves(icurv)%polyline%uv&
            (1:2,arc_sens(iarc),ifirst:ilast:(-1)**arc_sens(iarc))
       np = np + npi
    end do
    
  end subroutine get_contour_uv

  
  

  subroutine point_in_loop( &
       uv, &
       interdata, &
       lenloop, &
       arc2curve, &
       arc2split, &
       arc_sens, &
       inside )
    implicit none
    real(kind=fp), parameter                  :: M = 10._fp
    real(kind=fp),                intent(in)  :: uv(2)
    type(type_intersection_data), intent(in)  :: interdata
    integer,                      intent(in)  :: lenloop
    integer,                      intent(in)  :: arc2curve(lenloop)
    integer,                      intent(in)  :: arc2split(lenloop)
    integer,                      intent(in)  :: arc_sens(lenloop)
    logical,                      intent(out) :: inside
    real(kind=fp), pointer                    :: poly(:,:) => null()
    real(kind=fp)                             :: a, c, s, ub, vb
    integer                                   :: np
    integer                                   :: iarc, icurv, isplit, ifirst, ilast, i

    inside = .false.
    
    call random_number(a)
    a = 0.5_fp * CSTpi * a
    c = cos(a)
    s = sin(a)
    ub = uv(1) + M*c
    vb = uv(2) + M*s
    
    do iarc = 1,lenloop
       icurv = arc2curve(iarc)
       isplit = arc2split(iarc)

       ifirst = interdata%curves(icurv)%isplit(2,isplit)
       ilast  = interdata%curves(icurv)%isplit(2,isplit + (-1)**arc_sens(iarc))
       np = (-1)**arc_sens(iarc) * (ilast - ifirst) + 1

       poly => interdata%curves(icurv)%polyline%uv(:,arc_sens(iarc),ifirst:ilast:(-1)**arc_sens(iarc))
       do i = 1,np-1
          if ( ( ((poly(1,i) - uv(1))*s > (poly(2,i) - uv(2))*c) .neqv. &
               ((poly(1,i+1) - uv(1))*s > (poly(2,i+1) - uv(2))*c) ) .and. &
               ( ((uv(1) - poly(1,i))*(poly(2,i) - poly(2,i+1)) > (uv(2) - poly(2,i))*(poly(1,i) - poly(1,i+1))) .neqv. &
               ((ub - poly(1,i))*(poly(2,i) - poly(2,i+1)) > (vb - poly(2,i))*(poly(1,i) - poly(1,i+1))) ) ) then
             inside = .not.inside
          end if
       end do

       nullify(poly)
    end do
  end subroutine point_in_loop
  
  

  

  ! DCEL...

  
  subroutine get_v2f(dcel, iv, faces, nfaces)
    use mod_util
    ! returns the list of all the faces incident to a vertex, traversed in CCW order
    implicit none
    type(type_DCEL),      intent(in), target :: dcel
    integer,              intent(in)         :: iv
    integer, allocatable, intent(inout)      :: faces(:)
    integer,              intent(out)        :: nfaces
    type(type_DCELhalfedge), pointer         :: e => null()
    
    nfaces = 0
    e => DCEL%halfedges(DCEL%vertices(iv)%halfedge) ! outgoing
    do
       if ( nfaces > 0 ) then
          if ( e%face == faces(1) ) return
       end if
       call insert_after( &
            faces, &
            nfaces, &
            e%face, &
            nfaces )
       e => e%prev ! ingoing, same face
       if ( .not.associated(e%twin) ) return
       e => e%twin ! outgoing, different face
    end do
    
  end subroutine get_v2f


  subroutine get_v2e(dcel, iv, edges, nedges)
    use mod_util
    ! returns the list of all the faces incident to a vertex, traversed in CCW order
    implicit none
    type(type_DCEL),      intent(in), target :: dcel
    integer,              intent(in)         :: iv
    integer, allocatable, intent(inout)      :: edges(:)
    integer,              intent(out)        :: nedges
    type(type_DCELhalfedge), pointer         :: e => null()
    
    nedges = 0
    e => DCEL%halfedges(DCEL%vertices(iv)%halfedge) ! outgoing
    do
       call insert_after( &
            edges, &
            nedges, &
            e%edge, &
            nedges )
       e => e%prev ! ingoing, same face
       if ( .not.associated(e%twin) ) then
          call insert_after( &
               edges, &
               nedges, &
               e%edge, &
               nedges )
          return
       end if
       e => e%twin ! outgoing, different face
       if ( e%edge == edges(1) ) return
    end do
    
  end subroutine get_v2e
  
end module mod_brep
