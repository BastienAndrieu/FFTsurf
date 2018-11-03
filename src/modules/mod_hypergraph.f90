module mod_hypergraph

  implicit none

  integer, parameter :: PARAM_xtra_nhf = 10
  integer, parameter :: PARAM_xtra_nhe = 10


  type type_hyperface
     integer              :: nf = 0
     integer, allocatable :: faces(:)
  end type type_hyperface

  type type_hyperedge
     integer              :: ne = 0
     integer, allocatable :: halfedges(:,:)
     integer              :: verts(2) = 0
     integer              :: hyperfaces(2) = 0
  end type type_hyperedge

  type type_hypergraph
     integer                           :: nhe = 0, nhf = 0
     type(type_hyperface), allocatable :: hyperfaces(:)
     type(type_hyperedge), allocatable :: hyperedges(:)     
  end type type_hypergraph

contains


  subroutine make_hypergraph( &
       brep, &
       hypergraph, &
       feat_edge, &
       feat_vert )
    use mod_types_brep
    implicit none
    type(type_brep),       intent(inout) :: brep
    type(type_hypergraph), intent(inout) :: hypergraph
    logical, allocatable,  intent(inout) :: feat_edge(:)
    logical, allocatable,  intent(inout) :: feat_vert(:)
    integer                              :: valence(brep%nv)

    ! get feature edges
    if ( allocated(feat_edge) ) then
       if ( size(feat_edge) < brep%ne ) deallocate(feat_edge)
    end if
    if ( .not.allocated(feat_edge) ) allocate(feat_edge(brep%ne))
    call get_feature_edges( &
         brep, &
         feat_edge(1:brep%ne) )

    ! get hyperfaces
    hypergraph%nhf = 0
    call get_hyperfaces( &
         brep, &
         feat_edge(1:brep%ne), &
         hypergraph%hyperfaces, &
         hypergraph%nhf )

    ! get feature vertices
    if ( allocated(feat_vert) ) then
       if ( size(feat_vert) < brep%nv ) deallocate(feat_vert)
    end if
    if ( .not.allocated(feat_vert) ) allocate(feat_vert(brep%nv))
    call get_feature_vertices( &
         brep, &
         feat_edge(1:brep%ne), &
         feat_vert(1:brep%nv), &
         valence )

    ! get hyperedges
    hypergraph%nhe = 0
    call get_hyperedges( &
         brep, &
         feat_edge(1:brep%ne), &
         feat_vert(1:brep%nv), &
         valence, &
         hypergraph%hyperedges, &
         hypergraph%nhe )

  end subroutine make_hypergraph




  subroutine get_feature_edges( &
       brep, &
       feat_edge )
    use mod_types_brep
    use mod_brep
    ! returns the list of feature edges in a BREP, i.e. non smooth or boundary edges
    implicit none
    type(type_BREP), intent(in)  :: brep
    logical,         intent(out) :: feat_edge(brep%ne)
    integer                      :: ie

    do ie = 1,brep%ne
       feat_edge(ie) = ( .not.is_smooth(brep, ie) .or. is_boundary_edge(brep, ie) )
    end do

  end subroutine get_feature_edges


  subroutine get_feature_vertices( &
       brep, &
       feat_edge, &
       feat_vert, &
       valence )
    use mod_intersection
    use mod_types_brep
    use mod_brep
    implicit none
    real(kind=fp), parameter     :: EPScosangle = real(1.d-3, kind=fp)
    type(type_BREP), intent(in)  :: brep
    logical,         intent(in)  :: feat_edge(brep%ne)
    logical,         intent(out) :: feat_vert(brep%nv)
    integer,         intent(out) :: valence(brep%nv)
    real(kind=fp)                :: tangents(3,2,brep%ne)
    real(kind=fp)                :: duv_ds(2,2,2), dxyz_ds(3,2)
    integer                      :: stat
    real(kind=fp)                :: vec(3,2), cosangle
    integer                      :: ivert, iedge(2), jedge, iend, isplit, ipoint

    feat_vert(:) = .false.
    valence(:) = 0

    ! compute unit tangent vectors at feature edges' endpoints
    do jedge = 1,brep%ne ! <----------------------------------------------+
       if ( .not.feat_edge(jedge) ) cycle                                 !
       !                                                                  !
       do iend = 1,2 ! <---------------------------------------------+    !
          isplit = brep%edges(jedge)%isplit + iend - 1               !    !
          ! get polyline point index                                 !    !
          ipoint = brep%edges(jedge)%curve%isplit(2,isplit)          !    !
          ! compute tangent to intersection curve                    !    !
          ! ***** GERER CAS BOUNDARY EDGE! *****
          call diffgeom_intersection( &                              !    !
               brep%edges(jedge)%curve%surf, &                       !    !
               brep%edges(jedge)%curve%polyline%uv(:,:,ipoint), &    !    !
               duv_ds, &                                             !    !
               dxyz_ds, &                                            !    !
               stat )                                                !    !
          if ( stat /= 0 ) STOP 'get_feature_vertices : &            
               &feature edge is not a transversal intersection'      !    !
          tangents(1:3,iend,jedge) = dxyz_ds(1:3,1)                  !    !
       end do ! <----------------------------------------------------+    !
    end do ! <------------------------------------------------------------+

    do ivert = 1,brep%nv ! <-----------------------------------------------------+
       iedge = brep%verts(ivert)%halfedge ! outgoing                             !
       do ! <----------------------------------------------------------------+   !
          if ( feat_edge(iedge(1)) ) then ! <-----------------------------+  !   !
             valence(ivert) = valence(ivert) + 1                          !  !   !
             if ( valence(ivert) <= 2 ) then ! <-----------------------+  !  !   !
                vec(:,valence(ivert)) = tangents(:,iedge(2),iedge(1))  !  !  !   !
                if ( mod(valence(ivert) + iedge(2),2) == 0 ) &         !  !  !   !
                     vec(:,valence(ivert)) = -vec(:,valence(ivert))    !  !  !   !
             end if ! <------------------------------------------------+  !  !   !
          end if ! <------------------------------------------------------+  !   !
          !                                                                  !   !
          iedge = get_prev(brep, iedge) ! ingoing                            !   !
          iedge = get_twin(iedge)       ! outgoing                           !   !
          if ( brep%edges(iedge(1))%halfedges(iedge(2))%face < 1 .or. &      !   !
               iedge(1) == brep%verts(ivert)%halfedge(1) ) exit              !   !
       end do ! <------------------------------------------------------------+   !
       !                                                                         !
       if ( valence(ivert) == 0 ) cycle                                          !
       !                                                                         !
       if ( valence(ivert) == 2 ) then ! <-----------------------------------+   !
          ! check angle between left and right tangent directions            !   !
          cosangle = dot_product(vec(:,1), vec(:,2))                         !   !
          if ( abs(cosangle) > 1._fp - EPScosangle ) cycle                   !   !
       end if ! <------------------------------------------------------------+   !
       !                                                                         !
       feat_vert(ivert) = .true.                                                 !
    end do ! <-------------------------------------------------------------------+

  end subroutine get_feature_vertices



  subroutine get_hyperfaces( &
       brep, &
       feat_edge, &
       hyperfaces, &
       nhf )
    use mod_types_brep
    use mod_brep
    implicit none
    type(type_BREP),                   intent(inout) :: brep
    logical,                           intent(in)    :: feat_edge(brep%ne)
    type(type_hyperface), allocatable, intent(inout) :: hyperfaces(:)
    integer,                           intent(out)   :: nhf
    type(type_hyperface), allocatable                :: tmp(:)
    logical                                          :: visited(brep%nf)
    integer                                          :: ifa, ihf

    visited(:) = .false.
    nhf = 0
    do while ( any(.not.visited) )
       do ifa = 1,brep%nf
          ! look for an unvisited face to start a new hyperface
          if ( .not.visited(ifa) ) then

             if ( .not.allocated(hyperfaces) ) allocate(hyperfaces(PARAM_xtra_nhf))

             if ( nhf + 1 > size(hyperfaces) ) then
                ! reallocate hyperfaces with bigger size
                allocate(tmp(nhf))
                call transfer_hyperfaces( &
                     hyperfaces, &
                     tmp, &
                     nhf )
                deallocate(hyperfaces)

                allocate(hyperfaces(nhf + PARAM_xtra_nhf))
                call transfer_hyperfaces( &
                     tmp, &
                     hyperfaces, &
                     nhf )
                deallocate(tmp)
             end if

             nhf = nhf + 1
             hyperfaces(nhf)%nf = 0

             call build_hyperface( &
                  brep, &
                  ifa, &
                  feat_edge, &
                  hyperfaces(nhf), &
                  visited )
             !print *,''
          end if
       end do
    end do

    do ihf = 1,nhf
       do ifa = 1,hyperfaces(ihf)%nf
          brep%faces(hyperfaces(ihf)%faces(ifa))%hyperface = ihf
       end do
    end do

  end subroutine get_hyperfaces




  recursive subroutine build_hyperface( &
       brep, &
       ifa, &
       feat_edge, &
       hf, &
       visited )
    use mod_util
    use mod_types_brep
    use mod_brep
    implicit none
    type(type_BREP),      intent(in)    :: brep
    integer,              intent(in)    :: ifa
    logical,              intent(in)    :: feat_edge(brep%ne)
    type(type_hyperface), intent(inout) :: hf
    logical,              intent(inout) :: visited(brep%nf)
    integer                             :: iedge(2)
    integer                             :: iwire, istart, jfa

    if ( visited(ifa) ) return

    call insert_after( &
         hf%faces, &
         hf%nf, &
         ifa, &
         hf%nf )

    visited(ifa) = .true.

    do iwire = 1,1 + brep%faces(ifa)%ninner ! <---------------------------------------+
       if ( iwire == 1 ) then ! <--------------------------------------+              !
          ! traverse outer wire                                        !              !
          iedge = brep%faces(ifa)%outer                                !              !
       else ! ---------------------------------------------------------+              !
          ! traverse potential inner wires                             !              !
          iedge = brep%faces(ifa)%inner(:,iwire-1)                     !              !
       end if ! <------------------------------------------------------+              !
       istart = iedge(1)                                                              !
       do ! <-------------------------------------------------------------------+     !
          if ( .not.feat_edge(iedge(1)) ) then ! <--------------------------+   !     !
             jfa = brep%edges(iedge(1))%halfedges(1+mod(iedge(2),2))%face   !   !     !
             call build_hyperface(brep, jfa, feat_edge, hf, visited)        !   !     !
          end if ! <--------------------------------------------------------+   !     !
          iedge = brep%edges(iedge(1))%halfedges(iedge(2))%next                 !     !
          if ( iedge(1) == istart ) exit                                        !     !
       end do ! <---------------------------------------------------------------+     !
    end do ! <------------------------------------------------------------------------+

  end subroutine build_hyperface





  subroutine get_hyperedges( &
       brep, &
       feat_edge, &
       feat_vert, &
       valence, &
       hyperedges, &
       nhe )
    use mod_types_brep
    use mod_brep
    implicit none
    type(type_BREP),                   intent(inout) :: brep
    logical,                           intent(in)    :: feat_edge(brep%ne)
    logical,                           intent(in)    :: feat_vert(brep%nv)
    integer,                           intent(inout) :: valence(brep%nv)
    type(type_hyperedge), allocatable, intent(inout) :: hyperedges(:)
    integer,                           intent(out)   :: nhe
    logical                                          :: visited(brep%ne)
    integer                                          :: ivert, iedge, ihe, ihedg(2), ied

    visited(:) = .false.
    nhe = 0
    ! first, build open hyperedges (paths)
    do ivert = 1,brep%nv ! <-------------------------------------------------------------------------+
       if ( .not.feat_vert(ivert) ) cycle                                                            !
       do while ( valence(ivert) > 0 )                                                               !
          !                                                                                          !
          if ( .not.allocated(hyperedges) .or. &                                                     !
               nhe + 1 > size(hyperedges) ) then ! <-------+                                         !
             call reallocate_hyperedges( &                 !                                         !
                  hyperedges, &                            !                                         !
                  nhe + PARAM_xtra_nhe )                   !                                         !
          end if ! <---------------------------------------+                                         !
          !                                                                                          !
          nhe = nhe + 1                                                                              !
          PRINT *,'+1 PATH'                                                                          !
          hyperedges(nhe)%ne = 0                                                                     !
          call build_hyperedge( &                                                                    !
               brep, &                                                                               !
               ivert, &                                                                              !
               feat_edge, &                                                                          !
               visited, &                                                                            !
               feat_vert, &                                                                          !
               valence, &                                                                            !
               hyperedges(nhe) )                                                                     !
          !                                                                                          !
          ! set boundary vertices indices                                                            !
          hyperedges(nhe)%verts(1) = get_orig(brep, hyperedges(nhe)%halfedges(:,1))                  !
          hyperedges(nhe)%verts(2) = get_dest(brep, hyperedges(nhe)%halfedges(:,hyperedges(nhe)%ne)) !
          PRINT *,'   VERTS =',hyperedges(nhe)%verts                                                 !
       end do ! <------------------------------------------------------------------------------------+
    end do

    ! then, build closed hyperedges (cycles)
    do iedge = 1,brep%ne ! <------------------------------------+
       if ( visited(iedge) .or. .not.feat_edge(iedge) ) cycle   !
       !                                                        !
       if ( .not.allocated(hyperedges) .or. &                   !
            nhe + 1 > size(hyperedges) ) then ! <-------+       !
          call reallocate_hyperedges( &                 !       !
               hyperedges, &                            !       !
               nhe + PARAM_xtra_nhe )                   !       !
       end if ! <---------------------------------------+       !
       !                                                        !
       ivert = brep%edges(iedge)%halfedges(1)%orig              !
       nhe = nhe + 1                                            !
       PRINT *,'+1 CYCLE'                                       !
       hyperedges(nhe)%ne = 0                                   !
       call build_hyperedge( &                                  !
            brep, &                                             !
            ivert, &                                            !
            feat_edge, &                                        !
            visited, &                                          !
            feat_vert, &                                        !
            valence, &                                          !
            hyperedges(nhe) )                                   !
       !                                                        !
    end do ! <--------------------------------------------------+

    ! set incident hyperfaces
    do ihe = 1,nhe
       ihedg = hyperedges(ihe)%halfedges(:,1)
       hyperedges(ihe)%hyperfaces(1) = brep%faces(get_face(brep, ihedg))%hyperface
       ihedg = get_twin(ihedg)
       hyperedges(ihe)%hyperfaces(2) = brep%faces(get_face(brep, ihedg))%hyperface
    end do

    ! set edge -> hyperedge
    do ihe = 1,nhe
       do ied = 1,hyperedges(ihe)%ne
          brep%edges(hyperedges(ihe)%halfedges(1,ied))%hyperedge = ihe
       end do
    end do

  end subroutine get_hyperedges




  subroutine build_hyperedge( &
       brep, &
       istart, &
       feat_edge, &
       visited, &
       feat_vert, &
       valence, &
       he )
    use mod_util
    use mod_types_brep
    use mod_brep
    implicit none
    LOGICAL, PARAMETER :: DEBUG = .FALSE.
    type(type_BREP),      intent(in)    :: brep
    integer,              intent(in)    :: istart
    logical,              intent(in)    :: feat_edge(brep%ne)
    logical,              intent(inout) :: visited(brep%ne)
    logical,              intent(in)    :: feat_vert(brep%nv)
    integer,              intent(inout) :: valence(brep%nv)
    type(type_hyperedge), intent(inout) :: he
    integer                             :: ivert, ihedg(2)

    ivert = istart
    IF ( DEBUG ) PRINT *,'ISTART =',ISTART
    do
       IF ( DEBUG ) PRINT *,'IVERT =',IVERT
       ihedg = brep%verts(ivert)%halfedge
       do ! <-------------------------------------------------------------------+
          IF ( DEBUG ) PRINT *,'IHEDG =',IHEDG
          ! traverse halfedges around vertex to find an unvisited feature edge  !
          if ( feat_edge(ihedg(1)) .and. .not.visited(ihedg(1)) ) exit          !
          if ( he%ne > 0 ) then ! <-------------------------+                   !
             ! traverse halfedges CCW                       !                   !
             ihedg = get_twin(get_prev(brep, ihedg))        !                   !
          else ! -------------------------------------------+                   !
             ! traverse halfedges CW                        !                   !
             ihedg = get_next(brep, get_twin(ihedg))        !                   !
          end if ! <----------------------------------------+                   !
       end do ! <---------------------------------------------------------------+

       visited(ihedg(1)) = .true.
       call insert_column_after( &
            he%halfedges, &
            2, &
            he%ne, &
            ihedg, &
            he%ne )
       IF ( DEBUG ) PRINT *,'   + EDGE',IHEDG
       valence(ivert) = valence(ivert) - 1

       ivert = get_dest(brep, ihedg)
       valence(ivert) = valence(ivert) - 1

       if ( feat_vert(ivert) .or. ivert == istart ) exit
    end do

  end subroutine build_hyperedge



  subroutine reverse_hyperedge(self)
    use mod_types_brep
    use mod_brep
    implicit none
    type(type_hyperedge), intent(inout) :: self
    integer                             :: iedge

    self%verts(1:2) = self%verts([2,1])
    self%halfedges(1:2,1:self%ne) = self%halfedges(1:2,self%ne:1:-1)
    do iedge = 1,self%ne
       self%halfedges(1:2,iedge) = get_twin(self%halfedges(1:2,iedge))
    end do

  end subroutine reverse_hyperedge


  subroutine get_hyperedge_length( &
       self, &
       brep, &
       len_edg, &
       len_tot )
    use mod_types_brep
    use mod_brep
    implicit none
    type(type_hyperedge), intent(in)  :: self
    type(type_brep),      intent(in)  :: brep
    real(kind=fp),        intent(out) :: len_edg(self%ne)
    real(kind=fp),        intent(out) :: len_tot
    integer                           :: iedge

    len_tot = 0._fp
    do iedge = 1,self%ne
       len_edg(iedge) = get_polyline_length(brep, self%halfedges(:,iedge))
       len_tot = len_tot + len_edg(iedge)
    end do

  end subroutine get_hyperedge_length










  subroutine transfer_hyperfaces( &
       from, &
       to, &
       n )
    implicit none
    integer, intent(in) :: n
    type(type_hyperface), allocatable, intent(inout) :: from(:)
    type(type_hyperface), allocatable, intent(inout) :: to(:)
    integer                                          :: i

    do i = 1,n
       to(i)%nf = from(i)%nf
       call move_alloc(from(i)%faces, to(i)%faces)
    end do

  end subroutine transfer_hyperfaces


  subroutine transfer_hyperedges( &
       from, &
       to, &
       n )
    implicit none
    integer, intent(in) :: n
    type(type_hyperedge), allocatable, intent(inout) :: from(:)
    type(type_hyperedge), allocatable, intent(inout) :: to(:)
    integer                                          :: i

    do i = 1,n
       to(i)%ne = from(i)%ne
       call move_alloc(from(i)%halfedges, to(i)%halfedges)
       to(i)%verts = from(i)%verts
       to(i)%hyperfaces = from(i)%hyperfaces
    end do

  end subroutine transfer_hyperedges



  subroutine reallocate_hyperedges( &
       hyperedges, &
       n )
    implicit none
    type(type_hyperedge), allocatable, intent(inout) :: hyperedges(:)
    integer,                           intent(in)    :: n
    type(type_hyperedge), allocatable                :: tmp(:)
    integer                                          :: nhe

    if ( .not.allocated(hyperedges) ) then
       allocate(hyperedges(n))
       return
    end if

    nhe = size(hyperedges)

    ! reallocate hyperedges with bigger size
    allocate(tmp(nhe))
    call transfer_hyperedges( &
         hyperedges, &
         tmp, &
         nhe )
    deallocate(hyperedges)
    allocate(hyperedges(n))
    call transfer_hyperedges( &
         tmp, &
         hyperedges, &
         nhe )
    deallocate(tmp)

  end subroutine reallocate_hyperedges







  subroutine free_hypergraph( &
       self )
    implicit none
    type(type_hypergraph), intent(inout) :: self
    integer                              :: ihe, ihf

    do ihe = 1,self%nhe
       if ( allocated(self%hyperedges(ihe)%halfedges) ) &
            deallocate(self%hyperedges(ihe)%halfedges)
    end do

     do ihf = 1,self%nhf
       if ( allocated(self%hyperfaces(ihf)%faces) ) &
            deallocate(self%hyperfaces(ihf)%faces)
    end do

    if ( allocated(self%hyperedges) ) deallocate(self%hyperedges)
    if ( allocated(self%hyperfaces) ) deallocate(self%hyperfaces)

    self%nhe = 0
    self%nhf = 0
    
  end subroutine free_hypergraph






 
















  subroutine generate_brep_mesh( &
       brep, &
       hmin, &
       hmax, &
       tol, &
       feat_edge, &
       feat_vert, &
       hyperedges, &
       nhe, &
       mesh )
    use mod_brep!hypergraph
    use mod_mesh
    use mod_util
    use mod_polynomial
    use mod_halfedge
    implicit none
    type(type_brep),         intent(in)    :: brep
    real(kind=fp),           intent(in)    :: hmin
    real(kind=fp),           intent(in)    :: hmax
    real(kind=fp),           intent(in)    :: tol
    logical,                 intent(in)    :: feat_edge(brep%ne)
    logical,                 intent(in)    :: feat_vert(brep%nv)
    integer,                 intent(in)    :: nhe
    type(type_hyperedge),    intent(inout) :: hyperedges(nhe)
    type(type_surface_mesh), intent(inout) :: mesh
    integer                                :: finf, fpts, fedg
    integer                                :: bv2mv(brep%nv), be2mv(2,brep%ne)
    integer, allocatable                   :: loc2glob(:), idsf(:), typf(:)
    real(kind=fp), allocatable             :: uvf(:,:,:), xyzf(:,:), uvi(:,:)
    integer, allocatable                   :: trif(:,:)
    integer                                :: npf, np, ntrif
    integer                                :: ifirstedge, ifirstpoint, sens, ifirst, ilast
    type(type_path)                        :: pathtmp
    integer                                :: iface, iwire, ihedg(2), ivert, ihype, iedge, i
    character(3)                           :: strnum
    
    ! write 'info' file (common to all faces)
    call get_free_unit(finf)
    open(unit=finf, file='../tmp/info.dat', action='write')
    write (finf,*) hmin
    write (finf,*) hmax
    write (finf,*) tol
    close(finf)

    allocate(loc2glob(100), idsf(100), typf(100))
    mesh%nv = 0
    mesh%nt = 0
    bv2mv(:) = 0
    be2mv(:,:) = 0

    faces : do iface = 1,brep%nf
       ! write polynomial parametric surface 
       call write_polynomial(brep%faces(iface)%surface%x, '../tmp/c.cheb')
       !
       ! open boundary points & edges files
       call get_free_unit(fpts)
       open(unit=fpts, file='../tmp/bpts.dat', action='write')
       call get_free_unit(fedg)
       open(unit=fedg, file='../tmp/bedg.dat', action='write')
       !
       npf = 0
       wires : do iwire = 0,brep%faces(iface)%ninner ! <-----------------------------------+
          ! first halfedge of the wire                                                     !
          if ( iwire == 0 ) then ! <-------------------+                                   !
             ihedg = brep%faces(iface)%outer           !                                   !
          else ! --------------------------------------+                                   !
             ihedg = brep%faces(iface)%inner(:,iwire)  !                                   !
          end if ! <-----------------------------------+                                   !
          ifirstedge = ihedg(1)                                                            !
          !                                                                                !
          ifirstpoint = npf + 1                                                            !
          ! traverse the wire                                                              !
          halfedges : do ! <------------------------------------------------------------+  !
             ! get polyline points                                                      !  !
             call get_polyline_endpoints( &
                  brep, &
                  ihedg, &
                  ifirst, &
                  ilast, &
                  sens, &
                  np )
             !
             if ( .not.allocated(loc2glob) .or. &
                  size(loc2glob) < npf + np - 1 ) then ! <------+
                call reallocate_list(loc2glob, npf + np + 100)  !
             end if ! <-----------------------------------------+
             !
             if ( .not.allocated(idsf) .or. &
                  size(idsf) < npf + np - 1 ) then ! <----------+
                call reallocate_list(idsf, npf + np + 100)      !
             end if ! <-----------------------------------------+
             !
             if ( .not.allocated(typf) .or. &
                  size(typf) < npf + np - 1 ) then ! <----------+
                call reallocate_list(typf, npf + np + 100)      !
             end if ! <-----------------------------------------+
             !
             if ( allocated(uvf) .and. size(uvf,3) < np ) deallocate(uvf)
             if ( .not.allocated(uvf) ) allocate(uvf(2,2,np))
             ! polyligne uv dans le sens du contour de la face
             if ( sens == 2 ) then
                uvf(1:2,1:2,1:np) = brep%edges(ihedg(1))%curve%polyline%uv&
                     (1:2,[2,1],ifirst:ilast)
             else
                uvf(1:2,1:2,1:np) = brep%edges(ihedg(1))%curve%polyline%uv&
                     (1:2,1:2,ifirst:ilast:-1)
             end if
             !
             if ( allocated(xyzf) .and. size(xyzf,2) < np ) deallocate(xyzf)
             if ( .not.allocated(xyzf) ) allocate(xyzf(3,np))
             xyzf(1:3,1:np) = brep%edges(ihedg(1))%curve%polyline%xyz&
                  (1:3,ifirst:ilast:(-1)**sens)
             !                                                                          !  !
             ! get brep vertex index of first point
             ivert = get_orig(brep, ihedg)
             if ( bv2mv(ivert) == 0 ) then ! <-----------------------------------+
                ! the brep vertex is not yet present in the mesh                 !
                bv2mv(ivert) = mesh%nv + 1                                       !
                if ( feat_vert(ivert) ) then ! <------------------------------+  !
                   ! the mesh vertex is locked on a brep vertex               !  !
                   idsf(npf+1) = ivert                                        !  !
                   typf(npf+1) = 0                                            !  !
                else ! -------------------------------------------------------+  !
                   if ( feat_edge(ihedg(1)) ) then ! <---------------------+  !  !
                      ! the mesh vertex is locked on a hyperedge           !  !  !
                      idsf(npf+1) = ihedg(1)                               !  !  !
                      typf(npf+1) = 1                                      !  !  !
                   else ! -------------------------------------------------+  !  !
                      ! the mesh vertex is free to move along a hyperface  !  !  !
                      idsf(npf+1) = iface                                  !  !  !
                      typf(npf+1) = 2                                      !  !  !
                   end if ! <----------------------------------------------+  !  !
                end if ! <----------------------------------------------------+  !
                !                                                                !
                ! append new vertex to the mesh                                  !
                call append_vertices( &                                          !
                     mesh, &                                                     !
                     xyzf(1:3,1), &                                              !
                     uvf(1:2,1:2,1), &                                           !
                     idsf(npf+1), &                                              !
                     typf(npf+1), &                                              !
                     1 )                                                         !
             else ! -------------------------------------------------------------+
                if ( feat_edge(ihedg(1)) ) then ! <------------------------+     !
                   if ( mesh%typ(bv2mv(ivert)) == 2 ) then ! <----------+  !     !
                      mesh%ids(bv2mv(ivert)) = ihedg(1)                 !  !     !
                      mesh%typ(bv2mv(ivert)) = 1                        !  !     !
                      mesh%uv(1:2,1:2,bv2mv(ivert)) = uvf(1:2,[2,1],1)  !  !     !
                   end if ! <-------------------------------------------+  !     !
                end if ! <-------------------------------------------------+     !
             end if ! <----------------------------------------------------------+
             loc2glob(npf+1) = bv2mv(ivert)
             !
             if ( be2mv(1,ihedg(1)) == 0 ) then ! <-----------------------------------------+
                ! the brep edge is not yet present in the mesh                              !
                be2mv(1,ihedg(1)) = mesh%nv + 1                                             !
                be2mv(2,ihedg(1)) = mesh%nv + np - 2                                        !
                !                                                                           !
                if ( feat_edge(ihedg(1)) ) then ! <-----------------------+                 !
                   ! the mesh vertices are locked on a hyperedge          !                 !
                   idsf(npf+2:npf+np-1) = ihedg(1)                        !                 !
                   typf(npf+2:npf+np-1) = 1                               !                 !
                else ! ---------------------------------------------------+                 !
                   ! the mesh vertices are free to move along a hyperface !                 !
                   idsf(npf+2:npf+np-1) = iface                           !                 !
                   typf(npf+2:npf+np-1) = 2                               !                 !
                end if ! <------------------------------------------------+                 !
                !                                                                           !
                if ( sens == 1 ) then ! <-----------------------------------+               !
                   loc2glob(npf+2:npf+np-1) = mesh%nv + [(i, i=np-2,1,-1)]  !               !
                   call append_vertices( &                                  !               !
                        mesh, &                                             !               !
                        xyzf(1:3,np-1:2:-1), &                              !               !
                        uvf(1:2,1:2,np-1:2:-1), &                           !               !
                        idsf(npf+2:npf+np-1), &                             !               !
                        typf(npf+2:npf+np-1), &                             !               !
                        np-2 )                                              !               !
                else ! -----------------------------------------------------+               !
                   loc2glob(npf+2:npf+np-1) = mesh%nv + [(i, i=1,np-2)]     !               !
                   call append_vertices( &                                  !               !
                        mesh, &                                             !               !
                        xyzf(1:3,2:np-1), &                                 !               !
                        uvf(1:2,1:2,2:np-1), &                              !               !
                        idsf(npf+2:npf+np-1), &                             !               !
                        typf(npf+2:npf+np-1), &                             !               !
                        np-2 )                                              !               !
                end if ! <--------------------------------------------------+               !
             else ! ------------------------------------------------------------------------+
                if ( sens == 1 ) then ! <-----------------------------------+               !
                   loc2glob(npf+2:npf+np-1) = &                             !               !
                        [(i, i=be2mv(2,ihedg(1)),be2mv(1,ihedg(1)),-1)]     !               !
                else ! -----------------------------------------------------+               !
                   loc2glob(npf+2:npf+np-1) = &                             !               !
                        [(i, i=be2mv(1,ihedg(1)),be2mv(2,ihedg(1)))]        !               !
                end if ! <--------------------------------------------------+               !
             end if ! <---------------------------------------------------------------------+
             !
             ! write points...                            
             do i = 1,np-1 ! <---------------------------+
                write (fpts,*) uvf(:,1,i)                !
             end do ! <----------------------------------+
             ! ...and edges
             do i = 1,np-2 ! <---------------------------+
                write (fedg,*) npf + i, npf + i + 1      !
             end do ! <----------------------------------+
             !
             npf = npf + np - 1
             !
             ! move on to the next halfedge on the wire
             ihedg = get_next(brep, ihedg)
             !
             if ( ihedg(1) == ifirstedge ) then ! <------+
                ! the wire is complete                   !
                write (fedg,*) npf, ifirstpoint          !
                exit                                     !
             else ! -------------------------------------+
                write (fedg,*) npf, npf+1                !
             end if ! <----------------------------------+
             !           
          end do halfedges
       end do wires
       !
       ! close files
       close(fpts)
       close(fedg)
       !
       ! run mesher
       PRINT *,'MESH FACE #',IFACE
       call system( &!'/stck/bandrieu/Bureau/MeshGen/./meshgen.out &
            '/home/bastien/MeshGen/./meshgen.out &
            & ../tmp/c.cheb &
            & ../tmp/bpts.dat &
            & ../tmp/bedg.dat &
            & ../tmp/info.dat &
            & ../tmp/tri.dat &
            & ../tmp/uv.dat &
            & ../tmp/xyz.dat' )
       write (strnum,'(i3.3)') iface
       ! read face submesh                                                                     !
       call read_triangles( &                                                                  !
            '../tmp/tri.dat', &                                                                !
            trif, &                                                                            !
            ntrif )                                                                            !
       !                                                                                       !
       call read_points( &                                                                     !
            '../tmp/uv.dat', &                                                                 !
            uvi, &                                                                             !
            2, &                                                                               !
            np )                                                                               !
       if ( allocated(uvf) ) deallocate(uvf)                                                   !
       allocate(uvf(2,2,np-npf))                                                               !
       uvf(1:2,1,1:np-npf) = uvi(1:2,npf+1:np)                                                 !
       !                                                                                       !
       call read_points( &                                                                     !
            '../tmp/xyz.dat', &                                                                !
            xyzf, &                                                                            !
            3, &                                                                               !
            np )                                                                               !
       !                                                                                       !
       if ( size(loc2glob) < np ) call reallocate_list(loc2glob, np)                           !
       loc2glob(npf+1:np) = mesh%nv + [(i, i=1,np-npf)]                                        !
       !                                                                                       !
       if ( size(idsf) < np ) call reallocate_list(idsf, np)                                   !
       idsf(npf+1:np) = iface                                                                  !
       !                                                                                       !
       if ( size(typf) < np ) call reallocate_list(typf, np)                                   !
       typf(npf+1:np) = 2                                                                      !
       !                                                                                       !
       do i = 1,3 ! <----------------------------------+                                       !
          trif(i,1:ntrif) = loc2glob(trif(i,1:ntrif))  !                                       !
       end do ! <--------------------------------------+                                       !
       ! add new triangles                                                                     !
       call append_triangles( &                                                                !
            mesh, &                                                                            !
            trif(1:3,1:ntrif), &                                                               !
            [(brep%faces(iface)%hyperface, i=1,ntrif)], &                                      !
            ntrif )                                                                            !
       !                                                                                       !
       ! add new mesh vertices                                                                 !
       call append_vertices( &                                                                 !
            mesh, &                                                                            !
            xyzf(1:3,npf+1:np), &                                                              !
            uvf(1:2,1:2,1:np-npf), &                                                           !
            idsf(npf+1:np), &                                                                  !
            typf(npf+1:np), &                                                                  !
            np-npf )                                                                           !
       !                                                                                       !
    end do faces ! <---------------------------------------------------------------------------+

    ! create a path for each hyperedge
    do ihype = 1,nhe
       pathtmp%nv = 0
       pathtmp%hyperedge = ihype
       !
       if ( hyperedges(ihype)%verts(1) < 1 ) then ! <---+
          ihedg = hyperedges(ihype)%halfedges(:,1)      !
          ivert = get_orig(brep, ihedg)                 !
          hyperedges(ihype)%verts(1) = ivert            !
       end if ! <---------------------------------------+
       !
       do iedge = 1,hyperedges(ihype)%ne ! <---------------+
          ihedg = hyperedges(ihype)%halfedges(:,iedge)     !
          sens = 1 + mod(ihedg(2),2)                       !
          !                                                !
          call insert_after( &                             !
               pathtmp%verts, &                            !
               pathtmp%nv, &                               !
               bv2mv(get_orig(brep, ihedg)), &             !
               pathtmp%nv )                                !
          !                                                !
          np = be2mv(2,ihedg(1)) - be2mv(1,ihedg(1)) + 1   !
          if ( sens == 1 ) then ! <------+                 !
             ifirst = be2mv(2,ihedg(1))  !                 !
             ilast = be2mv(1,ihedg(1))   !                 !
          else ! ------------------------+                 !
             ifirst = be2mv(1,ihedg(1))  !                 !
             ilast = be2mv(2,ihedg(1))   !                 !
          end if ! <---------------------+                 !
          !                                                !
          call insert_n_after( &                           !
               pathtmp%verts, &                            !
               pathtmp%nv, &                               !
               np, &                                       !
               [(i, i=ifirst,ilast,(-1)**sens)], &         !
               pathtmp%nv )                                !
       end do ! <------------------------------------------+
       !
       if ( hyperedges(ihype)%verts(2) > 0 ) then ! <---+
          call insert_after( &                          !
               pathtmp%verts, &                         !
               pathtmp%nv, &                            !
               bv2mv(hyperedges(ihype)%verts(2)), &     !
               pathtmp%nv )                             !
       end if ! <---------------------------------------+
       !
       call append_path( &
            mesh, &
            pathtmp )
       !
       ! compute discrete curvilinear abscissa along the path
       allocate(mesh%paths(mesh%npaths)%s(mesh%paths(mesh%npaths)%nv))
       mesh%paths(mesh%npaths)%s(1) = 0._fp
       do i = 2,mesh%paths(mesh%npaths)%nv ! <-------------------------------+
          mesh%paths(mesh%npaths)%s(i) = mesh%paths(mesh%npaths)%s(i-1) + &  !
               norm2( &                                                      !
               mesh%xyz(:,mesh%paths(mesh%npaths)%verts(i)) - &              !
               mesh%xyz(:,mesh%paths(mesh%npaths)%verts(i-1)) )              !
       end do ! <------------------------------------------------------------+
       !
    end do
    
  end subroutine generate_brep_mesh








  
end module mod_hypergraph
