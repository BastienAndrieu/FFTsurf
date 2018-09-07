module mod_hypergraph

  implicit none
  
  type type_hyperface
     integer              :: nf = 0
     integer, allocatable :: faces(:)
  end type type_hyperface

  type type_hyperedge
     integer              :: ne = 0
     integer, allocatable :: edges(:,:)
     integer              :: verts(2) = 0
  end type type_hyperedge
  
contains

  subroutine get_feature_edges( &
       brep, &
       feat_edge )
    use mod_brep2
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
    use mod_brep2
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
    use mod_brep2
    implicit none
    integer, parameter                               :: PARAM_xtra_nhf = 10
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
          brep%faces(ifa)%hyperface = ihf
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
    use mod_brep2
    implicit none
    type(type_BREP),      intent(in)    :: brep
    integer,              intent(in)    :: ifa
    logical,              intent(in)    :: feat_edge(brep%ne)
    type(type_hyperface), intent(inout) :: hf
    logical,              intent(inout) :: visited(brep%nf)
    integer                             :: iedge(2)
    integer                             :: iloop, istart, jfa

    if ( visited(ifa) ) return
    
    call insert_after( &
         hf%faces, &
         hf%nf, &
         ifa, &
         hf%nf )
    
    visited(ifa) = .true.

    do iloop = 1,1 + brep%faces(ifa)%ninner ! <---------------------------------------+
       if ( iloop == 1 ) then ! <--------------------------------------+              !
          ! traverse outer loop                                        !              !
          iedge = brep%faces(ifa)%outer                                !              !
       else ! ---------------------------------------------------------+              !
          ! traverse potential inner loops                             !              !
          iedge = brep%faces(ifa)%inner(:,iloop-1)                     !              !
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
    use mod_brep2
    implicit none
    integer, parameter                               :: PARAM_xtra_nhe = 10
    type(type_BREP),                   intent(inout) :: brep
    logical,                           intent(in)    :: feat_edge(brep%ne)
    logical,                           intent(in)    :: feat_vert(brep%nv)
    integer,                           intent(inout) :: valence(brep%nv)
    type(type_hyperedge), allocatable, intent(inout) :: hyperedges(:)
    integer,                           intent(out)   :: nhe
    logical                                          :: visited(brep%ne)
    integer                                          :: ivert, iedge

    visited(:) = .false.
    nhe = 0
    ! first, build open hyperedges (paths)
    do ivert = 1,brep%nv ! <----------------------------------------------------------------------+
       if ( .not.feat_vert(ivert) ) cycle                                                         !
       do while ( valence(ivert) > 0 )                                                            !
          !                                                                                       !
          if ( .not.allocated(hyperedges) .or. &                                                  !
               nhe + 1 > size(hyperedges) ) then ! <-------+                                      !
             call reallocate_hyperedges( &                 !                                      !
                  hyperedges, &                            !                                      !
                  nhe + PARAM_xtra_nhe )                   !                                      !
          end if ! <---------------------------------------+                                      !
          !                                                                                       !
          nhe = nhe + 1                                                                           !
          PRINT *,'+1 PATH'                                                                       !
          hyperedges(nhe)%ne = 0                                                                  !
          call build_hyperedge( &                                                                 !
               brep, &                                                                            !
               ivert, &                                                                           !
               feat_edge, &                                                                       !
               visited, &                                                                         !
               feat_vert, &                                                                       !
               valence, &                                                                         !
               hyperedges(nhe) )                                                                  !
          !                                                                                       !
          ! set boundary vertices indices                                                         !
          hyperedges(nhe)%verts(1) = get_orig(brep, hyperedges(nhe)%edges(:,1))                   !
          hyperedges(nhe)%verts(2) = get_dest(brep, hyperedges(nhe)%edges(:,hyperedges(nhe)%ne))  !
          PRINT *,'   VERTS =',hyperedges(nhe)%verts                                              !
       end do ! <---------------------------------------------------------------------------------+
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
    use mod_brep2
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
            he%edges, &
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
       call move_alloc(from(i)%edges, to(i)%edges)
       to(i)%verts = from(i)%verts
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
       
  
end module mod_hypergraph
