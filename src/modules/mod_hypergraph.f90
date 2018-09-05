module mod_hypergraph

  implicit none
  
  type type_hyperface
     integer              :: nf = 0
     integer, allocatable :: faces(:)
  end type type_hyperface

  type type_hyperedge
     integer              :: ne = 0
     integer, allocatable :: edges(:)
     integer              :: vertices(2)
  end type type_hyperedge
  
contains

  subroutine get_feature_edges(brep, fe)
    use mod_brep
    ! returns the list of feature edges in a BREP, i.e. non smooth or boundary edges
    implicit none
    type(type_BREP), intent(in)  :: brep
    logical,         intent(out) :: fe(brep%ne)
    integer                      :: ie

    do ie = 1,brep%ne
       fe(ie) = ( .not.is_smooth(brep%edges(ie)) .or. .not.associated(brep%edges(ie)%halfedge%twin) )
    end do

  end subroutine get_feature_edges




  subroutine get_feature_vertices(brep, fe, fv, valence)
    use mod_brep
    ! returns the list of feature edges in a BREP, i.e. incident to 1 or more than 2 feature edges,
    ! or located at a tangent discontinuity between two adjacent feature edges
    implicit none
    real(kind=fp), parameter                 :: EPScosangle = real(1.d-3, kind=fp)
    type(type_BREP),      intent(in), target :: brep
    logical,              intent(in)         :: fe(brep%ne)
    logical,              intent(out)        :: fv(brep%dcel%nv)
    integer,              intent(out)        :: valence(brep%dcel%nv)
    real(kind=fp)                            :: cosangle
    integer, allocatable                     :: v2e(:)
    integer                                  :: nv2e
    integer                                  :: iv, jv, ie, je, pair(2,2), verts(2)

    fv(:) = .false.       
    outer_loop : do iv = 1,brep%dcel%nv ! <--------------------------------------------+
       valence(iv) = 0                                                                 !
       call get_v2e(brep%dcel, iv, v2e, nv2e)                                          !
       do je = 1,nv2e ! <---------------------------------------------------+          !
          ie = v2e(je)                                                      !          !
          if ( fe(ie) ) then ! <----------------------------------------+   !          !
             valence(iv) = valence(iv) + 1                              !   !          !
             if ( valence(iv) <= 2 ) then ! <-----------------------+   !   !          !
                verts = get_BREPedge_vertices(brep%edges(ie))       !   !   !          !
                do jv = 1,2 ! <-----------------+                   !   !   !          !
                   if ( verts(jv) == iv ) exit  !                   !   !   !          !
                end do ! <----------------------+                   !   !   !          !
                pair(valence(iv),1:2) = [ie, jv]                    !   !   !          !
             end if ! <---------------------------------------------+   !   !          !
          end if ! <----------------------------------------------------+   !          !
       end do ! <-----------------------------------------------------------+          !
       !                                                                               !      
       if ( valence(iv) == 0 ) cycle outer_loop                                        !
       !                                                                               !
       if ( valence(iv) == 2 ) then ! <------------------------------------------+     !
          ! check angle between left and right tangent directions at the vertex  !     !
          cosangle = dot_product( &                                              !     !
               brep%edges(pair(1,1))%tangents(:,pair(1,2)), &                    !     !
               brep%edges(pair(2,1))%tangents(:,pair(2,2)) )                     !     !
          if ( abs(cosangle) > 1._fp - EPScosangle ) cycle outer_loop            !     !
       end if ! <----------------------------------------------------------------+     !
       !                                                                               !
       fv(iv) = .true.                                                                 !
       !                                                                               !
    end do outer_loop ! <--------------------------------------------------------------+
    if ( allocated(v2e) ) deallocate(v2e)
    
  end subroutine get_feature_vertices



  

  subroutine get_hyperfaces(brep, feat_edg, hyperfaces, nhf)
    use mod_brep
    implicit none
    integer, parameter                               :: PARAM_xtra_nhf = 10
    type(type_BREP),                   intent(in)    :: brep
    logical,                           intent(in)    :: feat_edg(brep%ne)
    type(type_hyperface), allocatable, intent(inout) :: hyperfaces(:)
    integer,                           intent(out)   :: nhf
    type(type_hyperface), allocatable                :: tmp(:)
    logical                                          :: unvisited(brep%dcel%nf)
    integer                                          :: ifa
    
    unvisited(:) = .true.
    
    do while ( any(unvisited) )
       do ifa = 1,brep%dcel%nf
          ! look for an unvisited face to start a new hyperface
          if ( unvisited(ifa) ) exit
       end do

       if ( .not.allocated(hyperfaces) ) allocate(hyperfaces(PARAM_xtra_nhf))
       
       if ( nhf + 1 > size(hyperfaces) ) then
          ! reallocate hyperfaces with bigger size
          call move_alloc(from=hyperfaces, to=tmp)
          allocate(hyperfaces(nhf + PARAM_xtra_nhf), source=tmp)
       end if

       nhf = nhf + 1
       hyperfaces(nhf)%nf = 0
       
       call build_hyperface(brep, ifa, feat_edg, hyperfaces(nhf), unvisited)
       !print *,''
    end do
    
  end subroutine get_hyperfaces


  

  recursive subroutine build_hyperface(brep, ifa, feat_edg, hf, unvisited)
    use mod_util
    use mod_brep
    implicit none
    type(type_BREP),      intent(in), target :: brep
    integer,              intent(in)         :: ifa
    logical,              intent(in)         :: feat_edg(brep%ne)
    type(type_hyperface), intent(inout)      :: hf
    logical,              intent(inout)      :: unvisited(brep%dcel%nf)
    type(type_DCELhalfedge)                  :: e
    integer                                  :: i, jfa, vefirst

    if ( .not.unvisited(ifa) ) return
    
    call insert_after( &
         hf%faces, &
         hf%nf, &
         ifa, &
         hf%nf )
    
    unvisited(ifa) = .false.

    do i = 1,1 + brep%dcel%faces(ifa)%ninner ! <--------------------------------------+
       if ( i == 1 ) then ! <---------------------------------------+                 !
          ! traverse outer loop                                     !                 !
          e = brep%dcel%halfedges(brep%dcel%faces(ifa)%outer)       !                 !
       else ! ------------------------------------------------------+                 !
          ! traverse potential inner loops                          !                 !
          e = brep%dcel%halfedges(brep%dcel%faces(ifa)%inner(i-1))  !                 !
       end if ! <---------------------------------------------------+                 !
       vefirst = e%orig                                                               !
       while_loop : do ! <---------------------------------------------------------+  !
          if ( associated(e%twin) .and. .not.feat_edg(e%edge) ) then ! <--------+  !  !
             jfa = e%twin%face                                                  !  !  !
             call build_hyperface(brep, jfa, feat_edg, hf, unvisited)           !  !  !
          end if ! <------------------------------------------------------------+  !  !
          !                                                                        !  !
          if ( e%next%orig == vefirst ) exit while_loop                            !  !
          e = e%next                                                               !  !
       end do while_loop ! <-------------------------------------------------------+  !
    end do ! <------------------------------------------------------------------------+

  end subroutine build_hyperface





  subroutine get_hyperedges(brep, feat_edg, feat_ver, valence, hyperedges, nhe)
    use mod_util
    use mod_brep
    implicit none
    integer, parameter                               :: PARAM_xtra_nhe = 10
    type(type_BREP), target,           intent(in)    :: brep
    logical,                           intent(in)    :: feat_edg(brep%ne)
    logical,                           intent(in)    :: feat_ver(brep%dcel%nv)
    integer,                           intent(inout) :: valence(brep%dcel%nv)
    type(type_hyperedge), allocatable, intent(inout) :: hyperedges(:)
    integer,                           intent(out)   :: nhe
    type(type_hyperedge), allocatable                :: tmp(:)
    logical                                          :: unvisited(brep%ne)
    type(type_DCELhalfedge), pointer                 :: he => null()
    integer                                          :: iv
    
    ! first, do paths
    unvisited(:) = .true.
    do while ( any(pack(valence, mask=feat_ver) > 0) )
       do iv = 1,brep%dcel%nv
          ! look for a feature vertex with valence > 0
          if ( feat_ver(iv) .and. valence(iv) > 0 ) exit
       end do

       if ( .not.allocated(hyperedges) ) allocate(hyperedges(PARAM_xtra_nhe))

       if ( nhe + 1 > size(hyperedges) ) then
          ! reallocate hyperfaces with bigger size
          call move_alloc(from=hyperedges, to=tmp)
          allocate(hyperedges(nhe + PARAM_xtra_nhe), source=tmp)
       end if

       nhe = nhe + 1
       hyperedges(nhe)%ne = 0
       hyperedges(nhe)%vertices(1) = iv

       he => brep%dcel%halfedges(brep%dcel%vertices(iv)%halfedge)
       do
          call insert_after( &
               hyperedges(nhe)%edges, &
               hyperedges(nhe)%ne, &
               he%edge, &
               hyperedges(nhe)%ne )
       end do
       
    end do
    
  end subroutine get_hyperedges
  
  
end module mod_hypergraph
