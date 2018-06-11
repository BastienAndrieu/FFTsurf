module mod_types_dcel
  ! Defines derived types for implementing a Doubly-Connected Edge List (DCEL) data structure
  ! ( see "Computational geometry", De Berg et al. (2000), pp.29-33. )
  implicit none
  
  ! DERIVED TYPES ===============================================================================================================
  ! --------------------------------------------------------------------------------------------------------------------------- !
  type type_dcel_vertex                                                                                                         !
     type(type_dcel_halfedge), pointer     :: edge => null() ! pointer to one outgoing incident dcel_halfedge                   !
  end type type_dcel_vertex                                                                                                     !
  ! --------------------------------------------------------------------------------------------------------------------------- !


  ! --------------------------------------------------------------------------------------------------------------------------- !
  type type_dcel_halfedge                                                                                                       !
     type(type_dcel_face),     pointer     :: face => null() ! pointer to incident dcel_face (to the left)                      !
     type(type_dcel_halfedge), pointer     :: twin => null() ! pointer to twin (opposite) dcel_halfedge                         !
     type(type_dcel_halfedge), pointer     :: prev => null() ! pointer to previous dcel_halfedge                                !
     type(type_dcel_halfedge), pointer     :: next => null() ! pointer to next dcel_halfedge                                    !
     type(type_dcel_vertex),   pointer     :: orig => null() ! pointer to origin dcel_vertex                                    !
  end type type_dcel_halfedge                                                                                                   !
  ! --------------------------------------------------------------------------------------------------------------------------- !

  
  ! --------------------------------------------------------------------------------------------------------------------------- !
  ! Workaround for making arrays of pointers to dcel_halfedges                                                                  !
  type ptr_dcel_halfedge                                                                                                        !
     type(type_dcel_halfedge), pointer     :: ptr => null()                                                                     !
  end type ptr_dcel_halfedge                                                                                                    !
  ! --------------------------------------------------------------------------------------------------------------------------- !


  ! --------------------------------------------------------------------------------------------------------------------------- !
  type type_dcel_face                                                                                                           !
     type(type_dcel_halfedge), pointer     :: outer => null() ! pointer to one dcel_halfedge on the outer boundary              !
     type(ptr_dcel_halfedge), allocatable  :: inner(:)        ! pointers to dcel_halfedges on each potential inner boundary     !
  end type type_dcel_face                                                                                                       !
  ! --------------------------------------------------------------------------------------------------------------------------- !


  ! --------------------------------------------------------------------------------------------------------------------------- !
  type type_DCEL                                                                                                                !
     integer                               :: nv = 0   ! actual nb. of vertices                                                 !
     integer                               :: ne = 0   ! actual nb. of edges                                                    !
     integer                               :: nf = 0   ! actual nb. of faces                                                    !
     type(type_dcel_vertex),   allocatable :: verts(:) ! list of vertices                                                       !
     type(type_dcel_halfedge), allocatable :: edges(:) ! list of edges                                                          !
     type(type_dcel_face),     allocatable :: faces(:) ! list of faces                                                          !
  end type type_DCEL                                                                                                            !
  ! --------------------------------------------------------------------------------------------------------------------------- !
  ! =============================================================================================================================
end module mod_types_dcel
