module mod_types_brep

  use mod_diffgeom
  use mod_types_intersection
  
  implicit none

  type type_BREPface
     type(type_surface), pointer            :: surface => null()
     integer                                :: outer(2) = 0 ! [#edge, #halfedge]
     integer, allocatable                   :: inner(:,:)
     integer                                :: ninner = 0

     integer                                :: hyperface = 0
  end type type_BREPface

  type type_BREPhalfedge
     integer                                :: face = 0
     integer                                :: orig = 0
     integer                                :: prev(2) = 0 ! [#edge, #halfedge]
     integer                                :: next(2) = 0 ! [#edge, #halfedge]
  end type type_BREPhalfedge
  
  type type_BREPedge
     type(type_intersection_curve), pointer :: curve => null()
     integer                                :: isplit = 0
     type(type_BREPhalfedge)                :: halfedges(2)
     
     integer                                :: hyperedge = 0
  end type type_BREPedge

  
  type type_BREPvertex
     type(type_intersection_point), pointer :: point => null()
     integer                                :: halfedge(2) = 0 ! [#edge, #halfedge] (outgoing)
  end type type_BREPvertex


  type type_BREP
     ! Boundary Representation
     integer                                :: nv = 0, ne = 0, nf = 0
     type(type_BREPvertex), allocatable     :: verts(:)
     type(type_BREPedge), allocatable       :: edges(:)
     type(type_BREPface), allocatable       :: faces(:)
  end type type_BREP
  
end module mod_types_brep
