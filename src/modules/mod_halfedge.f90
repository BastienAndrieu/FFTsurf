module mod_halfedge

  use mod_mesh
  use mod_types_brep

  implicit none

  interface get_prev
     module procedure get_prev_mesh, get_prev_brep
  end interface get_prev

  interface get_next
     module procedure get_next_mesh, get_next_brep
  end interface get_next
  
  interface get_face
     module procedure get_face_mesh, get_face_brep
  end interface get_face

  interface get_twin
     module procedure get_twin_mesh, get_twin_brep
  end interface get_twin

  interface get_orig
     module procedure get_orig_mesh, get_orig_brep
  end interface get_orig

  interface get_dest
     module procedure get_dest_mesh, get_dest_brep
  end interface get_dest

contains
  
  function get_prev_mesh( &
       ihedg )
    ! returns the index of an halfedge's previous halfedge
    implicit none
    integer,                 intent(in) :: ihedg(2)
    integer                             :: get_prev_mesh(2)

    get_prev_mesh = [1 + mod(ihedg(1)+3-2,3), ihedg(2)] ! only for triangles
  end function get_prev_mesh


  function get_next_mesh( &
       ihedg )
    ! returns the index of an halfedge's next halfedge
    implicit none
    integer,                 intent(in) :: ihedg(2)
    integer                             :: get_next_mesh(2)

    get_next_mesh = [1 + mod(ihedg(1),3), ihedg(2)] ! only for triangles
  end function get_next_mesh
  

  function get_face_mesh( &
       ihedg )
    ! returns the index of an halfedge's incident face
    implicit none
    integer,         intent(in) :: ihedg(2)
    integer                     :: get_face_mesh

    get_face_mesh = ihedg(2)
  end function get_face_mesh

  
  function get_twin_mesh( &
       mesh, &
       ihedg )
    ! returns the index of an halfedge's twin halfedge
    implicit none
    type(type_surface_mesh), intent(in) :: mesh
    integer,                 intent(in) :: ihedg(2)
    integer                             :: get_twin_mesh(2)

    get_twin_mesh = mesh%twin(:,ihedg(1),ihedg(2))
  end function get_twin_mesh

  
  function get_orig_mesh( &
       mesh, &
       ihedg )
    ! returns the index of an halfedge's origin vertex
    implicit none
    type(type_surface_mesh), intent(in) :: mesh
    integer,                 intent(in) :: ihedg(2)
    integer                             :: get_orig_mesh

    get_orig_mesh = mesh%tri(ihedg(1),ihedg(2))
  end function get_orig_mesh


  function get_dest_mesh( &
       mesh, &
       ihedg )
    ! returns the index of an halfedge's destination vertex
    implicit none
    type(type_surface_mesh), intent(in) :: mesh
    integer,                 intent(in) :: ihedg(2)
    integer                             :: get_dest_mesh

    get_dest_mesh = get_orig_mesh(mesh, get_next_mesh(ihedg))
  end function get_dest_mesh




  

  function get_prev_brep( &
       brep, &
       ihedg )
    ! returns the index of an halfedge's previous halfedge
    implicit none
    type(type_BREP), intent(in) :: brep
    integer,         intent(in) :: ihedg(2)
    integer                     :: get_prev_brep(2)

    get_prev_brep = brep%edges(ihedg(1))%halfedges(ihedg(2))%prev
  end function get_prev_brep


  function get_next_brep( &
       brep, &
       ihedg )
    ! returns the index of an halfedge's next halfedge
    implicit none
    type(type_BREP), intent(in) :: brep
    integer,         intent(in) :: ihedg(2)
    integer                     :: get_next_brep(2)

    get_next_brep = brep%edges(ihedg(1))%halfedges(ihedg(2))%next
  end function get_next_brep

  
  function get_twin_brep( &
       ihedg )
    ! returns the index of an halfedge's twin halfedge
    implicit none
    integer,         intent(in) :: ihedg(2)
    integer                     :: get_twin_brep(2)

    get_twin_brep = [ihedg(1), 1 + mod(ihedg(2),2)]
  end function get_twin_brep
  

  function get_face_brep( &
       brep, &
       ihedg )
    ! returns the index of an halfedge's incident face
    implicit none
    type(type_BREP), intent(in) :: brep
    integer,         intent(in) :: ihedg(2)
    integer                     :: get_face_brep

    get_face_brep = brep%edges(ihedg(1))%halfedges(ihedg(2))%face
  end function get_face_brep

  
  function get_orig_brep( &
       brep, &
       ihedg )
    ! returns the index of an halfedge's origin vertex
    implicit none
    type(type_BREP), intent(in) :: brep
    integer,         intent(in) :: ihedg(2)
    integer                     :: get_orig_brep

    get_orig_brep = brep%edges(ihedg(1))%halfedges(ihedg(2))%orig
  end function get_orig_brep


  function get_dest_brep( &
       brep, &
       ihedg )
    ! returns the index of an halfedge's destination vertex
    implicit none
    type(type_BREP), intent(in) :: brep
    integer,         intent(in) :: ihedg(2)
    integer                     :: get_dest_brep

    get_dest_brep = get_orig_brep(brep, get_next_brep(brep, ihedg))
  end function get_dest_brep
  
end module mod_halfedge
