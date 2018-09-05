module mod_mesh

  use mod_math
  
  implicit none
  
  type type_path ! trouver nom + explicite...
     integer                    :: nv = 0
     integer, allocatable       :: verts(:)
     real(kind=fp), allocatable :: s(:)
     integer                    :: hyperedge = 0
  end type type_path
  
end module mod_mesh
