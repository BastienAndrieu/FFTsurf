program import_geometry

  use mod_util
  use mod_init
  use mod_options
  use mod_diffgeom
  use mod_import
  use mod_types_brep
  use mod_types_intersection
  use mod_hypergraph
  use mod_mesh
  use mod_brep

  implicit none

  type(type_options)                        :: options
  integer                                   :: narg, fid
  character(99)                             :: fileoptions
  type(type_surface), allocatable           :: surf(:)
  integer                                   :: nsurf
  type(type_surface_mesh)                   :: mesh
  type(type_brep)                           :: brep
  type(type_intersection_data), target      :: interdata
  type(type_hypergraph)                     :: hypergraph
  integer                                   :: i, j


  narg = command_argument_count()
  if ( narg < 1 ) then
     STOP 'you must provide an option file'
  else
     call get_command_argument(1, fileoptions)
  end if


  call init_from_mesh( &
       fileoptions, &
       options, &
       surf, &
       nsurf, &
       interdata, &
       brep, &
       hypergraph, &
       mesh )

  call write_mesh_files( &
       mesh, &
       'import_geometry/mesh/tri.dat', &
       'import_geometry/mesh/xyz.dat', &
       'import_geometry/mesh/uv.dat', &
       'import_geometry/mesh/idstyp.dat', &
       'import_geometry/mesh/paths.dat' )


  call write_brep_files( &
       brep, &
       'import_geometry/brep/verts.dat', &
       'import_geometry/brep/edges.dat', &
       'import_geometry/brep/faces.dat' )

  call get_free_unit( fid )
  open(unit=fid, file='import_geometry/brep/hyperfaces.dat', action='write')
  write (fid,*) hypergraph%nhf
  do i = 1,hypergraph%nhf
     write (fid,*) hypergraph%hyperfaces(i)%nf
     write (fid,*) hypergraph%hyperfaces(i)%faces(1:hypergraph%hyperfaces(i)%nf)
  end do
  close(fid)

  open(unit=fid, file='import_geometry/brep/hyperedges.dat', action='write')
  write (fid,*) hypergraph%nhe
  do i = 1,hypergraph%nhe
     write (fid,*) hypergraph%hyperedges(i)%ne
     write (fid,*) hypergraph%hyperedges(i)%verts
     write (fid,*) hypergraph%hyperedges(i)%hyperfaces
     do j = 1,hypergraph%hyperedges(i)%ne
        write (fid,*) hypergraph%hyperedges(i)%halfedges(1:2,j)
     end do
  end do
  close(fid)



  PRINT *,MESH%NV,' VERTICES,', MESH%NT, ' TRIANGLES'
  PRINT *,INTERDATA%NP,' POINTS,', INTERDATA%NC,' CURVES'
  PRINT *,BREP%NV,' BREP VERTS,',BREP%NE,' BREP EDGES,',BREP%NF,' BREP FACES'
  PRINT *,HYPERGRAPH%NHE,' HYPEREDGES, ',HYPERGRAPH%NHF,' HYPERFACES'

  
  
  ! *********************
  do i = 1,mesh%nv ! <-----------------------------------------------------+
     select case ( mesh%typ(i) ) ! <-----------------------------------+   !
     case (0) ! -------------------------------------------------------+   !
        mesh%xyz(:,i) = brep%verts(mesh%ids(i))%point%xyz              !   !
     case (1) ! -------------------------------------------------------+   !
        j = brep%edges(mesh%ids(i))%halfedges(2)%face                  !   !
        call eval( &                                                   !   !
             mesh%xyz(:,i), &                                          !   !
             brep%faces(j)%surface, &                                  !   !
             mesh%uv(:,1,i) )                                          !   !
     case (2) ! -------------------------------------------------------+   !
        j = mesh%ids(i)                                                !   !
        call eval( &                                                   !   !
             mesh%xyz(:,i), &                                          !   !
             brep%faces(j)%surface, &                                  !   !
             mesh%uv(:,1,i) )                                          !   !
     end select ! <----------------------------------------------------+   !
  end do ! <---------------------------------------------------------------+
  
  !call write_connectivity( &
  !     'import_geometry/mesh/', &
  !     mesh, &
  !     1 )
  !call write_face_ref( &
  !     'import_geometry/mesh/', &
  !     mesh, &
  !     1 )
  call write_xyz_positions( &
       'import_geometry/mesh/', &
       mesh, &
       1 )

end program import_geometry
