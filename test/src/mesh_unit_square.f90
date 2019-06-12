program mesh_unit_square

  use mod_util
  use mod_mesh
  use mod_mesh_optimization

  implicit none

  type(type_surface_mesh) :: mesh
  character(99)           :: argstr
  logical, allocatable    :: visited(:,:)
  real(kind=fp)           :: scale, edgelen, min_edgelen, max_edgelen
  integer                 :: twin(2)
  integer                 :: iv, it
  integer                 :: stat
  real(kind=fp)           :: longest_edgelen, shortest_edgelen
  integer, dimension(2)   :: longest_edge, shortest_edge
  integer                 :: iteration, i, j

  if ( command_argument_count() < 1 ) then
     scale = 1._fp
  else
     call get_command_argument(1, argstr)
     write (argstr, *) scale
  end if

  

  ! intial triangulation
  print *,'intial triangulation'
  allocate(mesh%tri(3,200))
  mesh%tri(1:3,1) = [1,2,3]
  mesh%tri(1:3,2) = [2,4,3]
  mesh%nt = 2
  allocate(mesh%ihf(mesh%nt))
  mesh%ihf(1:mesh%nt) = 0

  allocate(mesh%xyz(3,100))
  do j = 1,2
     do i = 1,2
        mesh%xyz(1:3,i + 2*(j-1)) = real([(-1)**i, (-1)**j, 0], kind=fp)
     end do
  end do
  mesh%nv = 4
  allocate(mesh%uv(2,2,mesh%nv), mesh%typ(mesh%nv), mesh%ids(mesh%nv))

  ! make halfedge data structure
  print *,'make halfedge data structure...'
  call make_halfedges(mesh)

  ! optimize mesh
  print *,'optimize mesh...'
  allocate(visited(3,mesh%nt+200))
  min_edgelen = 0.5_fp*scale
  max_edgelen = 2._fp*scale
  outer : do iteration = 1,10000
     PRINT *,''
     PRINT *,'ITERATION #',ITERATION
     if ( size(visited,2) < mesh%nt ) then
        deallocate(visited)
        allocate(visited(3,mesh%nt + MESH_xtra_nt))
     end if
     visited(:,:) = .false.

     longest_edge(1:2) = 0
     shortest_edge(1:2) = 0
     longest_edgelen = max_edgelen
     shortest_edgelen = min_edgelen
     do it = 1,mesh%nt
        do iv = 1,3
           if ( visited(iv,it) ) cycle
           visited(iv,it) = .true.
           twin = get_twin(mesh, [iv,it])
           if ( all(twin > 0) ) visited(twin(1), twin(2)) = .true.

           ! compute edge length
           edgelen = sum((&
                mesh%xyz(1:3,mesh%tri(iv,it)) - &
                mesh%xyz(1:3,mesh%tri(1+mod(iv,3),it))&
                )**2)

           if ( edgelen > longest_edgelen ) then
              longest_edge(1:2) = [iv,it]
              longest_edgelen = edgelen
           elseif ( edgelen < shortest_edgelen ) then
              shortest_edge(1:2) = [iv,it]
              shortest_edgelen = edgelen
           end if
        end do
     end do

     if ( all(shortest_edge > 0) ) then
        iv = shortest_edge(1)
        it = shortest_edge(2)
        print *,'collapse edge ',mesh%tri([iv, 1+mod(iv,3)],it)
        call collapse_edge( &
             mesh, &
             [iv,it], &
             stat )
        if ( stat > 0 ) shortest_edge(1:2) = 0
     end if

     if ( all(longest_edge > 0) ) then
        iv = longest_edge(1)
        it = longest_edge(2)
        print *,'split edge ',mesh%tri([iv, 1+mod(iv,3)],it)
        call split_edge( &
             mesh, &
             [iv,it] )
     end if

     if ( any(shortest_edge < 1) .and. any(longest_edge < 1) ) exit

  end do outer


   ! export optimized mesh
  print *,'export optimized mesh...'
  call write_obj_mesh( &
       mesh, &
       'mesh_unit_square/mesh.obj' )
  call write_mesh_files( &
       mesh, &
       filetri='mesh_unit_square/tri.dat', &
       filexyz='mesh_unit_square/xyz.dat' )
  call write_halfedges( &
       mesh, &
       'mesh_unit_square/v2h.dat', &
       'mesh_unit_square/twin.dat' )


  print *,'done.'

end program mesh_unit_square
