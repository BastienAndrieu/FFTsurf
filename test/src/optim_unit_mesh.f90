program optim_unit_mesh

  use mod_util
  use mod_mesh
  use mod_mesh_optimization

  implicit none
    
  logical :: DEBUG
  type(type_surface_mesh)    :: mesh
  integer                    :: tri(3)
  integer, allocatable       :: tritmp(:,:)
  real(kind=fp)              :: xyz(3)
  real(kind=fp), allocatable :: xyztmp(:,:)
  integer                    :: fid, err

  logical, allocatable       :: visited(:,:)
  real(kind=fp)              :: scale, edgelen, min_edgelen, max_edgelen
  integer                    :: twin(2)
  integer                    :: iv, it
  integer                    :: stat
  real(kind=fp)              :: longest_edgelen, shortest_edgelen
  integer, dimension(2)      :: longest_edge, shortest_edge
  integer                    :: iteration, max_iteration

  if ( command_argument_count() < 1 ) then
     DEBUG = .false.
  else
     DEBUG = .true.     
  end if


  ! load initial mesh data
  print *,'load triangles...'
  call get_free_unit(fid)
  allocate(mesh%tri(3,200))
  open(unit=fid, file='optim_unit_mesh/initial/tri.dat', action='read')
  do
     read (fid,*,iostat=err) tri
     if ( err /=0 ) exit
     if ( mesh%nt + 1 > size(mesh%tri,2) ) then
        call move_alloc(from=mesh%tri, to=tritmp)
        allocate(mesh%tri(3,mesh%nt + 200))
        mesh%tri(1:3,1:mesh%nt) = tritmp(1:3,1:mesh%nt)
        deallocate(tritmp)
     end if
     mesh%nt = mesh%nt + 1
     mesh%tri(1:3,mesh%nt) = tri(1:3)
  end do
  close(fid)
  allocate(mesh%ihf(mesh%nt))
  mesh%ihf(1:mesh%nt) = 0

  print *,'load vertex coordinates...'
  allocate(mesh%xyz(3,100))
  open(unit=fid, file='optim_unit_mesh/initial/xyz.dat', action='read')
  do
     read (fid,*,iostat=err) xyz
     if ( err /=0 ) exit
     if ( mesh%nv + 1 > size(mesh%xyz,2) ) then
        call move_alloc(from=mesh%xyz, to=xyztmp)
        allocate(mesh%xyz(3,mesh%nv + 100))
        mesh%xyz(1:3,1:mesh%nv) = xyztmp(1:3,1:mesh%nv)
        deallocate(xyztmp)
     end if
     mesh%nv = mesh%nv + 1
     mesh%xyz(1:3,mesh%nv) = xyz(1:3)
  end do
  close(fid)
  allocate(mesh%uv(2,2,mesh%nv), mesh%typ(mesh%nv), mesh%ids(mesh%nv))


  ! flatten z-direction
  mesh%xyz(3,1:mesh%nv) = 0._fp

  ! make halfedge data structure
  print *,'make halfedge data structure...'
  call make_halfedges(mesh)

  ! export initial mesh
  print *,'export initial mesh...'
  call write_obj_mesh( &
       mesh, &
       'optim_unit_mesh/initial_mesh.obj' )

  call write_mesh_files( &
       mesh, &
       filetri='optim_unit_mesh/initial/tri.dat', &
       filexyz='optim_unit_mesh/initial/xyz.dat' )
  call write_halfedges( &
       mesh, &
       'optim_unit_mesh/initial/v2h.dat', &
       'optim_unit_mesh/initial/twin.dat' )




  !*****
  scale = 0.5_fp
  scale = scale**2
  min_edgelen = 0.5_fp*scale
  max_edgelen = 2._fp*scale
  !*****

  allocate(visited(3,mesh%nt+200))
  ! optimize mesh
  print *,'optimize mesh...'
  max_iteration = 3000*mesh%nt
  outer : do iteration = 1,max_iteration
     PRINT *,''
     PRINT *,'ITERATION #',ITERATION
     !PRINT *,'NV =',MESH%NV,', NT =',MESH%NT

     !IF ( ITERATION == 16 ) DEBUG = .TRUE.
     
     IF ( DEBUG ) THEN
        call write_mesh_files( &
             mesh, &
             filetri='optim_unit_mesh/temporary/tri.dat', &
             filexyz='optim_unit_mesh/temporary/xyz.dat' )
        call write_halfedges( &
             mesh, &
             'optim_unit_mesh/temporary/v2h.dat', &
             'optim_unit_mesh/temporary/twin.dat' )
     END IF

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

     
     IF ( .FALSE. ) THEN
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
              !print *,min_edgelen, edgelen, max_edgelen

              if ( edgelen < min_edgelen ) then
                 ! too short -> collapse
                 print *,'collapse edge ',mesh%tri([iv, 1+mod(iv,3)],it)!, &
                 !'len =', edgelen
                 call collapse_edge( &
                      mesh, &
                      [iv,it], &
                      stat )
                 PRINT *,'NV =',MESH%NV,', NT =',MESH%NT
                 IF ( DEBUG ) GOTO 100
                 !cycle outer
              elseif ( edgelen > max_edgelen ) then
                 ! too long -> split
                 print *,'split edge ',mesh%tri([iv, 1+mod(iv,3)],it)!, &
                 !'len =', edgelen
                 call split_edge( &
                      mesh, &
                      [iv,it] )
                 PRINT *,'NV =',MESH%NV,', NT =',MESH%NT
                 IF ( DEBUG ) GOTO 100
                 !cycle outer
              end if
           end do
        end do
     END IF

     !exit

100  IF ( DEBUG ) THEN
        call write_mesh_files( &
             mesh, &
             filetri='optim_unit_mesh/optimized/tri.dat', &
             filexyz='optim_unit_mesh/optimized/xyz.dat' )
        call write_halfedges( &
             mesh, &
             'optim_unit_mesh/optimized/v2h.dat', &
             'optim_unit_mesh/optimized/twin.dat' )
        PAUSE
     END IF

  end do outer



  ! export optimized mesh
  print *,'export optimized mesh...'
  call write_obj_mesh( &
       mesh, &
       'optim_unit_mesh/optimized_mesh.obj' )
  call write_mesh_files( &
       mesh, &
       filetri='optim_unit_mesh/optimized/tri.dat', &
       filexyz='optim_unit_mesh/optimized/xyz.dat' )
  call write_halfedges( &
       mesh, &
       'optim_unit_mesh/optimized/v2h.dat', &
       'optim_unit_mesh/optimized/twin.dat' )


  print *,'done.'

end program optim_unit_mesh
