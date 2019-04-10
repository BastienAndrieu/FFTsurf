program optimmesh_topochange

  use mod_util
  use mod_mesh
  use mod_mesh_optimization

  implicit none
  
  character(100)             :: argstr
  type(type_surface_mesh)    :: mesh
  integer                    :: iface, ivert
  integer                    :: tri(3)
  integer, allocatable       :: tritmp(:,:)
  real(kind=fp)              :: xyz(3)
  real(kind=fp), allocatable :: xyztmp(:,:)
  integer                    :: fid, err, i
  integer                    :: stat
  
  ! load mesh data
  print *,'load triangles...'
  call get_free_unit(fid)
  allocate(mesh%tri(3,200))
  open(unit=fid, file='optimmesh/test_mesh_tri.dat', action='read')
  do
     read (fid,*,iostat=err) tri
     if ( err /=0 ) exit
     if ( mesh%nt + 1 > size(mesh%tri,2) ) then
        call move_alloc(from=mesh%tri, to=tritmp)
        allocate(mesh%tri(3,mesh%nt + 2000))
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
  open(unit=fid, file='optimmesh/test_mesh_xyz.dat', action='read')
  do
     read (fid,*,iostat=err) xyz
     if ( err /=0 ) exit
     if ( mesh%nv + 1 > size(mesh%xyz,2) ) then
        call move_alloc(from=mesh%xyz, to=xyztmp)
        allocate(mesh%xyz(3,mesh%nv + 1000))
        mesh%xyz(1:3,1:mesh%nv) = xyztmp(1:3,1:mesh%nv)
        deallocate(xyztmp)
     end if
     mesh%nv = mesh%nv + 1
     mesh%xyz(1:3,mesh%nv) = xyz(1:3)
  end do
  close(fid)
  allocate(mesh%uv(2,2,mesh%nv), mesh%typ(mesh%nv), mesh%ids(mesh%nv))


  ! make halfedge data structure
  print *,'make halfedge data structure...'
  call make_halfedges(mesh)


  ! export initial mesh
  print *,'export initial mesh...'
  call write_obj_mesh( &
       mesh, &
       'optimmesh/initial_mesh.obj' )

  call write_mesh_files( &
       mesh, &
       filetri='optimmesh/initial/tri.dat', &
       filexyz='optimmesh/initial/xyz.dat' )

  open(unit=fid, file='optimmesh/initial/v2h.dat', action='write')
  do i = 1,mesh%nv
     write (fid,*) mesh%v2h(:,i)
  end do
  close(fid)
  open(unit=fid, file='optimmesh/initial/twin.dat', action='write')
  do i = 1,mesh%nt
     write (fid,*) mesh%twin(:,:,i)
  end do
  close(fid)


  !! Read argument
  if ( command_argument_count() < 2 ) then
     STOP 'user must provide halfedge ID'
  else
     call get_command_argument(1, argstr)
     read (argstr,*) ivert
     call get_command_argument(2, argstr)
     read (argstr,*) iface
  end if

  ! collapse an edge
  print *,'collapse edge', [ivert,iface]
  call collapse_edge( &
       mesh, &
       [ivert,iface], &
       stat )
  !IF ( .true. ) THEN
  !   outer1 : do iface = 1,mesh%nt
  !      do ivert = 1,3
  !         if ( mesh%twin(2,ivert,iface) < 1 ) then
  !            print *,'collapse edge', [ivert,iface]!get_prev([ivert,iface])
  !            call collapse_edge( &
  !                 mesh, &
  !                 [ivert,iface] )!get_prev([ivert,iface]) )
  !            exit outer1
  !         end if
  !      end do
  !   end do outer1
  !END IF

  ! swap an edge
  IF ( .false. ) THEN
     outer2 : do iface = 1,mesh%nt
        do ivert = 1,3
           if ( mesh%twin(2,ivert,iface) > 0 ) then
              print *,'swap edge', ivert, iface
              call swap_edge( &
                   mesh, &
                   [ivert,iface], &
                   stat )
              exit outer2
           end if
        end do
     end do outer2
  END IF
  
  ! split an edge
  IF ( .false. ) THEN
     outer3 : do iface = 1,mesh%nt
        do ivert = 1,3
           if ( mesh%twin(2,ivert,iface) < 1 ) then
              print *,'split edge', ivert, iface
              call split_edge( &
                   mesh, &
                   [ivert,iface] )
              exit outer3
           end if
        end do
     end do outer3
  END IF

  
  ! export optimized mesh
  print *,'export optimized mesh...'
  call write_obj_mesh( &
       mesh, &
       'optimmesh/optimized_mesh.obj' )

  call write_mesh_files( &
       mesh, &
       filetri='optimmesh/optimized/tri.dat', &
       filexyz='optimmesh/optimized/xyz.dat' )

  open(unit=fid, file='optimmesh/optimized/v2h.dat', action='write')
  do i = 1,mesh%nv
     write (fid,*) mesh%v2h(:,i)
  end do
  close(fid)
  open(unit=fid, file='optimmesh/optimized/twin.dat', action='write')
  do i = 1,mesh%nt
     write (fid,*) mesh%twin(:,:,i)
  end do
  close(fid)

  print *,'done.'
  
end program optimmesh_topochange
