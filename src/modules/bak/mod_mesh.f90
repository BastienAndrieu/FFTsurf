module mod_mesh

  use mod_math
  
  implicit none

  integer, parameter :: PARAM_xtra_nt = 100
  integer, parameter :: PARAM_xtra_nv = 100
  integer, parameter :: PARAM_xtra_np = 10
  
  type type_path ! trouver nom + explicite...
     integer                    :: nv = 0
     integer, allocatable       :: verts(:)
     real(kind=fp), allocatable :: s(:)
     integer                    :: hyperedge = 0
  end type type_path


  type type_surface_mesh
     integer                      :: nv = 0, nt = 0
     integer, allocatable         :: tri(:,:)
     real(kind=fp), allocatable   :: xyz(:,:)
     real(kind=fp), allocatable   :: uv(:,:,:)
     integer, allocatable         :: ids(:)
     integer, allocatable         :: typ(:)
     integer                      :: npaths = 0
     type(type_path), allocatable :: paths(:)
  end type type_surface_mesh


contains

  subroutine write_inria_mesh( &
       mesh, &
       filename )
    use mod_util
    implicit none
    type(type_surface_mesh), intent(in) :: mesh
    character(*),            intent(in) :: filename
    integer                             :: fid, i, j

    call get_free_unit(fid)
    open(unit=fid, file=filename, action='write')
    
    write (fid,'(A22)') 'MeshVersionFormatted 2'
    write (fid,'(A9)') 'Dimension'
    write (fid,'(I1)') 3
    write (fid,'(A8)') 'Vertices'
    write (fid,'(I0)') mesh%nv
    do j = 1,mesh%nv
       do i = 1,3
          write (fid,'(ES22.15,1x)', advance='no') mesh%xyz(i,j)
       end do
       write (fid,'(I0)') mesh%typ(j)
    end do

    write (fid,'(A9)') 'Triangles'
    write (fid,'(I0)') mesh%nt
    do j = 1,mesh%nt
       write (fid,'(I0,1x,I0,1x,I0,1x,I0)') mesh%tri(1:3,j), 0
    end do
    
    write (fid,'(A3)') 'End'
    close(fid)
    
  end subroutine write_inria_mesh




  


  subroutine write_mesh_files( &
       mesh, &
       filetri, &
       filexyz, &
       fileuv, &
       fileidstyp, &
       filepaths )
    use mod_util
    implicit none
    type(type_surface_mesh), intent(in) :: mesh
    character(*),            intent(in) :: filetri, filexyz, fileuv, fileidstyp, filepaths
    integer                             :: fid, i, j

    call get_free_unit(fid)

    open(unit=fid, file=filetri, action='write')
    do i = 1,mesh%nt
       write (fid,*) mesh%tri(1:3,i)
    end do
    close(fid)

    open(unit=fid, file=filexyz, action='write')
    do i = 1,mesh%nv
       write (fid,*) mesh%xyz(1:3,i)
    end do
    close(fid)

    open(unit=fid, file=fileuv, action='write')
    do i = 1,mesh%nv
       write (fid,*) mesh%uv(1:2,1:2,i)
    end do
    close(fid)

    open(unit=fid, file=fileidstyp, action='write')
    do i = 1,mesh%nv
       write (fid,*) mesh%ids(i), mesh%typ(i)
    end do
    close(fid)

    open(unit=fid, file=filepaths, action='write')
    write (fid,*) mesh%npaths
    do i = 1,mesh%npaths
       write (fid,*) mesh%paths(i)%hyperedge
       write (fid,*) mesh%paths(i)%nv
       do j = 1,mesh%paths(i)%nv
          write (fid,*) mesh%paths(i)%verts(j), mesh%paths(i)%s(j)
       end do
    end do
    close(fid)
    
  end subroutine write_mesh_files








  

  subroutine append_triangles( &
       mesh, &
       tri, &
       nt )
    implicit none
    type(type_surface_mesh), intent(inout) :: mesh
    integer,                 intent(in)    :: nt
    integer,                 intent(in)    :: tri(3,nt)
    integer, allocatable                   :: tmp(:,:)


    if ( allocated(mesh%tri) ) then ! <-------------------------+
       if ( size(mesh%tri,2) < mesh%nt + nt ) then ! <-------+  !
          call move_alloc(from=mesh%tri, to=tmp)             !  !
          allocate(mesh%tri(3,mesh%nt + nt + PARAM_xtra_nt)) !  !
          mesh%tri(1:3,1:mesh%nt) = tmp(1:3,1:mesh%nt)       !  !
          deallocate(tmp)                                    !  !
       end if ! <--------------------------------------------+  !
    else ! -----------------------------------------------------+
       allocate(mesh%tri(3,nt + PARAM_xtra_nt))                 !
    end if ! <--------------------------------------------------+

    mesh%tri(1:3,mesh%nt+1:mesh%nt+nt) = tri(1:3,1:nt)
    mesh%nt = mesh%nt+nt
    
  end subroutine append_triangles



  
  subroutine append_vertices( &
       mesh, &
       xyz, &
       uv, &
       ids, &
       typ, &
       nv )
    implicit none
    type(type_surface_mesh), intent(inout) :: mesh
    integer,                 intent(in)    :: nv
    real(kind=fp),           intent(in)    :: xyz(3,nv)
    real(kind=fp),           intent(in)    :: uv(2,2,nv)
    integer,                 intent(in)    :: ids(nv)
    integer,                 intent(in)    :: typ(nv)
    real(kind=fp), allocatable             :: xyztmp(:,:), uvtmp(:,:,:)
    integer, allocatable                   :: itmp(:)

    if ( allocated(mesh%xyz) ) then ! <----------------------------+
       if ( size(mesh%xyz,2) < mesh%nv + nv ) then ! <----------+  !
          call move_alloc(from=mesh%xyz, to=xyztmp)             !  !
          allocate(mesh%xyz(3,mesh%nv + nv + PARAM_xtra_nv))    !  !
          mesh%xyz(1:3,1:mesh%nv) = xyztmp(1:3,1:mesh%nv)       !  !
          deallocate(xyztmp)                                    !  !
       end if ! <-----------------------------------------------+  !
    else ! --------------------------------------------------------+
       allocate(mesh%xyz(3,nv))                                    !
    end if ! <-----------------------------------------------------+

    if ( allocated(mesh%uv) ) then ! <-----------------------------+
       if ( size(mesh%uv,3) < mesh%nv + nv ) then ! <-----------+  !
          call move_alloc(from=mesh%uv, to=uvtmp)               !  !
          allocate(mesh%uv(2,2,mesh%nv + nv + PARAM_xtra_nv))   !  !
          mesh%uv(1:2,1:2,1:mesh%nv) = uvtmp(1:2,1:2,1:mesh%nv) !  !
          deallocate(uvtmp)                                     !  !
       end if ! <-----------------------------------------------+  !
    else ! --------------------------------------------------------+
       allocate(mesh%uv(2,2,nv))                                   !
    end if ! <-----------------------------------------------------+

    if ( allocated(mesh%ids) ) then ! <----------------------------+
       if ( size(mesh%ids) < mesh%nv + nv ) then ! <------------+  !
          call move_alloc(from=mesh%ids, to=itmp)               !  !
          allocate(mesh%ids(mesh%nv + nv + PARAM_xtra_nv))      !  !
          mesh%ids(1:mesh%nv) = itmp(1:mesh%nv)                 !  !
          deallocate(itmp)                                      !  !
       end if ! <-----------------------------------------------+  !
    else ! --------------------------------------------------------+
       allocate(mesh%ids(nv))                                      !
    end if ! <-----------------------------------------------------+

    if ( allocated(mesh%typ) ) then ! <----------------------------+
       if ( size(mesh%typ) < mesh%nv + nv ) then ! <------------+  !
          call move_alloc(from=mesh%typ, to=itmp)               !  !
          allocate(mesh%typ(mesh%nv + nv + PARAM_xtra_nv))      !  !
          mesh%typ(1:mesh%nv) = itmp(1:mesh%nv)                 !  !
          deallocate(itmp)                                      !  !
       end if ! <-----------------------------------------------+  !
    else ! --------------------------------------------------------+
       allocate(mesh%typ(nv + PARAM_xtra_nv))                      !
    end if ! <-----------------------------------------------------+

    mesh%xyz(1:3,mesh%nv+1:mesh%nv+nv) = xyz(1:3,1:nv)
    mesh%uv(1:2,1:2,mesh%nv+1:mesh%nv+nv) = uv(1:2,1:2,1:nv)
    mesh%ids(mesh%nv+1:mesh%nv+nv) = ids(1:nv)
    mesh%typ(mesh%nv+1:mesh%nv+nv) = typ(1:nv)
    mesh%nv = mesh%nv + nv
    
  end subroutine append_vertices




  subroutine append_path( &
       mesh, &
       path )
    implicit none
    type(type_surface_mesh), intent(inout) :: mesh
    type(type_path),         intent(in)    :: path
    type(type_path), allocatable           :: tmp(:)

    if ( allocated(mesh%paths) ) then ! <-----------------------+
       if ( size(mesh%paths) < mesh%npaths + 1 ) then ! <----+  !
          allocate(tmp(mesh%npaths))                         !  !
          call transfer_paths( &                             !  !
               mesh%paths, &                                 !  !
               tmp, &                                        !  !
               mesh%npaths )                                 !  !
          deallocate(mesh%paths)                             !  !
          allocate(mesh%paths(mesh%npaths + PARAM_xtra_np))  !  !
          call transfer_paths( &                             !  !
               tmp, &                                        !  !
               mesh%paths, &                                 !  !
               mesh%npaths )                                 !  !
          deallocate(tmp)                                    !  !
       end if ! <--------------------------------------------+  !
    else ! -----------------------------------------------------+
       allocate(mesh%paths(PARAM_xtra_np))                      !
    end if ! <--------------------------------------------------+

    mesh%paths(mesh%npaths+1)%nv = path%nv
    if ( allocated(path%verts) ) then
       allocate(mesh%paths(mesh%npaths+1)%verts(path%nv))
       mesh%paths(mesh%npaths+1)%verts(1:path%nv) = path%verts(1:path%nv)
    end if
    if ( allocated(path%s) ) then
       allocate(mesh%paths(mesh%npaths+1)%s(path%nv))
       mesh%paths(mesh%npaths+1)%s(1:path%nv) = path%s(1:path%nv)
    end if
    mesh%paths(mesh%npaths+1)%hyperedge = path%hyperedge
    mesh%npaths = mesh%npaths + 1
    
  end subroutine append_path


  
  subroutine transfer_paths( &
       from, &
       to, &
       n )
    implicit none
    integer,         intent(in)    :: n
    type(type_path), intent(inout) :: from(:)
    type(type_path), intent(inout) :: to(:)
    integer                        :: i

    do i = 1,n
       to(i)%nv = from(i)%nv
       call move_alloc(from(i)%verts, to(i)%verts)
       call move_alloc(from(i)%s, to(i)%s)
       to(i)%hyperedge = from(i)%hyperedge
    end do

  end subroutine transfer_paths

  
end module mod_mesh
