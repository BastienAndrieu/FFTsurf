module mod_mesh

  use mod_math
  
  implicit none

  integer, parameter :: MESH_xtra_nt = 100
  integer, parameter :: MESH_xtra_nv = 100
  integer, parameter :: MESH_xtra_np = 10
  
  type type_path ! trouver nom + explicite...
     integer                    :: nv = 0
     integer, allocatable       :: verts(:)
     real(kind=fp), allocatable :: s(:)
     integer                    :: hyperedge = 0
  end type type_path


  type type_surface_mesh
     integer                      :: nv = 0, nt = 0
     integer, allocatable         :: v2h(:,:), twin(:,:,:)
     integer, allocatable         :: tri(:,:) ! f2v...
     real(kind=fp), allocatable   :: xyz(:,:)
     real(kind=fp), allocatable   :: uv(:,:,:)
     integer, allocatable         :: ids(:)
     integer, allocatable         :: typ(:)
     integer                      :: npaths = 0
     type(type_path), allocatable :: paths(:)
  end type type_surface_mesh


contains
  

  subroutine make_halfedges( &
       mesh )
    implicit none
    type(type_surface_mesh), intent(inout) :: mesh
    integer, allocatable                   :: h(:,:)
    integer                                :: pairv(2)
    integer                                :: ifa, ied, ih

    allocate(mesh%twin(2,3,mesh%nt))
    mesh%twin(:,:,:) = 0
    
    allocate(h(2,mesh%nv*(mesh%nv - 1) / 2))
    h(:,:) = 0
    do ifa = 1,mesh%nt
       do ied = 1,3
          pairv = mesh%tri([ied, 1+mod(ied,3)],ifa)
          ih = hash_integer_pair(pairv)

          if ( h(1,ih) > 0 ) then
             ! collision -> we found the twin
             mesh%twin(1:2,ied,ifa) = h(1:2,ih)
             mesh%twin(1:2,h(1,ih),h(2,ih)) = [ied,ifa]
          else
             h(1:2,ih) = [ied,ifa]
          end if
       end do
    end do
    deallocate(h)

    ! vertex -> outgoing half-edge incidences
    allocate(mesh%v2h(2,mesh%nv))
    mesh%v2h(:,:) = 0
    do ifa = 1,mesh%nt
       do ied = 1,3
          if ( mesh%twin(1,ied,ifa) > 0 ) then
             if ( mesh%v2h(1,mesh%tri(ied,ifa)) > 0 ) cycle
          end if

          mesh%v2h(1:2,mesh%tri(ied,ifa)) = [ied,ifa]
       end do
    end do
    
  end subroutine make_halfedges


  
  function hash_integer_pair( &
       pq )
    implicit none
    integer, intent(in) :: pq(2)
    integer             :: hash_integer_pair
    integer             :: r, s

    r = min(pq(1), pq(2))
    if ( pq(1) == pq(2) .or. r < 1 ) then
       hash_integer_pair = 0
    else
       s = max(pq(1), pq(2))
       hash_integer_pair = r + (s-1)*(s-2) / 2
    end if

  end function hash_integer_pair


  

  

  subroutine write_tecplot_mesh( &
       mesh, &
       filename, &
       zonename )
    use mod_util
    implicit none
    type(type_surface_mesh), intent(in) :: mesh
    character(*),            intent(in) :: filename
    character(*),            intent(in) :: zonename
    integer                             :: fid, i, j

    call get_free_unit(fid)
    open(unit=fid, file=filename, action='write')

    write (fid,*) 'VARIABLES = "X" "Y" "Z"'

    write (fid,*) 'ZONE T="' // zonename // '"'
    write (fid,'(A2,I5,A3,I5)') 'N=',mesh%nv,' E=',mesh%nt
    write (fid,*) 'ZONETYPE=FETriangle'
    write (fid,*) 'DATAPACKING=BLOCK'

    write (fid,*) 'VARLOCATION = ([1-3]=NODAL)'
    write (fid,*) 'DT=(DOUBLE DOUBLE DOUBLE)'

    do i = 1,3
       write (fid,*) ''
       do j = 1,mesh%nv
          write (fid,'(ES22.15)') mesh%xyz(i,j)
       end do
    end do

    write (fid,*) ''

    do j = 1,mesh%nt
       write (fid,'(3I7)') mesh%tri(:,j)
    end do

    close(fid)
    
  end subroutine write_tecplot_mesh



  

  

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




  subroutine write_obj_mesh( &
       mesh, &
       filename )
    use mod_util
    implicit none
    character(7), parameter :: fmt = 'ES13.6'
    type(type_surface_mesh), intent(in) :: mesh
    character(*),            intent(in) :: filename
    integer                             :: fid, i, j

    call get_free_unit(fid)
    open(unit=fid, file=filename, action='write')
    
    do j = 1,mesh%nv
       write (fid,'(A1,1x)', advance='no') 'v'
       do i = 1,3
          write (fid,'('//fmt//',1x)', advance='no') mesh%xyz(i,j)
       end do
       write (fid,*)
    end do

    do j = 1,mesh%nt
       write (fid,'(A1,1x,I0,1x,I0,1x,I0)') 'f', mesh%tri(1:3,j)
    end do
    
    close(fid)

  end subroutine write_obj_mesh
  


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
          allocate(mesh%tri(3,mesh%nt + nt + MESH_xtra_nt))  !  !
          mesh%tri(1:3,1:mesh%nt) = tmp(1:3,1:mesh%nt)       !  !
          deallocate(tmp)                                    !  !
       end if ! <--------------------------------------------+  !
    else ! -----------------------------------------------------+
       allocate(mesh%tri(3,nt + MESH_xtra_nt))                  !
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
          allocate(mesh%xyz(3,mesh%nv + nv + MESH_xtra_nv))     !  !
          mesh%xyz(1:3,1:mesh%nv) = xyztmp(1:3,1:mesh%nv)       !  !
          deallocate(xyztmp)                                    !  !
       end if ! <-----------------------------------------------+  !
    else ! --------------------------------------------------------+
       allocate(mesh%xyz(3,nv))                                    !
    end if ! <-----------------------------------------------------+

    if ( allocated(mesh%uv) ) then ! <-----------------------------+
       if ( size(mesh%uv,3) < mesh%nv + nv ) then ! <-----------+  !
          call move_alloc(from=mesh%uv, to=uvtmp)               !  !
          allocate(mesh%uv(2,2,mesh%nv + nv + MESH_xtra_nv))    !  !
          mesh%uv(1:2,1:2,1:mesh%nv) = uvtmp(1:2,1:2,1:mesh%nv) !  !
          deallocate(uvtmp)                                     !  !
       end if ! <-----------------------------------------------+  !
    else ! --------------------------------------------------------+
       allocate(mesh%uv(2,2,nv))                                   !
    end if ! <-----------------------------------------------------+

    if ( allocated(mesh%ids) ) then ! <----------------------------+
       if ( size(mesh%ids) < mesh%nv + nv ) then ! <------------+  !
          call move_alloc(from=mesh%ids, to=itmp)               !  !
          allocate(mesh%ids(mesh%nv + nv + MESH_xtra_nv))       !  !
          mesh%ids(1:mesh%nv) = itmp(1:mesh%nv)                 !  !
          deallocate(itmp)                                      !  !
       end if ! <-----------------------------------------------+  !
    else ! --------------------------------------------------------+
       allocate(mesh%ids(nv))                                      !
    end if ! <-----------------------------------------------------+

    if ( allocated(mesh%typ) ) then ! <----------------------------+
       if ( size(mesh%typ) < mesh%nv + nv ) then ! <------------+  !
          call move_alloc(from=mesh%typ, to=itmp)               !  !
          allocate(mesh%typ(mesh%nv + nv + MESH_xtra_nv))       !  !
          mesh%typ(1:mesh%nv) = itmp(1:mesh%nv)                 !  !
          deallocate(itmp)                                      !  !
       end if ! <-----------------------------------------------+  !
    else ! --------------------------------------------------------+
       allocate(mesh%typ(nv + MESH_xtra_nv))                       !
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
          allocate(mesh%paths(mesh%npaths + MESH_xtra_np))   !  !
          call transfer_paths( &                             !  !
               tmp, &                                        !  !
               mesh%paths, &                                 !  !
               mesh%npaths )                                 !  !
          deallocate(tmp)                                    !  !
       end if ! <--------------------------------------------+  !
    else ! -----------------------------------------------------+
       allocate(mesh%paths(MESH_xtra_np))                       !
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
