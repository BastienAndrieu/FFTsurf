module mod_mesh_optimization

  use mod_math
  use mod_halfedge
  use mod_mesh

contains

  subroutine swap_edge( &
       mesh, &
       ihedg, &
       stat )
    implicit none
    type(type_surface_mesh), intent(inout) :: mesh
    integer,                 intent(in)    :: ihedg(2)
    integer,                 intent(out)   :: stat
    integer                                :: halfedges(2,2)
    integer                                :: ijk(3,2)
    integer                                :: opposite_vert(2)
    integer                                :: twin_j(2,2)
    integer                                :: ih, jh, iv, it, jt
    
    halfedges(1:2,1) = ihedg
    halfedges(1:2,2) = get_twin(mesh, ihedg)

    do ih = 1,2
       ijk(1,ih) = halfedges(1,ih) ! local index of origin vertex in incident triangle
       it = halfedges(2,ih) ! global index of incident triangle
       if ( it < 1 ) then
          ! cannot swap a boundary edge
          stat = 1
          return
       end if
       ! local indices of next and previous vertices in incident triangle
       ijk(2,ih) = 1 + mod(ijk(1,ih),3)   ! next
       ijk(3,ih) = 1 + mod(ijk(1,ih)+1,3) ! previous
       !
       opposite_vert(ih) = mesh%tri(ijk(3,ih),it)
       twin_j(1:2,ih) = mesh%twin(1:2,ijk(2,ih),it)
    end do
    
    stat = 0

    do ih = 1,2
       it = halfedges(2,ih)
       jh = 1 + mod(ih,2)
       jt = halfedges(2,jh)
       ! change tri (f2v)
       mesh%tri(ijk(2,ih),it) = opposite_vert(jh)
       ! change 'interior' twins
       mesh%twin(1:2,ijk(2,ih),it) = [ijk(2,ih), jt]
       ! change 'diagonal' twins
       mesh%twin(1:2,ijk(1,ih),it) = twin_j(1:2,jh)
       if ( all(twin_j(1:2,jh) > 0) ) then
          mesh%twin(1:2,twin_j(1,jh),twin_j(2,jh)) = [ijk(1,ih),it]
       end if
       ! change v2h if necessary
       iv = mesh%tri(ijk(1,ih),it)
       if ( mesh%v2h(2,iv) == halfedges(2,jh) ) mesh%v2h(1:2,iv) = [ijk(1,ih),it]
    end do

  end subroutine swap_edge






  subroutine split_edge( &
       mesh, &
       ihedg )
    use mod_util
    implicit none
    LOGICAL :: DEBUG = .FALSE.
    type(type_surface_mesh), intent(inout) :: mesh
    integer,                 intent(in)    :: ihedg(2)
    integer                                :: halfedges(2,2)
    integer, dimension(2)                  :: verts, faces, faces_new
    integer                                :: ijk(3,2)
    integer                                :: twin_j(2,2)
    integer                                :: n_new_faces
    real(kind=fp)                          :: uv_new(2,2,1)
    real(kind=fp)                          :: xyz_new(3,1)
    integer                                :: typ_new(1)
    integer                                :: ids_new(1)
    integer                                :: tri_new(3,2)
    integer                                :: ihf_new(2)
    integer, allocatable                   :: twin_tmp(:,:,:)
    integer                                :: nv_tmp
    integer                                :: ih, jh, it, iv
    
    verts(1) = get_orig(mesh, ihedg)
    verts(2) = get_dest(mesh, ihedg)

    halfedges(1:2,1) = ihedg
    halfedges(1:2,2) = get_twin(mesh, ihedg)
    
    do ih = 1,2
       ijk(1,ih) = halfedges(1,ih) ! local index of origin vertex in incident triangle
       it = halfedges(2,ih) ! global index of incident triangle
       faces(ih) = it
       ! local indices of next and previous vertices in incident triangle
       ijk(2,ih) = 1 + mod(ijk(1,ih),3)   ! next
       ijk(3,ih) = 1 + mod(ijk(1,ih)+1,3) ! previous
       !
       if ( it > 0 ) twin_j(1:2,ih) = mesh%twin(1:2,ijk(2,ih),it)
    end do
    n_new_faces = count(faces > 0)

    ! insert new triangles
    do ih = 1,n_new_faces
       tri_new(1,ih) = mesh%nv + 1
       tri_new(2:3,ih) = mesh%tri(ijk(2:3,ih),faces(ih))
       ihf_new(ih) = mesh%ihf(faces(ih))
       faces_new(ih) = mesh%nt + ih
    end do
    
    call append_triangles( &
         mesh, &
         tri_new(1:3,1:n_new_faces), &
         ihf_new(1:n_new_faces), &
         n_new_faces )

    ! insert new vertex
    ! uv, typ, ids***
    xyz_new(1:3,1) = 0.5_fp*(mesh%xyz(1:3,verts(1)) + mesh%xyz(1:3,verts(2))) !***
    call append_vertices( &
         mesh, &
         xyz_new, &
         uv_new, &
         ids_new, &
         typ_new, &
         1 )

    ! edit tri
    do ih = 1,n_new_faces
       mesh%tri(ijk(2,ih),faces(ih)) = mesh%nv
    end do

    ! edit twins
    if ( size(mesh%twin,3) < mesh%nt ) then
       call move_alloc(from=mesh%twin, to=twin_tmp)
       allocate(mesh%twin(2,3,mesh%nt + MESH_xtra_nt))
       mesh%twin(1:2,1:3,1:size(twin_tmp,3)) = twin_tmp(1:2,1:3,1:size(twin_tmp,3))
       deallocate(twin_tmp)
    end if
    do ih = 1,n_new_faces
       mesh%twin(1:2,ijk(2,ih),faces(ih)) = [3, faces_new(ih)]
       mesh%twin(1:2,3,faces_new(ih)) = [ijk(2,ih), faces(ih)]
       mesh%twin(1:2,2,faces_new(ih)) = twin_j(1:2,ih)
       if ( all(twin_j(1:2,ih) > 0) ) then
          mesh%twin(1:2,twin_j(1,ih),twin_j(2,ih)) = [2,faces_new(ih)]
       end if
       if ( n_new_faces > 1 ) then
          jh = 1 + mod(ih,2)
          mesh%twin(1:2,ijk(1,ih),faces(ih)) = [1, faces_new(jh)]
          mesh%twin(1:2,1,faces_new(ih)) = [ijk(1,jh), faces(jh)]    
       else
          mesh%twin(1:2,1,faces_new(ih)) = [0,0]
       end if
    end do

    ! edit v2h
    nv_tmp = mesh%nv - 1
    iv = nv_tmp
    call insert_column_after( &
      mesh%v2h, &
      2, &
      nv_tmp, &
      [1, faces_new(1)], &
      iv )
    do ih = 1,n_new_faces
       iv = verts(1+mod(ih,2))
       if ( mesh%v2h(2,iv) == faces(ih) ) then
          IF ( DEBUG ) PRINT *, IV, ', ', mesh%v2h(1:2,iv), '--->', [2, faces_new(ih)]
          mesh%v2h(1:2,iv) = [2, faces_new(ih)]
       end if
    end do
    
  end subroutine split_edge






  subroutine collapse_edge( &
       mesh, &
       ihedg, &
       stat )
    implicit none
    LOGICAL :: DEBUG = .FALSE.
    type(type_surface_mesh), intent(inout) :: mesh
    integer,                 intent(in)    :: ihedg(2)
    integer,                 intent(out)   :: stat
    integer                                :: halfedges(2,2)
    integer, dimension(2)                  :: verts, faces, opposite_vert
    integer                                :: twin_jk(2,2,2)
    integer                                :: i_new_vert, i_removed_vert
    integer                                :: vert_id_new(mesh%nv)
    integer                                :: face_id_new(mesh%nt)
    integer                                :: ih, iv, jv, kv, it, i

    verts(1) = get_orig(mesh, ihedg)
    verts(2) = get_dest(mesh, ihedg)

    ! prevent degenerate cases ---------------------------------------------+
    if ( is_boundary_vertex(mesh, verts(1)) .and. &                         !
         is_boundary_vertex(mesh, verts(2)) .and. &                         !
         .not.is_boundary_halfedge(mesh, ihedg) ) then                      !
       ! the edge collapse would lead to a non-manifold or degenerate mesh  !
       PRINT *,'edge collapse would lead to a degenerate mesh'
       stat = 1                                                             !
       return                                                               !
    end if                                                                  !
    !                                                                       !
    if ( is_boundary_halfedge(mesh, ihedg) .and. &                          !
         is_boundary_halfedge(mesh, get_prev(ihedg)) .and. &                !
         is_boundary_halfedge(mesh, get_next(ihedg)) ) then                 !
       ! do not collapse an edge of an isolated triangle                    !
       PRINT *,'cannot collapse edge of an isolated triangle'
       stat = 1                                                             !
       return                                                               !
    end if                                                                  !
    ! ----------------------------------------------------------------------+

    stat = 0
    
    halfedges(1:2,1) = ihedg
    halfedges(1:2,2) = get_twin(mesh, ihedg)

    if ( is_boundary_vertex(mesh, verts(1)) .and. &
         .not.is_boundary_vertex(mesh, verts(2)) ) then
       ! permute configuration
       IF ( DEBUG ) PRINT *,'PERMUTE CONFIGURATION'
       halfedges(1:2,1:2) = halfedges(1:2,[2,1])
       verts(1:2) = verts([2,1])
    end if
    
    do ih = 1,2
       iv = halfedges(1,ih) ! local index of origin vertex in incident triangle
       it = halfedges(2,ih) ! global index of incident triangle
       faces(ih) = it
       if ( it > 0 ) then
          jv = 1 + mod(iv,3) ! local index of next vertex in incident triangle
          kv = 1 + mod(iv+1,3) ! local index of previous vertex in incident triangle
          ! store current twins of previous and next halfedges in current triangle
          twin_jk(1:2,1,ih) = mesh%twin(1:2,jv,it)
          twin_jk(1:2,2,ih) = mesh%twin(1:2,kv,it)
          !
          opposite_vert(ih) = mesh%tri(kv,it) ! global index of opposite vertex
       end if
    end do
    IF ( DEBUG ) PRINT *,'FACES =',FACES

    ! new vertex resulting from collapse
    i_new_vert = minval(verts)
    i_removed_vert = maxval(verts)
    ! uv, typ, ids, xyz***
    if ( is_boundary_vertex(mesh, verts(2)) .and. &
         .not.is_boundary_vertex(mesh, verts(1)) ) then
       mesh%xyz(1:3,i_new_vert) = mesh%xyz(1:3,verts(2))
    else
       mesh%xyz(1:3,i_new_vert) = 0.5_fp*(mesh%xyz(1:3,verts(1)) + mesh%xyz(1:3,verts(2))) !***
    end if ! *******

    ! change v2h of new vertex
    if ( is_boundary_vertex(mesh, verts(2)) .and. &
         mesh%v2h(2,verts(2)) /= faces(1) ) then
       mesh%v2h(1:2,i_new_vert) = mesh%v2h(1:2,verts(2))
    else
       mesh%v2h(1:2,i_new_vert) = twin_jk(1:2,2,1)
    end if
    IF ( DEBUG ) PRINT *,'NEW VERT:',i_new_vert,', V2H:', mesh%v2h(1:2,i_new_vert)

    ! change v2h of opposite vertices
    do ih = 1,2
       if ( faces(ih) < 1 ) cycle
       iv = opposite_vert(ih)
       if ( mesh%v2h(2,iv) == faces(ih) ) then
          IF ( DEBUG ) PRINT *,'OPP. VERT:',IV,', V2H:',mesh%v2h(1:2,iv),' -->',twin_jk(1:2,1,ih)
          mesh%v2h(1:2,iv) = twin_jk(1:2,1,ih)
       end if
    end do
    
    ! change twins
    do ih = 1,2
       if ( faces(ih) < 1 ) cycle
       do iv = 1,2
          jv = 1 + mod(iv,2)
          if ( all(twin_jk(1:2,iv,ih) > 0) ) then
             mesh%twin(1:2,twin_jk(1,iv,ih),twin_jk(2,iv,ih)) = twin_jk(1:2,jv,ih)
          end if
       end do
    end do

    ! new face indices
    if ( faces(1) > faces(2) ) faces = faces([2,1]) ! sort by ascending face index
    if ( faces(1) > 0 ) then
       face_id_new(1:faces(1)-1) = [(i, i=1,faces(1)-1)]
       face_id_new(faces(1)+1:faces(2)-1) = [(i, i=faces(1),faces(2)-2)]
       face_id_new(faces(2)+1:mesh%nt) = [(i, i=faces(2)-1,mesh%nt-2)]
    else
       face_id_new(1:faces(2)-1) = [(i, i=1,faces(2)-1)]
       face_id_new(faces(2)+1:mesh%nt) = [(i, i=faces(2),mesh%nt-1)]
    end if
    
    ! remove triangles incident to collapsed edge
    if ( faces(1) > 0 ) then
       mesh%tri(1:3,faces(1):faces(2)-2) = mesh%tri(1:3,faces(1)+1:faces(2)-1)
       mesh%tri(1:3,faces(2)-1:mesh%nt-2) = mesh%tri(1:3,faces(2)+1:mesh%nt)
       mesh%twin(1:2,1:3,faces(1):faces(2)-2) = mesh%twin(1:2,1:3,faces(1)+1:faces(2)-1)
       mesh%twin(1:2,1:3,faces(2)-1:mesh%nt-2) = mesh%twin(1:2,1:3,faces(2)+1:mesh%nt)
       mesh%ihf(faces(1):faces(2)-2) = mesh%ihf(faces(1)+1:faces(2)-1)
       mesh%ihf(faces(2)-1:mesh%nt-2) = mesh%ihf(faces(2)+1:mesh%nt)
       mesh%nt = mesh%nt - 2
    else
       mesh%tri(1:3,faces(2):mesh%nt-1) = mesh%tri(1:3,faces(2)+1:mesh%nt)
       mesh%twin(1:2,1:3,faces(2):mesh%nt-1) = mesh%twin(1:2,1:3,faces(2)+1:mesh%nt)
       mesh%ihf(faces(2):mesh%nt-1) = mesh%ihf(faces(2)+1:mesh%nt)
       mesh%nt = mesh%nt - 1
    end if

    ! new vertex indices
    vert_id_new(1:i_removed_vert-1) = [(i, i=1,i_removed_vert-1)]
    vert_id_new(i_removed_vert) = i_new_vert
    vert_id_new(i_removed_vert+1:mesh%nv) = [(i, i=i_removed_vert,mesh%nv-1)]

    ! apply re-ordering
    do it = 1,mesh%nt
       mesh%tri(1:3,it) = vert_id_new(mesh%tri(1:3,it))
       do iv = 1,3
          if ( mesh%twin(2,iv,it) > 0 ) then
             mesh%twin(2,iv,it) = face_id_new(mesh%twin(2,iv,it))
          end if
       end do          
    end do

    ! remove vertex
    ! uv, typ, ids***
    mesh%xyz(1:3,i_removed_vert:mesh%nv-1) = mesh%xyz(1:3,i_removed_vert+1:mesh%nv)
    mesh%v2h(1:2,i_removed_vert:mesh%nv-1) = mesh%v2h(1:2,i_removed_vert+1:mesh%nv)
    mesh%nv = mesh%nv - 1
    ! apply re-ordering
    do iv = 1,mesh%nv
       mesh%v2h(2,iv) = face_id_new(mesh%v2h(2,iv))
    end do
    ! possibly remove from path ...***
  
  end subroutine collapse_edge

end module mod_mesh_optimization
