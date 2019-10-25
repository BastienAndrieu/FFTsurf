module mod_mesh_optimization

  use mod_math
  use mod_halfedge
  use mod_mesh

  implicit none

contains

  subroutine flip_edge( &
       mesh, &
       ihedg, &
       stat )
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! TO DO: CHECK IF FLIPPING EDGE DOESN'T MAKE INVERTED NEW TRIANGLES
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    LOGICAL :: DEBUG = .true.
    real(kind=fp), parameter               :: EPSflip = 0.05_fp
    type(type_surface_mesh), intent(inout) :: mesh
    integer,                 intent(in)    :: ihedg(2)
    integer,                 intent(out)   :: stat
    integer                                :: halfedges(2,2)
    integer                                :: ijk(3,2)
    integer                                :: opposite_vert(2)
    integer                                :: twin_j(2,2)
    integer                                :: ih, jh, iv, it, jt
    real(kind=fp), dimension(3)            :: aBefore, aAfter

    IF ( DEBUG ) THEN
       PRINT *, '--- FLIP EDGE ---'
       PRINT *, '   EDGE #', IHEDG(1), ' OF FACE #', IHEDG(2)
       PRINT *, '   VERTS #', GET_ORIG(MESH, IHEDG), GET_DEST(MESH, IHEDG)
       PRINT *, '   XYZ_MID =', 0.5_FP*(MESH%XYZ(:,GET_ORIG(MESH, IHEDG)) + &
            MESH%XYZ(:,GET_DEST(MESH, IHEDG)))
    END IF

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

    ! check if edge flip produces (near-)degenerate triangles
    do ih = 1,2
      it = halfedges(2,ih)
      aBefore = cross( &
         mesh%xyz(1:3,mesh%tri(2,it)) - mesh%xyz(1:3,mesh%tri(1,it)), &
         mesh%xyz(1:3,mesh%tri(3,it)) - mesh%xyz(1:3,mesh%tri(1,it)) )
      !
      jh = 1 + mod(ih,2)
      aAfter = cross( &
         mesh%xyz(1:3,opposite_vert(jh)) - mesh%xyz(1:3,mesh%tri(ijk(1,ih),it)), &
         mesh%xyz(1:3,mesh%tri(ijk(3,ih),it)) - mesh%xyz(1:3,mesh%tri(ijk(1,ih),it)) )
      !
      if ( sign(1._fp, dot_product(aAfter, aBefore)) * &
         sqrt(sum(aAfter**2)/sum(aBefore**2)) < EPSflip ) then
         stat = 2
         return
      end if
    end do
    stat = 0

    do ih = 1,2
       it = halfedges(2,ih)
       jh = 1 + mod(ih,2)
       jt = halfedges(2,jh)
       ! change tri (f2v)
       mesh%tri(ijk(2,ih),it) = opposite_vert(jh)
       ! change 'interior' twins
       mesh%twin(1:2,ijk(2,ih),it) = [ijk(2,jh), jt]![ijk(2,ih), jt]
       ! change 'diagonal' twins
       mesh%twin(1:2,ijk(1,ih),it) = twin_j(1:2,jh)
       if ( all(twin_j(1:2,jh) > 0) ) then
          mesh%twin(1:2,twin_j(1,jh),twin_j(2,jh)) = [ijk(1,ih),it]
       end if
       ! change v2h if necessary
       iv = mesh%tri(ijk(1,ih),it)
       if ( mesh%v2h(2,iv) == halfedges(2,jh) ) mesh%v2h(1:2,iv) = [ijk(1,ih),it]
    end do

  end subroutine flip_edge






  subroutine split_edge( &
       mesh, &
       brep, &
       hypg, &
       ihedg, &
       resulting_vertex )
    use mod_util
    use mod_types_brep
    use mod_hypergraph
    implicit none
    LOGICAL :: DEBUG = .true.
    !LOGICAL :: DEBUG_PROJ = .FALSE.
    type(type_surface_mesh), intent(inout) :: mesh
    type(type_brep),         intent(in)    :: brep
    type(type_hypergraph),   intent(in)    :: hypg
    integer,                 intent(in)    :: ihedg(2)
    integer, optional,       intent(out)   :: resulting_vertex
    integer                                :: halfedges(2,2)
    integer, dimension(2)                  :: verts, faces, faces_new
    integer                                :: ijk(3,2)
    integer                                :: twin_j(2,2)
    integer                                :: n_new_faces
    real(kind=fp)                          :: uv_new(2,2)
    real(kind=fp)                          :: xyz_new(3)!,1)
    integer                                :: typ_new!(1)
    integer                                :: ids_new!(1)
    integer                                :: tri_new(3,2)
    integer                                :: ihf_new(2)
    integer, allocatable                   :: twin_tmp(:,:,:)
    integer                                :: nv_tmp
    integer                                :: ih, jh, it, iv
    integer                                :: stat_proj
    integer                                :: order(2), cavity(4)

    verts(1) = get_orig(mesh, ihedg)
    verts(2) = get_dest(mesh, ihedg)

    IF ( DEBUG ) THEN
       PRINT *, ''
       PRINT *, '--- SPLIT EDGE ---'
       PRINT *, '   BEFORE SPLIT, MESH%NV, MESH%NT =', MESH%NV, MESH%NT
       PRINT *, '   EDGE #', IHEDG(1), ' OF FACE #', IHEDG(2)
       PRINT *, '   VERTS #', VERTS
       PRINT *, '   XYZ_MID =', 0.5_FP*(MESH%XYZ(:,VERTS(1)) + MESH%XYZ(:,VERTS(2)))
    END IF

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
    IF ( DEBUG ) THEN
      PRINT *,'FACES =',FACES
      DO IT = 1,2
         PRINT *, FACES(IT), '  VERTS:', MESH%TRI(1:3,FACES(IT))
      END DO
    END IF

    ! insert new triangles
    IF ( DEBUG ) PRINT *,'NEW FACE(S):'
    do ih = 1,n_new_faces
       tri_new(1,ih) = mesh%nv + 1
       tri_new(2:3,ih) = mesh%tri(ijk(2:3,ih),faces(ih))
       ihf_new(ih) = mesh%ihf(faces(ih))
       faces_new(ih) = mesh%nt + ih
       IF ( DEBUG ) PRINT *,faces_new(ih), ', VERTS:', tri_new(1:3,IH)
    end do

    call append_triangles( &
         mesh, &
         tri_new(1:3,1:n_new_faces), &
         ihf_new(1:n_new_faces), &
         n_new_faces )

    ! insert new vertex
    ! 1) sort verts in descending typ order
    !if ( mesh%typ(verts(1)) <  mesh%typ(verts(2))) verts = verts([2,1])
    if ( mesh%typ(verts(1)) <  mesh%typ(verts(2))) then
      order = [2,1]
    else
      order = [1,2]
    end if
    !dxyz = 0.5_fp*(mesh%xyz(1:3,verts(2)) - mesh%xyz(1:3,verts(1)))
    typ_new = mesh%typ(verts(order(1)))

    ! 2) find midpoint
    select case ( mesh%typ(verts(order(1))) ) 
    case (1)
      call get_midpoint_hyperedge( &
         mesh, &
         verts(order), &
         brep, &
         hypg, &
         stat_proj, &
         ids_new, &
         uv_new, &
         xyz_new )
      if ( stat_proj > 0 ) then
         print *,'split_edge: failed to project on hyperedge'
         PAUSE
      end if
      !
    case (2)
      IF ( .TRUE. ) THEN
         cavity(1:2) = verts(order(1:2))
         cavity(3) = mesh%tri(ijk(3,1),faces(1))
         cavity(4) = mesh%tri(ijk(3,2),faces(2))
         call get_centroid_hyperface( &
            mesh, &
            cavity, &
            4, &
            brep, &
            stat_proj, &
            ids_new, &
            uv_new(1:2,1), &
            xyz_new )
      ELSE
         call get_midpoint_hyperface( &
            mesh, &
            verts(order), &
            brep, &
            stat_proj, &
            ids_new, &
            uv_new(1:2,1), &
            xyz_new )
      END IF
      if ( stat_proj > 0 ) then
         print *,'split_edge: failed to project on hyperface'
         PAUSE
      end if
      !
    end select

    ! uv, typ, ids***
    !xyz_new(1:3,1) = 0.5_fp*(mesh%xyz(1:3,verts(1)) + mesh%xyz(1:3,verts(2))) !***
    call append_vertices( &
         mesh, &
         xyz_new, &
         uv_new, &
         [ids_new], &
         [typ_new], &
         1 )
    if ( present(resulting_vertex) ) resulting_vertex = mesh%nv
    if ( allocated(mesh%hTargetV) ) mesh%hTargetV(mesh%nv) = 0.5_fp*sum(mesh%hTargetV(verts))

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
    !PRINT *,'BEFORE INSERT COLUMN, SIZE V2H =', SIZE(MESH%V2H,2)
    call insert_column_after( &
         mesh%v2h, &
         2, &
         nv_tmp, &
         [1, faces_new(1)], &
         iv )
    !PRINT *, 'NEW VERT, V2H =', MESH%V2H(1:2,MESH%NV), ' (MESH%NT =', MESH%NT, ')'
    IF ( NV_TMP /= MESH%NV ) THEN
      PRINT *,'split_edge: ERROR INSERT COLUMN V2H'
      PAUSE
    END IF
    
    !IF ( MESH%V2H(2,MESH%NV) > MESH%NT ) PAUSE
    do ih = 1,n_new_faces
       iv = verts(1+mod(ih,2))
       if ( mesh%v2h(2,iv) == faces(ih) ) then
          IF ( DEBUG ) PRINT *, IV, ', ', mesh%v2h(1:2,iv), '--->', [2, faces_new(ih)]
          mesh%v2h(1:2,iv) = [2, faces_new(ih)]
       end if
    end do

    !DO IV = 1,MESH%NV
    !  IF ( MESH%V2H(2,IV) > MESH%NT ) THEN
    !     PRINT *,'IV =', IV, ', V2H =', MESH%V2H(2,IV), ' (MESH%NT =', MESH%NT, ')'
    !     PAUSE
    !  END IF
    !END DO

  end subroutine split_edge






  subroutine collapse_edge( &
       mesh, &
       brep, &
       hypg, &
       ihedg, &
       stat, &
       resulting_vertex )
    use mod_util
    use mod_types_brep
    use mod_hypergraph
    implicit none
    LOGICAL :: DEBUG = .true.
    type(type_surface_mesh), intent(inout) :: mesh
    type(type_brep),         intent(in)    :: brep
    type(type_hypergraph),   intent(in)    :: hypg
    integer,                 intent(in)    :: ihedg(2)
    integer,                 intent(out)   :: stat
    integer, optional,       intent(out)   :: resulting_vertex
    integer                                :: halfedges(2,2)
    integer, dimension(2)                  :: verts, faces, opposite_vert
    integer                                :: twin_jk(2,2,2)
    integer                                :: i_new_vert, i_removed_vert
    integer                                :: vert_id_new(mesh%nv)
    integer                                :: face_id_new(mesh%nt)
    integer                                :: ih, iv, jv, kv, it, i
    real(kind=fp)                          :: uv_new(2,2)
    real(kind=fp)                          :: xyz_new(3)
    integer                                :: typ_new
    integer                                :: ids_new
    integer                                :: stat_proj
    integer                                :: ihype, ipath
    integer, allocatable                   :: cavity(:)
    integer                                :: ncavity

    verts(1) = get_orig(mesh, ihedg)
    verts(2) = get_dest(mesh, ihedg)

    IF ( DEBUG ) THEN
      PRINT *, ''
      PRINT *, '--- COLLAPSE EDGE ---'
      PRINT *, '   BEFORE COLLAPSE, MESH%NV, MESH%NT =', MESH%NV, MESH%NT
      PRINT *, '   EDGE #', IHEDG(1), ' OF FACE #', IHEDG(2)
      PRINT *, '   VERTS #', VERTS
      PRINT *, '   XYZ_MID =', 0.5_FP*(MESH%XYZ(:,VERTS(1)) + MESH%XYZ(:,VERTS(2)))
    END IF

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
    IF ( DEBUG ) THEN
      PRINT *,'FACES =',FACES
      DO IT = 1,2
         PRINT *, FACES(IT), '  VERTS:', MESH%TRI(1:3,FACES(IT))
      END DO
    END IF

    ! new vertex resulting from collapse
    i_new_vert = minval(verts)
    i_removed_vert = maxval(verts)
    ! uv, typ, ids, xyz***!****
    if ( is_boundary_vertex(mesh, verts(2)) .and. &
         .not.is_boundary_vertex(mesh, verts(1)) ) then
       mesh%xyz(1:3,i_new_vert) = mesh%xyz(1:3,verts(2))
    else
       !mesh%xyz(1:3,i_new_vert) = 0.5_fp*(mesh%xyz(1:3,verts(1)) + mesh%xyz(1:3,verts(2))) !***
      ! 1) sort verts in ascending typ order
      if ( mesh%typ(verts(1)) > mesh%typ(verts(2)) ) verts = verts([2,1])
      typ_new = mesh%typ(verts(1))
      !
      if ( mesh%typ(verts(1)) < mesh%typ(verts(2)) ) then
         ids_new = mesh%ids(verts(1))
         uv_new = mesh%uv(1:2,1:2,verts(1))
         xyz_new = mesh%xyz(1:3,verts(1))
      else ! mesh%typ(verts(1)) == mesh%typ(verts(2))
         select case ( mesh%typ(verts(1)) )
         case (0)
            stat = 2
            return
         case (1)
            call get_midpoint_hyperedge( &
               mesh, &
               verts, &
               brep, &
               hypg, &
               stat_proj, &
               ids_new, &
               uv_new, &
               xyz_new )
            if ( stat_proj > 0 ) then
               print *,'collapse_edge: failed to project on hyperedge'
               stat = stat_proj
               PAUSE
               return
            end if
            !
         case (2)
            IF ( .TRUE. ) THEN
               call get_cavity_verts_around_edge( &
                  mesh, &
                  verts, &
                  cavity, &
                  ncavity )

               call get_centroid_hyperface( &
                  mesh, &
                  cavity(1:ncavity), &
                  ncavity, &
                  brep, &
                  stat_proj, &
                  ids_new, &
                  uv_new(1:2,1), &
                  xyz_new )
            ELSE
               call get_midpoint_hyperface( &
                  mesh, &
                  verts, &
                  brep, &
                  stat_proj, &
                  ids_new, &
                  uv_new(1:2,1), &
                  xyz_new )
            END IF
            if ( stat_proj > 0 ) then
               print *,'collapse_edge: failed to project on hyperface'
               stat = stat_proj
               PAUSE
               return
            end if
            !
         end select
      end if
    end if ! *******

    mesh%typ(i_new_vert) = typ_new
    mesh%ids(i_new_vert) = ids_new
    mesh%uv(1:2,1:2,i_new_vert) = uv_new
    mesh%xyz(1:3,i_new_vert) = xyz_new
    if ( allocated(mesh%hTargetV) ) mesh%hTargetV(i_new_vert) = 0.5_fp*sum(mesh%hTargetV(verts))

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
       IF ( DEBUG ) PRINT *,'OPP. VERT:',IV,', V2H:',mesh%v2h(1:2,iv), '  |', faces(ih)
       if ( mesh%v2h(2,iv) == faces(ih) ) then
          !IF ( DEBUG ) PRINT *,'OPP. VERT:',IV,', V2H:',mesh%v2h(1:2,iv),' -->',twin_jk(1:2,1,ih)
          IF ( DEBUG ) PRINT *,'    <--',twin_jk(1:2,1,ih)
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

    !DO IT = 1,MESH%NT
    !  !IF ( IT == FACES(1) .OR. FACE(2))
    !  IF ( face_id_new(IT) < 0 .OR. face_id_new(IT) > MESH%NT ) THEN
    !     PRINT *,'RENUM FACES :', IT, ' -->', face_id_new(IT), '/', MESH%NT
    !  END IF
    !END DO

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
    mesh%typ(i_removed_vert:mesh%nv-1) = mesh%typ(i_removed_vert+1:mesh%nv)
    mesh%ids(i_removed_vert:mesh%nv-1) = mesh%ids(i_removed_vert+1:mesh%nv)
    mesh%uv(1:2,1:2,i_removed_vert:mesh%nv-1) = mesh%uv(1:2,1:2,i_removed_vert+1:mesh%nv)
    mesh%xyz(1:3,i_removed_vert:mesh%nv-1) = mesh%xyz(1:3,i_removed_vert+1:mesh%nv)
    mesh%v2h(1:2,i_removed_vert:mesh%nv-1) = mesh%v2h(1:2,i_removed_vert+1:mesh%nv)
    if ( allocated(mesh%hTargetV) ) then
      mesh%hTargetV(i_removed_vert:mesh%nv-1) = mesh%hTargetV(i_removed_vert+1:mesh%nv)
    end if

    mesh%nv = mesh%nv - 1

    ! apply re-ordering
    do iv = 1,mesh%nv
       IF ( mesh%v2h(2,iv) > mesh%nt+2 ) THEN
         PRINT *, 'collapse_edge: v2h(2,', iv, '/', mesh%nv, ') =', &
         mesh%v2h(2,iv), '> mesh%nt (=', mesh%nt, ')'
         !PAUSE
       END IF
       mesh%v2h(2,iv) = face_id_new(mesh%v2h(2,iv))
    end do
    ! possibly remove from path
    if ( mesh%typ(i_removed_vert) == 1 ) then
      ihype = brep%edges(mesh%ids(i_removed_vert))%hyperedge
      do ipath = 1,mesh%npaths
         if ( mesh%paths(ipath)%hyperedge == ihype ) then
            call remove_from_list( &
                 i_removed_vert, &
                 mesh%paths(ipath)%verts, &
                 mesh%paths(ipath)%nv )
            exit
         end if
      end do
    end if

    if ( present(resulting_vertex) ) resulting_vertex = i_new_vert

  end subroutine collapse_edge







  subroutine surface_mesh_optimization( &
       mesh, &
       brep, &
       hypg, &
       hmin, &
       hmax, &
       tolchord )
    use mod_types_brep
    use mod_hypergraph
    use mod_optimmesh
    implicit none
    integer, parameter :: min_valence = 5, max_valence = 7, npass = 10
    real(kind=fp), parameter :: tolH = 0.25_fp
    real(kind=fp), parameter :: tolHMax = (1._fp + tolH)**2!2._fp!
    real(kind=fp), parameter :: tolHMin = (1._fp - tolH)**2!0.5_fp!
    LOGICAL :: ALLOW_SPLIT     = .true.
    LOGICAL :: ALLOW_COLLAPSE  = .true.
    LOGICAL :: ALLOW_FLIP      = .true.
    LOGICAL :: ALLOW_SMOOTHING = .true.
    LOGICAL :: DEBUG = .TRUE.
    type(type_surface_mesh), intent(inout) :: mesh
    type(type_brep),         intent(in)    :: brep
    type(type_hypergraph),   intent(in)    :: hypg
    real(kind=fp),           intent(in)    :: hmin
    real(kind=fp),           intent(in)    :: hmax
    real(kind=fp),           intent(in)    :: tolchord
    real(kind=fp)                          :: htarget(2*mesh%nv), h
    logical                                :: visited(3,2*mesh%nt)
    integer                                :: valence(2*mesh%nv)
    integer                                :: verts(4)
    integer                                :: valence_minmax(2,2), valence_tmp(4)
    integer                                :: stat_split, split_count
    integer                                :: stat_collapse, collapse_count
    integer                                :: stat_flip, flip_count
    integer                                :: it, jt, ie, iv, ih(2), ipass, nt0
    CHARACTER(3) :: strnum   

    PRINT *, 'PRESCRIBED H MIN, MAX =', HMIN, HMAX
    nt0 = mesh%nt
    do ipass = 1,npass
      ! 1) split long edges/collapse short edges
      split_count = 0
      collapse_count = 0
      if ( ALLOW_SPLIT .or. ALLOW_COLLAPSE ) then
         call set_hTargetV( &
            mesh, &
            hmin, &
            hmax, &
            tolchord, &
            opt_adaptScale=.false., &
            opt_brep=brep, &
            !opt_nPassesGradation=20, &
            !opt_gradationMax=1.5_fp, &
            opt_fractionPrevious=0._fp )
         htarget(1:mesh%nv) = mesh%hTargetV(1:mesh%nv)
         PRINT *,'HTARGET MIN,MAX =', MINVAL(htarget(1:mesh%nv)), MAXVAL(htarget(1:mesh%nv))
         !PAUSE
      end if

      IF ( DEBUG ) THEN
         write (strnum, '(i3.3)') ipass-1
         if ( allocated(mesh%hTargetV) ) then
            call write_vtk_mesh_sol( &
               mesh, &
               '../debug/surface_mesh_optimization/pass_'//strnum//'.vtk', &
               solv=mesh%hTargetV(1:mesh%nv), &
               solv_label='hTarget' )
         else
            call write_vtk_mesh( &
               mesh, &
               '../debug/surface_mesh_optimization/pass_'//strnum//'.vtk' )
         end if
      END IF

      edge_splits_collapses : do
         !call target_edge_lengths( &
         !   mesh, &
         !   hmin, &
         !   hmax, &
         !   htarget(1:mesh%nv), &
         !   .false., &
         !   brep )
         IF ( split_count > 3*nt0 ) EXIT edge_splits_collapses
         IF ( mesh%nt < 5 ) EXIT edge_splits_collapses
         visited(1:3,1:mesh%nt) = .false.
         do it = 1,mesh%nt
            do ie = 1,3
               if ( visited(ie,it) ) cycle
               visited(ie,it) = .true.
               ih = get_twin(mesh, [ie,it])
               if ( ih(2) < 1 ) cycle
               visited(ih(1),ih(2)) = .true.

               verts(1) = mesh%tri(ie,it)
               verts(2) = mesh%tri(1+mod(ie,3),it)
               h = sum((mesh%xyz(1:3,verts(1)) - mesh%xyz(1:3,verts(2)))**2)
               !if ( h > 4._fp*maxval(htarget(verts(1:2)))**2 .and. ALLOW_SPLIT ) then
               if ( ALLOW_SPLIT .and. &
                  h > tolHMax*maxval(htarget(verts(1:2)))**2  .and. &
                  h > 4._fp*minval(htarget(verts(1:2)))**2 ) then
                  !PRINT *, 'H = ', SQRT(H), ' | ', SQRT(tolHMax)*htarget(verts(1:2))
                  call split_edge( &
                     mesh, &
                     brep, &
                     hypg, &
                     [ie,it] )
                  stat_split = 0
                  if ( stat_split == 0 ) then
                     !PRINT *, '*** HNEW =', MESH%hTargetV(MESH%NV)
                     split_count = split_count + 1
                     if ( mesh%nt > size(visited,2) ) stop
                     if ( mesh%nv > size(htarget) ) stop
                     cycle edge_splits_collapses
                  end if
               !elseif ( h < 0.5_fp*minval(htarget(verts(1:2)))**2 .and. ALLOW_COLLAPSE ) then
               elseif( ALLOW_COLLAPSE .and. h < tolHMin*minval(htarget(verts(1:2)))**2 ) then
                  !PRINT *, 'H = ', SQRT(H), ' | ', SQRT(tolHMin)*htarget(verts(1:2))
                  call collapse_edge( &
                     mesh, &
                     brep, &
                     hypg, &
                     [ie,it], &
                     stat_collapse )
                     if ( stat_collapse == 0 ) then
                        !PRINT *, '*** HNEW =', MESH%hTargetV(MINVAL(VERTS(1:2)))
                        collapse_count = collapse_count + 1
                        cycle edge_splits_collapses
                     end if
               end if
            end do
         end do

         ! we get there if and only if no more edge splits/collapses are performed
         exit edge_splits_collapses
      end do edge_splits_collapses

      PRINT *, split_count, ' EDGE SPLIT(S) (mesh originally had', nt0, ' triangles)'
      PRINT *, collapse_count, ' EDGE COLLAPSE(S) (mesh originally had', nt0, ' triangles)'
      !PAUSE

       ! 2) compute valence
       do iv = 1,mesh%nv
          valence(iv) = 1
          ih = mesh%v2h(1:2,iv)
          it = get_face(ih)
          do
             ih = get_twin(mesh, get_prev(ih))
             jt = get_face(ih)
             if ( jt < 1 .or. jt == it ) exit
             valence(iv) = valence(iv) + 1
          end do
       end do

       ! 3) valence improving edge flipping
       if ( ALLOW_FLIP ) then
         flip_count = 0
         valence_improving_flips : do 
            IF ( FLIP_COUNT > 3*MESH%NT ) EXIT  valence_improving_flips
            visited(1:3,1:mesh%nt) = .false.
            do it = 1,mesh%nt
               do ie = 1,3
                  if ( visited(ie,it) ) cycle
                  visited(ie,it) = .true.
                  ih = get_twin(mesh, [ie,it])
                  if ( ih(2) < 1 ) cycle
                  visited(ih(1),ih(2)) = .true.
                  ! edge vertices
                  verts(1) = mesh%tri(ie,it)
                  verts(2) = mesh%tri(1+mod(ie,3),it)
                  if ( maxval(mesh%typ(verts(1:2))) < 2 ) cycle
                  ! opposite vertices
                  verts(3) = mesh%tri(1+mod(ie+1,3),it)
                  verts(4) = mesh%tri(1+mod(ih(1)+1,3),ih(2))
                  ! compute new valences if the edge gets flipped
                  valence_tmp = valence(verts)
                  valence_tmp(1:2) = valence_tmp(1:2) - 1
                  valence_tmp(3:4) = valence_tmp(3:4) + 1
                  ! compute min/max valence in pair of triangles
                  valence_minmax(1,1:2) = 1000
                  valence_minmax(2,1:2) = 0
                  do iv = 1,4
                     valence_minmax(1,1) = min(valence_minmax(1,1), valence(verts(iv)))
                     valence_minmax(2,1) = max(valence_minmax(2,1), valence(verts(iv)))
                     valence_minmax(1,2) = min(valence_minmax(1,2), valence_tmp(iv))
                     valence_minmax(2,2) = max(valence_minmax(2,2), valence_tmp(iv))
                  end do
                  !
                  if ( valence_minmax(1,1) >= min_valence .and. &
                        valence_minmax(2,1) <= max_valence ) cycle
                  ! check if edge flip would improve valence uniformity
                  if ( valence_minmax(2,2) - valence_minmax(1,2) < &
                        valence_minmax(2,1) - valence_minmax(1,1)) then
                     PRINT *,'BEFORE: MIN/MAX VAL =', valence_minmax(1:2,1)
                     PRINT *,'AFTER:  MIN/MAX VAL =', valence_minmax(1:2,2)
                     ! flip egde
                     call flip_edge( &
                           mesh, &
                           [ie,it], &
                           stat_flip )
                     PRINT *,'STAT_FLIP =', STAT_FLIP
                     PRINT *,''
                     !PAUSE
                     if ( stat_flip == 0 ) then
                        flip_count = flip_count + 1
                        do iv = 1,4
                           valence(verts(iv)) = valence_tmp(iv)
                        end do
                        cycle valence_improving_flips
                     end if
                  end if
               end do
            end do

            ! we get there if and only if no more edge flips are performed
            exit valence_improving_flips ! 

         end do valence_improving_flips

         PRINT *, flip_count, ' EDGE FLIP(S) (mesh contains', mesh%nt, ' triangles)'
      end if

       ! 4) variational mesh optimization
      if ( ALLOW_SMOOTHING ) then
       call optim_jiao( &
            brep, &
            hypg%hyperedges(1:hypg%nhe), &
            hypg%nhe, &
            mesh, &
            1._fp, &
            0.7_fp, &
            2, &
            4, &
            6, &
            hmin, &
            hmax )
      end if
       
      if ( ALLOW_SPLIT .and. split_count < 1 .and. &
           ALLOW_COLLAPSE .and. collapse_count < 1 .and. & 
           ALLOW_FLIP .and. flip_count < 1 ) then
         exit
      end if
       ! 3) energy reduction edge flipping

    end do

    IF ( DEBUG ) THEN
      write (strnum, '(i3.3)') ipass-1
      call set_hTargetV( &
         mesh, &
         hmin, &
         hmax, &
         tolchord, &
         opt_adaptScale=.false., &
         opt_brep=brep, &
         opt_fractionPrevious=0._fp )
      if ( allocated(mesh%hTargetV) ) then
         call write_vtk_mesh_sol( &
            mesh, &
            '../debug/surface_mesh_optimization/pass_'//strnum//'.vtk', &
            solv=mesh%hTargetV(1:mesh%nv), &
            solv_label='hTarget' )
      else
         call write_vtk_mesh( &
            mesh, &
            '../debug/surface_mesh_optimization/pass_'//strnum//'.vtk' )
      end if
   END IF

  end subroutine surface_mesh_optimization











  !subroutine compute_mesh_gradation( &
  !     mesh, &
  !     htarget, &
  !     gradation, &
  !     nedg )
  !  ! cf. "Geometric surface mesh optimization", Frey and Borouchaki (1998)
  !  implicit none
  !  type(type_surface_mesh), intent(in)  :: mesh
  !  real(kind=fp),           intent(in)  :: htarget(mesh%nv)
  !  real(kind=fp),           intent(out) :: gradation(:)
  !  integer, optional,       intent(out) :: nedg
  !  real(kind=fp)                        :: norm_len
  !  integer                              :: it, ie, iv, jv
  !  k = 0
  !  do it = 1,mesh%nt
  !     do ie = 1,3
  !        iv = mesh%tri(ie,it)
  !        jv = mesh%tri(1+mod(ie,3),it)
  !        if ( iv < jv ) then
  !           k = k + 1
  !           norm_len = norm2(mesh%xyz(1:3,iv) - mesh%xyz(1:3,jv))
  !           if ( abs(htarget(iv) - htarget(jv)) < EPSfp ) then
  !              norm_len = norm_len/htarget(iv)
  !           else
  !              norm_len = norm_len*(htarget(iv) - htarget(jv)) / &
  !                   (htarget(iv)*htarget(jv)*(log(htarget(iv)) - log(htarget(jv)))
  !           end if
  !           gradation(k) = max(htarget(iv)/htarget(jv), htarget(jv)/htarget(iv))**(1._fp/norm_len)
  !        end if            
  !     end do
  !  end do
  !  if ( present(nedg) ) nedg = k  
  !end subroutine compute_mesh_gradation

subroutine get_midpoint_hyperedge( &
   mesh, &
   verts, &
   brep, &
   hypg, &
   stat, &
   ids_mid, &
   uv_mid, &
   xyz_mid )
   use mod_types_brep
   use mod_hypergraph
   use mod_diffgeom
   use mod_intersection
   use mod_projection
   implicit none
   LOGICAL :: DEBUG_PROJ = .false.
   type(type_surface_mesh), intent(in)  :: mesh
   integer,                 intent(in)  :: verts(2)
   type(type_brep),         intent(in)  :: brep
   type(type_hypergraph),   intent(in)  :: hypg
   integer,                 intent(out) :: stat
   integer,                 intent(out) :: ids_mid
   real(kind=fp),           intent(out) :: uv_mid(2,2)
   real(kind=fp),           intent(out) :: xyz_mid(3)
   integer                              :: iedge, ihype
   real(kind=fp)                        :: duv_ds(2,2,2), dxyz_ds(3,2)
   real(kind=fp)                        :: ds, duv(2,2), dxyz(3)

   dxyz = 0.5_fp*(mesh%xyz(1:3,verts(2)) - mesh%xyz(1:3,verts(1)))

   iedge = mesh%ids(verts(1))
   ihype = brep%edges(iedge)%hyperedge
   call diffgeom_intersection( &
      brep%edges(iedge)%curve%surf, &
      mesh%uv(1:2,1:2,verts(1)), &
      duv_ds, &
      dxyz_ds, &
      stat )
   if ( stat > 0 ) then
      print *,'get_midpoint_hyperedge: not a tangential intersection curve'
      return
   end if
   ds = dot_product(dxyz, dxyz_ds(1:3,1))
   duv(1:2,1:2) = ds * duv_ds(1:2,1,1:2)
   dxyz = ds * dxyz_ds(1:3,1)
   !
   call projection_hyperedge( &
      brep, &
      hypg%hyperedges(ihype), &
      iedge, &
      mesh%uv(1:2,1:2,verts(1)), &
      mesh%xyz(1:3,verts(1)), &
      duv, &
      dxyz, &
      ids_mid, &
      uv_mid, &
      xyz_mid, &
      DEBUG_PROJ, &
      stat )
   if ( stat > 0 ) then
      print *,'get_midpoint_hyperedge: failed to project on hyperedge'
      return
   end if

end subroutine get_midpoint_hyperedge





subroutine get_midpoint_hyperface( &
   mesh, &
   verts, &
   brep, &
   stat, &
   ids_mid, &
   uv_mid, &
   xyz_mid )
   use mod_types_brep
   use mod_diffgeom
   use mod_intersection
   use mod_projection
   implicit none
   LOGICAL :: DEBUG_PROJ = .false.
   type(type_surface_mesh), intent(in)  :: mesh
   integer,                 intent(in)  :: verts(2)
   type(type_brep),         intent(in)  :: brep
   integer,                 intent(out) :: stat
   integer,                 intent(out) :: ids_mid
   real(kind=fp),           intent(out) :: uv_mid(2)
   real(kind=fp),           intent(out) :: xyz_mid(3)
   integer                              :: iface, ivar
   real(kind=fp)                        :: dxyz_duv(3,2)
   real(kind=fp)                        :: duv(2), dxyz(3)
   logical                              :: singular

   dxyz = 0.5_fp*(mesh%xyz(1:3,verts(2)) - mesh%xyz(1:3,verts(1)))

   iface = mesh%ids(verts(1))
   do ivar = 1,2 ! <-------------------+
      call evald1( &                   !
         dxyz_duv(1:3,ivar), &         !
         brep%faces(iface)%surface, &  !
         mesh%uv(1:2,1,verts(1)), &    !
         ivar )                        !
   end do ! <--------------------------+
   call solve_NxN( &
      duv(1:2), &
      matmul(transpose(dxyz_duv), dxyz_duv), &
      matmul(transpose(dxyz_duv), dxyz), &
      singular )
   dxyz = matmul(dxyz_duv, duv(1:2))
   !
   call projection_hyperface( &
      brep, &
      iface, &
      mesh%uv(1:2,1,verts(1)), &
      mesh%xyz(1:3,verts(1)), &
      duv, &
      dxyz, &
      ids_mid, &
      uv_mid, &
      DEBUG_PROJ, &
      stat )
   if ( stat > 0 ) then
      print *,'get_midpoint_hyperface: failed to project on hyperedge'
      return
   else
      call eval( &
         xyz_mid, &
         brep%faces(ids_mid)%surface, &
         uv_mid ) 
   end if

end subroutine get_midpoint_hyperface














subroutine get_centroid_hyperface( &
   mesh, &
   verts, &
   nverts, &
   brep, &
   stat, &
   ids_mid, &
   uv_mid, &
   xyz_mid )
   use mod_types_brep
   use mod_diffgeom
   use mod_intersection
   use mod_projection
   implicit none
   LOGICAL :: DEBUG_PROJ = .false.
   type(type_surface_mesh), intent(in)  :: mesh
   integer,                 intent(in)  :: nverts
   integer,                 intent(in)  :: verts(nverts)
   type(type_brep),         intent(in)  :: brep
   integer,                 intent(out) :: stat
   integer,                 intent(out) :: ids_mid
   real(kind=fp),           intent(out) :: uv_mid(2)
   real(kind=fp),           intent(out) :: xyz_mid(3)
   integer                              :: ivert, iface, ivar
   real(kind=fp)                        :: dxyz_duv(3,2)
   real(kind=fp)                        :: duv(2), dxyz(3)
   logical                              :: singular

   xyz_mid(1:3) = 0._fp
   do ivert = 1,nverts
      xyz_mid(1:3) = xyz_mid(1:3) + mesh%xyz(1:3,verts(ivert))
   end do
   xyz_mid = xyz_mid/real(nverts, kind=fp)
   !dxyz = 0.5_fp*(mesh%xyz(1:3,verts(2)) - mesh%xyz(1:3,verts(1)))
   dxyz = xyz_mid - mesh%xyz(1:3,verts(1))

   iface = mesh%ids(verts(1))
   do ivar = 1,2 ! <-------------------+
      call evald1( &                   !
         dxyz_duv(1:3,ivar), &         !
         brep%faces(iface)%surface, &  !
         mesh%uv(1:2,1,verts(1)), &    !
         ivar )                        !
   end do ! <--------------------------+
   call solve_NxN( &
      duv(1:2), &
      matmul(transpose(dxyz_duv), dxyz_duv), &
      matmul(transpose(dxyz_duv), dxyz), &
      singular )
   dxyz = matmul(dxyz_duv, duv(1:2))
   !
   call projection_hyperface( &
      brep, &
      iface, &
      mesh%uv(1:2,1,verts(1)), &
      mesh%xyz(1:3,verts(1)), &
      duv, &
      dxyz, &
      ids_mid, &
      uv_mid, &
      DEBUG_PROJ, &
      stat )
   if ( stat > 0 ) then
      print *,'get_midpoint_hyperface: failed to project on hyperedge'
      return
   else
      call eval( &
         xyz_mid, &
         brep%faces(ids_mid)%surface, &
         uv_mid ) 
   end if

end subroutine get_centroid_hyperface



subroutine get_cavity_verts_around_edge( &
   mesh, &
   edge_verts, &
   cavity, &
   ncavity )
   use mod_util
   implicit none
   type(type_surface_mesh), intent(in)    :: mesh
   integer,                 intent(in)    :: edge_verts(2)
   integer, allocatable,    intent(inout) :: cavity(:)
   integer,                 intent(out)   :: ncavity
   logical                                :: in_cavity
   integer                                :: iv, jv, kv, it, ih(2)
   INTEGER :: CNT

   if ( allocated(cavity) ) then
      if ( size(cavity) < 2 ) deallocate(cavity)
   end if
   allocate(cavity(6))

   cavity(1:2) = edge_verts(1:2)
   ncavity = 2

   !PRINT *,'GET CAVITY AROUND EDGE (VERTS=', edge_verts, ')'
   !PRINT *,'CAVITY =', CAVITY(1:NCAVITY)

   do iv = 1,2
      ih = mesh%v2h(1:2,edge_verts(iv))
      it = ih(2)
      CNT = 0
      adjacent_verts: do 
         CNT = CNT + 1
         IF ( CNT > PARAM_max_cycles_around_mesh_vertex ) THEN
            PRINT *,'get_cavity_verts_around_edge: wrong cycling from vert#', edge_verts(iv)
            PRINT *,'CAVITY =', CAVITY(1:NCAVITY)
            PRINT *,'XYZ =', MESH%XYZ(1:3,edge_verts(iv))
            call write_vtk_mesh( &
               mesh, &
               '../debug/get_cavity_verts_around_edge.vtk' )
            PAUSE
         END IF
         jv = get_dest(mesh, ih)
         in_cavity = .false.
         check_in_cavity: do kv = 1,ncavity
            if ( cavity(kv) == jv ) then
               in_cavity = .true.
               exit check_in_cavity
            end if
         end do check_in_cavity

         if ( .not. in_cavity ) then
            call insert_after( &
               cavity, &
               ncavity, &
               jv, &
               ncavity )
            !PRINT *,'CAVITY =', CAVITY(1:NCAVITY)
         end if

         ih = get_twin(mesh, get_prev(ih))
         if ( ih(2) == it .or. ih(2) < 1 ) exit adjacent_verts
      end do adjacent_verts
   end do

end subroutine get_cavity_verts_around_edge

end module mod_mesh_optimization
