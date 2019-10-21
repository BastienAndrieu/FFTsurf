module mod_mesh_optimization

  use mod_math
  use mod_halfedge
  use mod_mesh

contains

  subroutine flip_edge( &
       mesh, &
       ihedg, &
       stat )
    implicit none
    LOGICAL :: DEBUG = .true.
    type(type_surface_mesh), intent(inout) :: mesh
    integer,                 intent(in)    :: ihedg(2)
    integer,                 intent(out)   :: stat
    integer                                :: halfedges(2,2)
    integer                                :: ijk(3,2)
    integer                                :: opposite_vert(2)
    integer                                :: twin_j(2,2)
    integer                                :: ih, jh, iv, it, jt

    IF ( DEBUG ) THEN
       PRINT *, '--- SWAP EDGE ---'
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
       ihedg )
    use mod_util
    use mod_types_brep
    use mod_hypergraph
    use mod_diffgeom
    use mod_intersection
    use mod_projection
    implicit none
    LOGICAL :: DEBUG = .true.
    LOGICAL :: DEBUG_PROJ = .FALSE.
    type(type_surface_mesh), intent(inout) :: mesh
    type(type_brep),         intent(in)    :: brep
    type(type_hypergraph),   intent(in)    :: hypg
    integer,                 intent(in)    :: ihedg(2)
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
    !
    real(kind=fp)                          :: dxyz(3), duv(2,2), ds
    real(kind=fp)                          :: dxyz_ds(3,2), duv_ds(2,2,2)
    real(kind=fp)                          :: dxyz_duv(3,2)
    integer                                :: stat_tangent, stat_proj
    logical                                :: singular
    integer                                :: iedge, ihype, iface, ivar

    verts(1) = get_orig(mesh, ihedg)
    verts(2) = get_dest(mesh, ihedg)

    IF ( DEBUG ) THEN
       PRINT *, '--- SPLIT EDGE ---'
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
    ! 1) sort verts in descending typ order
    if ( mesh%typ(verts(1)) <  mesh%typ(verts(2))) verts = verts([2,1])
    dxyz = 0.5_fp*(mesh%xyz(1:3,verts(2)) - mesh%xyz(1:3,verts(1)))
    typ_new = mesh%typ(verts(1))

    select case ( mesh%typ(verts(1)) ) 
    case (1)
       iedge = mesh%ids(verts(1))
       ihype = brep%edges(iedge)%hyperedge
       call diffgeom_intersection( &
            brep%edges(iedge)%curve%surf, &
            mesh%uv(1:2,1:2,verts(1)), &
            duv_ds, &
            dxyz_ds, &
            stat_tangent )
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
            ids_new, &
            uv_new, &
            xyz_new, &
            DEBUG_PROJ, &
            stat_proj )
       if ( stat_proj > 0 ) then
          print *,'split_edge: failed to project on hyperedge'
          PAUSE
       end if
       !
    case (2)
       iface = mesh%ids(verts(1))
       do ivar = 1,2 ! <---------------------+
          call evald1( &                     !
               dxyz_duv(1:3,ivar), &         !
               brep%faces(iface)%surface, &  !
               mesh%uv(1:2,1,verts(1)), &    !
               ivar )                        !
       end do ! <----------------------------+
       call solve_NxN( &
            duv(1:2,1), &
            matmul(transpose(dxyz_duv), dxyz_duv), &
            matmul(transpose(dxyz_duv), dxyz), &
            singular )
       dxyz = matmul(dxyz_duv, duv(1:2,1))
       !
       call projection_hyperface( &
            brep, &
            iface, &
            mesh%uv(1:2,1,verts(1)), &
            mesh%xyz(1:3,verts(1)), &
            duv, &
            dxyz, &
            ids_new, &
            uv_new(1:2,1), &
            DEBUG_PROJ, &
            stat_proj )
       if ( stat_proj > 0 ) then
          print *,'split_edge: failed to project on hyperedge'
          PAUSE
       else
          call eval( &
               xyz_new, &
               brep%faces(ids_new)%surface, &
               uv_new(1:2,1) ) 
       end if
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







  subroutine surface_mesh_optimization( &
       mesh, &
       brep, &
       hypg, &
       hmin, &
       hmax &
       )
    use mod_types_brep
    use mod_hypergraph
    use mod_optimmesh
    implicit none
    integer, parameter :: min_valence = 5, max_valence = 7, npass = 1
    type(type_surface_mesh), intent(inout) :: mesh
    type(type_brep),         intent(in)    :: brep
    type(type_hypergraph),   intent(in)    :: hypg
    real(kind=fp),           intent(in)    :: hmin
    real(kind=fp),           intent(in)    :: hmax
    logical                                :: visited(3,mesh%nt)
    integer                                :: valence(mesh%nv)
    integer                                :: verts(4)
    integer                                :: valence_minmax(2,2), valence_tmp(4)
    integer                                :: stat_flip, flip_count
    integer                                :: it, jt, ie, iv, ih(2), ipass

    

    do ipass = 1,npass
       ! 0) compute valence
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

       ! 1) valence improving edge flipping
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

       ! 2) variational mesh optimization
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
       
       ! 3) energy reduction edge flipping
    end do
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





end module mod_mesh_optimization
