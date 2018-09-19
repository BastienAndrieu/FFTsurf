module mod_optimmesh

  use mod_math

  implicit none

contains


  subroutine contract_small_edges( &
       brep, &
       hyperedges, &
       nhe, &
       mesh, &
       hmin )
    use mod_types_brep
    use mod_mesh
    use mod_hypergraph
    implicit none
    type(type_brep),         intent(in)    :: brep
    integer,                 intent(in)    :: nhe
    type(type_hyperedge),    intent(in)    :: hyperedges(nhe)
    type(type_surface_mesh), intent(inout) :: mesh
    real(kind=fp),           intent(in)    :: hmin
    real(kind=fp)                          :: hminsqr, h
    integer                                :: typij(2), pairv(2)
    integer                                :: iface, iedge, ivert, jvert, ihype, jhype
    integer :: k

    hminsqr = hmin**2

    k = 0
    outer_loop : do
       k = k + 1
       PRINT *,''
       PRINT *,'OUTER, K=',k
       !IF ( K > 3 ) EXIT
       do iface = 1,mesh%nt
          halfedges : do iedge = 1,3
             ivert = mesh%tri(iedge,iface)
             jvert = mesh%tri(1+mod(iedge,3),iface)
             h = sum((mesh%xyz(:,jvert) - mesh%xyz(:,ivert))**2)
             if ( h < hminsqr ) then
                ! check if this edge can be contracted
                typij = mesh%typ([ivert,jvert])
                if ( typij(1) < typij(2) ) then
                   pairv = [ivert, jvert]
                else
                   pairv = [jvert,ivert]
                   typij = typij([2,1])
                end if

                if ( typij(2) < 1 ) cycle halfedges ! this edge cannot be contracted
                if ( typij(1) == 0 ) then
                   if ( typij(2) == 1 ) then
                      jhype = brep%edges(mesh%ids(pairv(2)))%hyperedge
                      if ( all(hyperedges(jhype)%verts /= mesh%ids(pairv(1))) ) cycle halfedges
                   end if
                elseif ( typij(1) == 1 ) then
                   if ( typij(2) == 1 ) then
                      ihype = brep%edges(mesh%ids(pairv(1)))%hyperedge
                      jhype = brep%edges(mesh%ids(pairv(2)))%hyperedge
                      if ( ihype /= jhype ) cycle halfedges
                   end if
                end if
                
                PRINT *,'CONTRACT EDGE',PAIRV
                call edge_contraction( &
                     brep, &
                     hyperedges, &
                     nhe, &
                     mesh, &
                     [iedge,iface] )
                !RETURN
                cycle outer_loop
             end if
          end do halfedges
       end do

       exit
    end do outer_loop
    PRINT *,'DONE'

  end subroutine contract_small_edges





  subroutine edge_contraction( &
       brep, &
       hyperedges, &
       nhe, &
       mesh, &
       ihedg )
    use mod_util
    use mod_linalg
    use mod_diffgeom
    use mod_types_brep
    use mod_mesh
    use mod_hypergraph
    use mod_projection
    use mod_intersection
    use mod_halfedge
    implicit none
    type(type_brep),         intent(in)    :: brep
    integer,                 intent(in)    :: nhe
    type(type_hyperedge),    intent(in)    :: hyperedges(nhe)
    type(type_surface_mesh), intent(inout) :: mesh
    integer,                 intent(in)    :: ihedg(2)
    real(kind=fp)                          :: xyznew(3), uvnew(2,2)
    real(kind=fp)                          :: tng(3,2), ds, duv_ds(2,2,2), duv(2,2), dxyz(3)
    logical                                :: singular
    integer                                :: idsnew, typnew, stat
    integer                                :: verts(2), faces(2), i, j, k, ivar, inew
    integer                                :: it, iv, ip, ihype
    integer                                :: v_old2new(mesh%nv), t_old2new(mesh%nt)

    verts = mesh%tri([ihedg(1), 1+mod(ihedg(1),3)],ihedg(2))

    ! geometry
    if ( mesh%typ(verts(1)) == mesh%typ(verts(2)) ) then
       ! midpoint
       typnew = mesh%typ(verts(1))
       if ( mesh%ids(verts(1)) == mesh%ids(verts(2)) ) then
          uvnew  = 0.5 * sum(mesh%uv(:,:,verts), 3)
          idsnew = mesh%ids(verts(1))
          if ( typnew == 1 ) then
             call eval( &
                  xyznew, &
                  brep%edges(idsnew)%curve%surf(2)%ptr, &
                  uvnew(:,1) )
          else
             call eval( &
                  xyznew, &
                  brep%faces(idsnew)%surface, &
                  uvnew(:,1) )
          end if
       else
          dxyz = 0.5_fp * (mesh%xyz(:,verts(2)) - mesh%xyz(:,verts(1)))
          if ( typnew == 1 ) then
             call diffgeom_intersection( &
                  brep%edges(idsnew)%curve%surf, &
                  mesh%uv(:,:,verts(1)), &
                  duv_ds, &
                  tng, &
                  stat )
             ds = dot_product(dxyz, tng(:,1))
             duv = ds * duv_ds(:,1,:)
             dxyz = ds * tng(:,1)
             call projection_hyperedge( &
                  brep, &
                  hyperedges(brep%edges(mesh%ids(verts(1)))%hyperedge), &
                  mesh%ids(verts(1)), &
                  mesh%uv(:,:,verts(1)), &
                  mesh%xyz(:,verts(1)), &
                  duv, &
                  dxyz, &
                  idsnew, &
                  uvnew, &
                  .true. )
             call eval( &
                  xyznew, &
                  brep%edges(idsnew)%curve%surf(2)%ptr, &
                  uvnew(:,1) )
          else
             do ivar = 1,2
                call evald1( &
                     tng(:,ivar), &
                     brep%faces(mesh%ids(verts(1)))%surface, &
                     mesh%uv(:,1,verts(1)), &
                     ivar )
             end do
             call solve_NxN( &
                  duv(:,1), &
                  matmul(transpose(tng), tng), &
                  matmul(transpose(tng), dxyz), &
                  singular )
             call projection_hyperface( &
                  brep, &
                  mesh%ids(verts(1)), &
                  mesh%uv(:,1,verts(1)), &
                  mesh%xyz(:,verts(1)), &
                  duv(:,1), &
                  dxyz, &
                  idsnew, &
                  uvnew(:,1), &
                  .true. )
          end if
       end if
    else
       !i = get_orig(mesh, ihedg)
       i = verts(minloc(mesh%typ(verts),1))
       ! contract towards i
       xyznew = mesh%xyz(:,i)
       uvnew  = mesh%uv(:,:,i)
       idsnew = mesh%ids(i)
       typnew = mesh%typ(i)
    end if

    PRINT *,'    IDS =',mesh%ids(verts)
    PRINT *,'    TYP =',mesh%typ(verts)

    ! new topology
    inew = minval(verts)
    mesh%ids(inew) = idsnew
    mesh%typ(inew) = typnew
    mesh%xyz(:,inew) = xyznew
    mesh%uv(:,:,inew) = uvnew

    PRINT *,'  I_NEW =',inew
    PRINT *,'XYZ_NEW =',xyznew
    PRINT *,' UV_NEW =',uvnew
    PRINT *,'TYP_NEW =',typnew
    PRINT *,'IDS_NEW =',idsnew

    i = maxval(verts)
    v_old2new(1:i) = [(j, j=1,i)]
    v_old2new(i) = inew
    v_old2new(i+1:mesh%nv) = [(j, j=i,mesh%nv-1)]
    if ( mesh%typ(i) == 1 ) then
       ihype = brep%edges(mesh%ids(i))%hyperedge
       do ip = 1,mesh%npaths
          if ( mesh%paths(ip)%hyperedge == ihype ) then
             call remove_from_list( &
                  i, &
                  mesh%paths(ip)%verts, &
                  mesh%paths(ip)%nv )
             exit
          end if
       end do
    end if
    
    ! remove vertex
    mesh%typ(i:mesh%nv-1) = mesh%typ(i+1:mesh%nv)
    mesh%ids(i:mesh%nv-1) = mesh%ids(i+1:mesh%nv)
    mesh%xyz(1:3,i:mesh%nv-1) = mesh%xyz(1:3,i+1:mesh%nv)
    mesh%uv(1:2,1:2,i:mesh%nv-1) = mesh%uv(1:2,1:2,i+1:mesh%nv)
    mesh%v2h(1:2,i:mesh%nv-1) = mesh%v2h(1:2,i+1:mesh%nv)
    mesh%nv = mesh%nv - 1

    ! remove triangles
    faces(1) = ihedg(2)
    triangles : do it = 1,mesh%nt
       do iv = 1,3
          if ( mesh%tri(1+mod(iv,3),it) == verts(1) .and. mesh%tri(iv,it) == verts(2) ) then
             faces(2) = it
             exit triangles
          end if
       end do
    end do triangles

    i = minval(faces)
    j = maxval(faces)
    t_old2new(1:i-1) = [(k, k=1,i-1)]
    t_old2new(i+1:j-1) = [(k, k=i,j-2)]
    t_old2new(j+1:mesh%nt) = [(k, k=j-1,mesh%nt-2)]

    mesh%tri(1:3,i:j-2) = mesh%tri(1:3,i+1:j-1)
    mesh%tri(1:3,j-1:mesh%nt-2) = mesh%tri(1:3,j+1:mesh%nt)
    mesh%nt = mesh%nt - 2

    do it = 1,mesh%nt
       mesh%tri(:,it) = v_old2new(mesh%tri(:,it))
    end do

    do ip = 1,mesh%npaths
       mesh%paths(ip)%verts(1:mesh%paths(ip)%nv) = v_old2new(mesh%paths(ip)%verts(1:mesh%paths(ip)%nv))
    end do

  end subroutine edge_contraction



  subroutine smooth_function_mesh( &
       mesh, &
       fv, &
       npass )
    use mod_mesh
    use mod_halfedge
    implicit none
    type(type_surface_mesh), intent(in)    :: mesh
    real(kind=fp),           intent(inout) :: fv(mesh%nv)
    integer,                 intent(in)    :: npass
    real(kind=fp)                          :: tmp(mesh%nv), w, sumw
    integer                                :: ipass, ivert, jvert, ihedg(2), iface, jface, n

    do ipass = 1,npass
       tmp(1:mesh%nv) = 0._fp
       do ivert = 1,mesh%nv
          ihedg = mesh%v2h(:,ivert)
          iface = get_face(ihedg)
          sumw = 0._fp
          n = 0
          adjacent_verts : do
             n = n + 1
             jvert = get_dest(mesh, ihedg)
             w = 1._fp / norm2(mesh%xyz(:,jvert) - mesh%xyz(:,ivert))
             tmp(ivert) = tmp(ivert) + w * fv(jvert)
             sumw = sumw + w

             ihedg = get_prev(ihedg)
             ihedg = get_twin(mesh, ihedg)
             jface = get_face(ihedg)
             if ( jface < 1 .or. jface == iface ) exit adjacent_verts
          end do adjacent_verts
          w = sumw / real(n, kind=fp)
          tmp(ivert) = tmp(ivert) + w * fv(ivert)
          sumw = sumw + w
          tmp(ivert) = tmp(ivert) / sumw
       end do
       fv(1:mesh%nv) = tmp(1:mesh%nv)
    end do

  end subroutine smooth_function_mesh


  subroutine discrete_minimum_curvature_radius( &
       mesh, &
       rad )
    use mod_mesh
    use mod_halfedge
    implicit none
    type(type_surface_mesh), intent(in)  :: mesh
    real(kind=fp),           intent(out) :: rad(mesh%nv)
    real(kind=fp)                        :: vec(3), norv(3,mesh%nv), veci(3), veck(3), denom
    integer                              :: it, jt, iv, jv, kv, lv, ih(2), ip

    norv(:,:) = 0._fp
    ! compute triangle pseudo-normals
    do it = 1,mesh%nt
       vec = cross( &
            mesh%xyz(:,mesh%tri(2,it)) - mesh%xyz(:,mesh%tri(1,it)), &
            mesh%xyz(:,mesh%tri(3,it)) - mesh%xyz(:,mesh%tri(1,it)) )
       do jv = 1,3
          iv = mesh%tri(jv,it)
          norv(:,iv) = norv(:,iv) + vec
       end do
    end do

    ! compute discrete minimum curvature radius at vertices
    rad(:) = huge(1._fp)
    do iv = 1,mesh%nv
       norv(:,iv) = norv(:,iv) / norm2(norv(:,iv))

       ih = mesh%v2h(:,iv)
       it = get_face(ih)
       do
          jv = get_dest(mesh, ih)
          vec = mesh%xyz(:,jv) - mesh%xyz(:,iv)
          rad(iv) = min(rad(iv), sum(vec**2)/abs(dot_product(vec, norv(:,iv))))
          ih = get_prev(ih)
          ih = get_twin(mesh, ih)
          jt = get_face(ih)
          if ( jt < 1 .or. jt == it ) exit
       end do
    end do

    do ip = 1,mesh%npaths
       do lv = 1,mesh%paths(ip)%nv
          iv = mesh%paths(ip)%verts(1 + mod(lv + mesh%paths(ip)%nv - 2,mesh%paths(ip)%nv))
          jv = mesh%paths(ip)%verts(lv)
          kv = mesh%paths(ip)%verts(1 + mod(lv,mesh%paths(ip)%nv))

          veci = mesh%xyz(:,iv) - mesh%xyz(:,jv)
          veck = mesh%xyz(:,kv) - mesh%xyz(:,jv)
          
          denom = 1._fp - dot_product(veck, veci)**2 / (sum(veck**2)*sum(veci**2))
          !PRINT *,'DENOM =',DENOM
          if ( denom < EPSfp ) then
             rad(jv) = huge(1._fp)
          else
             rad(jv) = norm2(mesh%xyz(:,kv) - mesh%xyz(:,iv)) / sqrt(denom)
          end if

          !PRINT *,'TRIANGLE VERTICES :'
          !PRINT *,MESH%XYZ(:,IV)
          !PRINT *,MESH%XYZ(:,JV)
          !PRINT *,MESH%XYZ(:,KV)
          !PRINT *,'CIRCUMRADIUS =',rad(mesh%paths(ip)%verts(jv))
       end do
    end do

  end subroutine discrete_minimum_curvature_radius





  subroutine optim_jiao( &
       brep, &
       hyperedges, &
       nhe, &
       mesh, &
       passmax, &
       frac_conf, &
       hmin, &
       hmax )
    USE MOD_UTIL
    use mod_diffgeom
    use mod_intersection
    use mod_hypergraph
    use mod_types_brep
    use mod_brep
    use mod_mesh
    use mod_halfedge
    use mod_linalg
    use mod_projection
    USE MOD_TOLERANCES
    implicit none
    integer, parameter :: ipass1 = 5, ipass2 = 15
    real(kind=fp), parameter :: frac_conf1 = 1._fp, frac_conf2 = 0.75_fp
    type(type_brep),         intent(in)    :: brep
    integer,                 intent(in)    :: nhe
    type(type_hyperedge),    intent(in)    :: hyperedges(nhe)
    type(type_surface_mesh), intent(inout) :: mesh
    integer,                 intent(in)    :: passmax
    real(kind=fp),           intent(in)    :: frac_conf
    real(kind=fp),           intent(in)    :: hmin
    real(kind=fp),           intent(in)    :: hmax
    real(kind=fp)                          :: hve(mesh%nv)
    real(kind=fp)                          :: wei(mesh%nt)
    real(kind=fp)                          :: ener(mesh%nt)
    real(kind=fp)                          :: grad(3,mesh%nv)
    real(kind=fp)                          :: hess(3,3,mesh%nv)
    real(kind=fp)                          :: duv(2,2,mesh%nv), ds
    real(kind=fp)                          :: tng(3,2), duv_ds(2,2,2)
    integer                                :: stat
    logical                                :: singular
    real(kind=fp), dimension(4)            :: lowerb, upperb
    real(kind=fp)                          :: uvtmp(2,2), dxyz(3)
    integer, dimension(mesh%nv)            :: idsnew, typnew
    real(kind=fp)                          :: uvnew(2,2,mesh%nv)
    integer                                :: ivert, jvert, ihedg(2)
    integer                                :: iedge, ihype
    integer                                :: iface, jface
    integer                                :: ipass, ivar
    logical                                :: check_change
    real(kind=fp)                          :: frac_conf_ramp
    CHARACTER(2) :: STRNUM
    INTEGER :: FID

    CALL GET_FREE_UNIT(FID)

    lowerb(:) = -2._fp
    upperb(:) = 2._fp

    do ipass = 1,passmax
       PRINT *,'OPTIM PASS',IPASS,'/',PASSMAX
       ! compute target edge lengths at vertices
       if ( .false. ) then!ipass < 6 ) then
          hve(:) = 1._fp
       else
          call discrete_minimum_curvature_radius( &
               mesh, &
               hve )
          !hve = min(minval(hve)*hmax/hmin, hve)
          PRINT *,'MIN(HVE) =', real(minval(hve)), minloc(hve,1)
          hve = hve/minval(hve)
          hve = min(hmax/hmin, hve)
          !hve = min(10._fp, hve)
          ! normalize to [0,1] for plotting
          hve = (hve - minval(hve)) / (maxval(hve) - minval(hve))

          call smooth_function_mesh( &
               mesh, &
               hve, &
               10 )
       end if

       IF ( .TRUE. ) THEN
          WRITE (STRNUM, '(I2.2)') IPASS
          call write_tecplot_mesh_solv( &
               mesh, &
               'optimmesh/optim_'//strnum//'.dat', &
               'pass_'//strnum, &
               hve )
       END IF

       ! ...
       ! compute triangle weights
       wei(:) = 1._fp
       ! ...
       ! compute energy, gradient and hessian
       if ( ipass < ipass1 ) then
          frac_conf_ramp = frac_conf1
       elseif ( ipass > ipass2 ) then
          frac_conf_ramp = frac_conf2
       else
          frac_conf_ramp = frac_conf1 + ( real(ipass - ipass1, kind=fp) / real(ipass2 - ipass1, kind=fp) ) * &
               (frac_conf2 - frac_conf1) 
       end if
       PRINT *,'FRAC_CONF_RAMP =',frac_conf_ramp

       call energy_jiao( &
            mesh%nt, &
            mesh%tri(1:3,1:mesh%nt), &
            wei, &
            mesh%nv, &
            mesh%xyz(1:3,1:mesh%nv), &
            hve, &
            frac_conf_ramp, &
            ener, &
            grad, &
            hess )
       !
       ! compute vertex displacements
       !OPEN(UNIT=FID, FILE='Jouke/meshgen/brepmesh/debug_dxyz.dat', ACTION='WRITE')
       compute_duv : do ivert = 1,mesh%nv
          !PRINT *,'IVERT =', IVERT
          select case ( mesh%typ(ivert) ) ! <-------------------------------+
          case (0) ! -------------------------------------------------------+
             uvnew(:,:,ivert) = mesh%uv(:,:,ivert)                          !
             idsnew(ivert) = mesh%ids(ivert)
             typnew(ivert) = 0
             DXYZ(:) = 0._FP
          case (1) ! -------------------------------------------------------+
             call diffgeom_intersection( &                                  !
                  brep%edges(mesh%ids(ivert))%curve%surf, &                 !
                  mesh%uv(:,:,ivert), &                                     !
                  duv_ds, &                                                 !
                  tng, &                                                    !
                  stat )                                                    !
             ds = - dot_product(tng(:,1), grad(:,ivert)) / &                !
                  dot_product(tng(:,1), matmul(hess(:,:,ivert), tng(:,1)))  !
             duv(:,:,ivert) = ds * duv_ds(:,1,:)                            !
             dxyz = ds * tng(:,1)
             !!
             typnew(ivert) = 1
             uvnew(:,:,ivert) = mesh%uv(:,:,ivert)!uvtmp
             iedge = mesh%ids(ivert)   ! brep edge index
             ihype = brep%edges(iedge)%hyperedge ! hyperedge index
             !PRINT *,''; PRINT *,'';
             !PRINT *, IVERT!, ihype
             call projection_hyperedge( &
                  brep, &
                  hyperedges(ihype), &
                  iedge, &
                  mesh%uv(:,:,ivert), &
                  mesh%xyz(:,ivert), &
                  duv(:,:,ivert), &
                  dxyz, &
                  idsnew(ivert), &
                  uvtmp, &
                  .false. )!(ivert == 2045) )
             uvnew(:,:,ivert) = uvtmp
             !
          case (2) ! -------------------------------------------------------
             do ivar = 1,2 ! <-------------------------------+
                call evald1( &                               !
                     tng(:,ivar), &                          !
                     brep%faces(mesh%ids(ivert))%surface, &  !
                     mesh%uv(:,1,ivert), &                   !
                     ivar )                                  !
             end do ! <--------------------------------------+
             call solve_NxN( &
                  duv(:,1,ivert), &
                  matmul(matmul(transpose(tng), hess(:,:,ivert)), tng), &
                  -matmul(transpose(tng), grad(:,ivert)), &
                  singular )
             dxyz = matmul(tng, duv(:,1,ivert))
             !
             ! handle passing to an adjacent face (crossing of a smooth edge...)
             idsnew(ivert) = mesh%ids(ivert)
             typnew(ivert) = 2
             uvnew(:,1,ivert) = mesh%uv(:,1,ivert) + duv(:,1,ivert)
             uvnew(:,2,ivert) = 0._fp
             ihedg = mesh%v2h(:,ivert) ! mesh halfedge index
             iface = mesh%ids(ivert)   ! brep face index
             !
             check_change = ( maxval(abs(duv(:,1,ivert))) > &
                  1._fp + EPSuv - maxval(abs(mesh%uv(:,1,ivert))) )
             if ( .not.check_change ) then
                adjacent_verts2 : do ! <-----------------------------------------------------+
                   jvert = get_dest(mesh, ihedg) ! mesh vertex index                         !
                   jface = mesh%ids(jvert)       ! brep face index                           !
                   ! (quasi) necessary conditions for a change of supporting brep face:      !
                   ! - at least one adjacent vertex is supported by a different brep face;   !
                   ! - this vertex is in the halfspace pointed by the xyz displacement.      !
                   if ( mesh%typ(jvert) /= 2 .or. jface /= iface ) then ! <---------------+  !
                      if ( dot_product(dxyz, &                                            !  !
                           mesh%xyz(:,jvert) - mesh%xyz(:,ivert)) > 0._fp ) then ! <---+  !  !
                         check_change = .true.                                         !  !  !
                         exit adjacent_verts2                                          !  !  !
                      end if ! <-------------------------------------------------------+  !  !
                   end if ! <-------------------------------------------------------------+  !
                   ! move on to next adjacent vertex                                         !
                   ihedg = get_prev(ihedg)       ! outgoing mesh halfedge                    !
                   ihedg = get_twin(mesh, ihedg) ! ingoing mesh halfedge                     !
                   if ( ihedg(2) < 1 .or. ihedg(2) ==  mesh%v2h(2,ivert) ) exit              !
                end do adjacent_verts2 ! <---------------------------------------------------+
             end if

             if ( check_change ) then
                !PRINT *,'IVERT =',IVERT
                call projection_hyperface( &
                     brep, &
                     iface, &
                     mesh%uv(:,1,ivert), &
                     mesh%xyz(:,ivert), &
                     duv(:,1,ivert), &
                     dxyz, &
                     idsnew(ivert), &
                     uvtmp(:,1), &
                     .false. )!(ivert == 32406) )!
                uvnew(:,1,ivert) = uvtmp(:,1)
             end if
             !
          end select ! <----------------------------------------------------+
          !
          IF ( MAXVAL(ABS(UVNEW(:,:,IVERT))) > 1._FP + EPSUV ) PRINT *, IVERT, TYPNEW(IVERT), IDSNEW(IVERT), UVNEW(:,:,IVERT)
          !WRITE (FID,*) DXYZ
       end do compute_duv
       !CLOSE(FID)

       ! update vertices coordinates
       mesh%ids(1:mesh%nv) = idsnew(1:mesh%nv)
       mesh%typ(1:mesh%nv) = typnew(1:mesh%nv)
       mesh%uv(1:2,1:2,1:mesh%nv) = uvnew(1:2,1:2,1:mesh%nv)
       update_uvxyz : do ivert = 1,mesh%nv ! <----------------------------------+
          select case ( mesh%typ(ivert) ) ! <-------------------------------+   !
          case (0) ! -------------------------------------------------------+   !
             mesh%xyz(:,ivert) = brep%verts(mesh%ids(ivert))%point%xyz      !   !
          case (1) ! -------------------------------------------------------+   !
             iface = brep%edges(mesh%ids(ivert))%halfedges(2)%face          !   !
             call eval( &                                                   !   !
                  mesh%xyz(:,ivert), &                                      !   !
                  brep%faces(iface)%surface, &                              !   !
                  mesh%uv(:,1,ivert) )                                      !   !
          case (2) ! -------------------------------------------------------+   !
             iface = mesh%ids(ivert)                                        !   !
             call eval( &                                                   !   !
                  mesh%xyz(:,ivert), &                                      !   !
                  brep%faces(iface)%surface, &                              !   !
                  mesh%uv(:,1,ivert) )                                      !   !
          end select ! <----------------------------------------------------+   !
       end do update_uvxyz ! <--------------------------------------------------+


       IF ( .TRUE. ) THEN
          WRITE (STRNUM, '(I2.2)') IPASS
          !call write_inria_mesh( &
          !     mesh, &
          !     'Jouke/meshgen/brepmesh/brepmesh_optim_'//strnum//'.mesh' )
          call write_tecplot_mesh( &
               mesh, &
               'Jouke/meshgen/brepmesh/brepmesh_optim_'//strnum//'.dat', &
               'pass_'//strnum )
          !call write_obj_mesh( &
          !     mesh, &
          !     'Jouke/meshgen/brepmesh/brepmesh_optim_'//strnum//'.obj' )
       END IF

       !PAUSE
    end do

  end subroutine optim_jiao













  subroutine energy_jiao( &
       ntr, &
       tri, &
       wei, &
       nve, &
       xyz, &
       hve, &
       frac_conf, &
       ener, &
       grad, &
       hess )
    implicit none
    integer,       intent(in)  :: ntr            ! number of triangles
    integer,       intent(in)  :: tri(3,ntr)     ! triangles (vertex indices)
    real(kind=fp), intent(in)  :: wei(ntr)       ! triangle weights
    integer,       intent(in)  :: nve            ! number of vertices
    real(kind=fp), intent(in)  :: xyz(3,nve)     ! vertices (xyz coordinates)
    real(kind=fp), intent(in)  :: hve(nve)       ! target edge length at vertices
    real(kind=fp), intent(in)  :: frac_conf      ! fraction of conformal energy
    real(kind=fp), intent(out) :: ener(ntr)      ! triangle combined energy
    real(kind=fp), intent(out) :: grad(3,nve)    ! gradient at vertices
    real(kind=fp), intent(out) :: hess(3,3,nve)  ! Hessian matrix at vertices
    real(kind=fp)              :: frac_isom      ! fraction of isometric energy
    real(kind=fp)              :: Id4(3,3)       ! 4*identity(3,3)
    real(kind=fp)              :: nor_tri(3,ntr) ! triangle normals
    real(kind=fp)              :: area_act(ntr)  ! (twice) triangle area (actual)
    real(kind=fp)              :: area_ide(ntr)  ! (twice) triangle area (ideal)
    real(kind=fp)              :: htri           ! target edge length in a triangle
    real(kind=fp)              :: edg_vec(3,3)   ! edge vectors in a triangle
    real(kind=fp)              :: edg_len(3)     ! squared edge lengths in a triangle
    real(kind=fp)              :: edg_90d(3,3)   ! 90 degree CCW rotated edge vectors
    real(kind=fp)              :: ener_conf      ! conformal energy (angle-preserving)
    real(kind=fp)              :: ener_isom      ! isometric energy (area-preserving)
    real(kind=fp)              :: grad_conf(3)   ! gradient of conformal energy
    real(kind=fp)              :: grad_isom(3)   ! gradient of isometric energy
    real(kind=fp)              :: hess_conf(3,3) ! Hessian of conformal energy
    real(kind=fp)              :: hess_isom(3,3) ! Hessian of isometric energy
    integer                    :: itr, ive, i

    frac_isom = 1._fp - frac_conf
    Id4 = 4._fp * identity_matrix(3)

    ! compute triangle normals, actual and target areas
    do itr = 1,ntr
       ! (pseudo-)normal vector
       nor_tri(:,itr) = cross( &
            xyz(:,tri(2,itr)) - xyz(:,tri(1,itr)), &
            xyz(:,tri(3,itr)) - xyz(:,tri(1,itr)) )
       ! actual triangle area
       area_act(itr) = norm2(nor_tri(:,itr))
       nor_tri(:,itr) = nor_tri(:,itr) / area_act(itr) ! unit normal
       ! target edge length for current triangle
       htri = sum(hve(tri(:,itr)))! / 3._fp
       ! ideal triangle area
       area_ide(itr) = htri**2
    end do
    ! normalize ideal triangle areas
    area_ide = area_ide * sum(area_act) / sum(area_ide)


    grad(:,:) = 0._fp
    hess(:,:,:) = 0._fp
    do itr = 1,ntr ! <-------------------------------------------------------------------------------------+
       ! edge vectors                                                                                      !
       edg_vec = xyz(:,tri([3,1,2],itr)) - xyz(:,tri([2,3,1],itr))                                         !
       ! squared edge lengths                                                                              !
       edg_len = sum(edg_vec**2, 1)                                                                        !
       ! 90 degree counter-clockwise rotated edge vectors                                                  !
       do i = 1,3 ! <-----------------------------------------+                                            !
          edg_90d(:,i) = cross(nor_tri(:,itr), edg_vec(:,i))  !                                            !
       end do ! <---------------------------------------------+                                            !
       ! triangle conformal energy                                                                         !
       ener_conf = sum(edg_len)/area_act(itr)                                                              !
       ! triangle isometric energy                                                                         !
       ener_isom = area_act(itr) / area_ide(itr)                                                           !
       ener_isom = ener_isom + 1._fp/ener_isom                                                             !
       ! triangle combined energy                                                                          !
       ener(itr) = frac_conf*ener_conf + frac_isom*ener_isom                                               !
       !                                                                                                   !
       ! vertex-based gradient and Hessian matrix of total energy                                          !
       do i = 1,3 ! <----------------------------------------------------------------------------------+   !
          ive = tri(i,itr)                                                                             !   !
          ! gradient of conformal energy                                                               !   !
          grad_conf = 2._fp*(edg_vec(:,1+mod(i,3)) - edg_vec(:,1+mod(i+1,3))) - ener_conf*edg_90d(:,i) !   !
          grad_conf = grad_conf / area_act(itr)                                                        !   !
          ! gradient of isometric energy                                                               !   !
          grad_isom = (area_act(itr)**2 - area_ide(itr)**2)*edg_90d(:,i)                               !   !
          grad_isom = grad_isom / ( area_ide(itr) * area_act(itr)**2 )                                 !   !
          ! gradient of total energy                                                                   !   !
          grad(:,ive) = grad(:,ive) + wei(itr) * ( frac_conf*grad_conf + frac_isom*grad_isom )         !   !
          !                                                                                            !   !
          ! Hessian matrix of conformal energy                                                         !   !
          hess_conf = outer_product(grad_conf, edg_90d(:,i))                                           !   !
          hess_conf = hess_conf + transpose(hess_conf)                                                 !   !
          hess_conf = (Id4 - hess_conf) / area_act(itr)                                                !   !
          ! Hessian matrix of isometric energy                                                         !   !
          hess_isom = 2._fp * area_ide(itr) * outer_product(edg_90d(:,i), edg_90d(:,i))                !   !
          hess_isom = hess_isom / area_act(itr)**3                                                     !   !
          ! Hessian matrix of combined energy                                                          !   !
          hess(:,:,ive) = hess(:,:,ive) + wei(itr) * ( frac_conf*hess_conf + frac_isom*hess_isom )     !   !
       end do ! <--------------------------------------------------------------------------------------+   !
       !                                                                                                   !
    end do ! <---------------------------------------------------------------------------------------------+

  end subroutine energy_jiao


end module mod_optimmesh
