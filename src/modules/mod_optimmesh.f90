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
       !PRINT *,'OUTER, K=',k
       !IF ( K > 10 ) EXIT
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
                !PAUSE
                !IF ( K > 9 ) EXIT OUTER_LOOP
                call edge_contraction( &
                     brep, &
                     hyperedges, &
                     nhe, &
                     mesh, &
                     [iedge,iface] )
                
                !call write_connectivity( &
                !     '../debug/', &
                !     mesh, &
                !     1 )
                !call write_xyz_positions( &
                !     '../debug/', &
                !     mesh, &
                !     1 )
                !IF ( ANY(TYPIJ < 2) ) PAUSE
                !IF ( ANY(TYPIJ < 2) ) EXIT OUTER_LOOP
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
    integer                                :: verts(2), faces(2), i, j, k, ivar
    integer                                :: ikp, jkp, irm, jrm, inew, imx
    integer                                :: it, iv, ip, ihype
    integer                                :: v_old2new(mesh%nv), t_old2new(mesh%nt)
    logical :: debugproj

    verts = mesh%tri([ihedg(1), 1+mod(ihedg(1),3)],ihedg(2))
    PRINT *,'VERTS =',VERTS
    
    !PRINT *,'XYZ ='
    !CALL PRINT_MAT(TRANSPOSE(MESH%XYZ(1:3,VERTS)))

    ! geometry
    if ( mesh%typ(verts(1)) == mesh%typ(verts(2)) ) then
       ! midpoint
       jkp = minloc(verts, 1)
       ikp = verts(jkp)
       typnew = mesh%typ(verts(1))
       if ( mesh%ids(verts(1)) == mesh%ids(verts(2)) ) then
          uvnew  = 0.5 * sum(mesh%uv(:,:,verts), 3)
          idsnew = mesh%ids(verts(1))
          if ( typnew == 1 ) then
             dxyz = 0.5_fp * (mesh%xyz(:,verts(2)) - mesh%xyz(:,verts(1)))
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
                  .true., &
                  stat )
             call eval( &
                  xyznew, &
                  brep%edges(idsnew)%curve%surf(1)%ptr, &
                  uvnew(:,1) )
             !call eval( &
             !     xyznew, &
             !     brep%edges(idsnew)%curve%surf(2)%ptr, &
             !     uvnew(:,1) )
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
                  .true., &
                  stat )
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
             debugproj = .false.
             call projection_hyperface( &
                  brep, &
                  mesh%ids(verts(1)), &
                  mesh%uv(:,1,verts(1)), &
                  mesh%xyz(:,verts(1)), &
                  duv(:,1), &
                  dxyz, &
                  idsnew, &
                  uvnew(:,1), &
                  debugproj, &
                  stat )
             !PRINT *,'STAT =',STAT
             if ( stat > 0 ) PAUSE
             call eval( &
                  xyznew, &
                  brep%faces(idsnew)%surface, &
                  uvnew(:,1) )
          end if
       end if
    else
       jkp = minloc(mesh%typ(verts),1)
       ikp = verts(jkp)
       ! contract towards ikp
       xyznew = mesh%xyz(:,ikp)
       uvnew  = mesh%uv(:,:,ikp)
       idsnew = mesh%ids(ikp)
       typnew = mesh%typ(ikp)
    end if

    jrm = 1+mod(jkp,2)
    irm = verts(jrm)

    PRINT *,'    IDS =',mesh%ids(verts)
    PRINT *,'    TYP =',mesh%typ(verts)

    inew = minval(verts)
    imx = maxval(verts)

    v_old2new(1:imx-1) = [(j, j=1,imx-1)]
    v_old2new(imx) = inew
    v_old2new(imx+1:mesh%nv) = [(j, j=imx,mesh%nv-1)]
    !v_old2new([irm,ikp]) = inew
 
    if ( mesh%typ(irm) == 1 ) then
       ihype = brep%edges(mesh%ids(irm))%hyperedge
       do ip = 1,mesh%npaths
          if ( mesh%paths(ip)%hyperedge == ihype ) then
             call remove_from_list( &
                  imx, &
                  mesh%paths(ip)%verts, &
                  mesh%paths(ip)%nv )
             !do iv = 2,mesh%paths(ip)%nv-1
             !   if ( mesh%paths(ip)%verts(iv) == irm ) exit
             !end do
             !mesh%paths(ip)%verts(iv:mesh%paths(ip)%nv-1) = mesh%paths(ip)%verts(iv+1:mesh%paths(ip)%nv)
             !mesh%paths(ip)%nv = mesh%paths(ip)%nv - 1
             exit
          end if
       end do
    end if

    ! new topology
    mesh%ids(inew) = idsnew
    mesh%typ(inew) = typnew
    mesh%xyz(:,inew) = xyznew
    mesh%uv(:,:,inew) = uvnew

    PRINT *,'   I_RM =',irm
    PRINT *,'   I_KP =',ikp
    PRINT *,'  I_NEW =',inew
    PRINT *,'XYZ_NEW =',xyznew
    PRINT *,' UV_NEW =',uvnew
    PRINT *,'TYP_NEW =',typnew
    PRINT *,'IDS_NEW =',idsnew
    
    ! remove vertex
    mesh%typ(imx:mesh%nv-1) = mesh%typ(imx+1:mesh%nv)
    mesh%ids(imx:mesh%nv-1) = mesh%ids(imx+1:mesh%nv)
    mesh%xyz(1:3,imx:mesh%nv-1) = mesh%xyz(1:3,imx+1:mesh%nv)
    mesh%uv(1:2,1:2,imx:mesh%nv-1) = mesh%uv(1:2,1:2,imx+1:mesh%nv)
    !mesh%v2h(1:2,i:mesh%nv-1) = mesh%v2h(1:2,i+1:mesh%nv)
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
    mesh%ihf(i:j-2) = mesh%ihf(i+1:j-1)
    mesh%ihf(j-1:mesh%nt-2) = mesh%ihf(j+1:mesh%nt)
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
    use mod_tolerances
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
             w = sum((mesh%xyz(:,jvert) - mesh%xyz(:,ivert))**2)
             if ( w > EPSxyzsqr ) then
                n = n + 1
                w = 1._fp / sqrt(w)
                tmp(ivert) = tmp(ivert) + w * fv(jvert)
                sumw = sumw + w
             else
                PRINT *,'/!\ ADJACENT VERTICES TOO CLOSE, VERTS =',IVERT,JVERT,', DISTANCE =',SQRT(W)
             end if

             ihedg = get_prev(ihedg)
             ihedg = get_twin(mesh, ihedg)
             jface = get_face(ihedg)
             if ( jface < 1 .or. jface == iface ) exit adjacent_verts
          end do adjacent_verts
          if ( sumw < EPSfp ) then
             w = 1._fp
          else
             w = sumw / real(n, kind=fp)
          end if
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
    ! compute pseudo-normals at vertices
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
          jv = mesh%paths(ip)%verts(lv)
          if ( mesh%typ(jv) == 0 ) cycle
          iv = mesh%paths(ip)%verts(1 + mod(lv + mesh%paths(ip)%nv - 2,mesh%paths(ip)%nv))
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
       frac_conf1, &
       frac_conf2, &
       ipass1, &
       ipass2, &
       passmax, &
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
    real(kind=fp), parameter               :: EPSdxyz = 5.d-2 ! *min(h)
    real(kind=fp), parameter               :: EPSdxyzsqr = EPSdxyz**2
    type(type_brep),         intent(in)    :: brep
    integer,                 intent(in)    :: nhe
    type(type_hyperedge),    intent(in)    :: hyperedges(nhe)
    type(type_surface_mesh), intent(inout) :: mesh
    real(kind=fp),           intent(in)    :: frac_conf1
    real(kind=fp),           intent(in)    :: frac_conf2
    integer,                 intent(in)    :: ipass1
    integer,                 intent(in)    :: ipass2
    integer,                 intent(in)    :: passmax
    real(kind=fp),           intent(in)    :: hmin
    real(kind=fp),           intent(in)    :: hmax
    real(kind=fp)                          :: hve(mesh%nv), rc(mesh%nv)
    real(kind=fp)                          :: wei(mesh%nt)
    real(kind=fp)                          :: hminsqr
    real(kind=fp)                          :: ener(mesh%nt)
    real(kind=fp)                          :: grad(3,mesh%nv)
    real(kind=fp)                          :: hess(3,3,mesh%nv)
    real(kind=fp)                          :: duv(2,2,mesh%nv), ds
    real(kind=fp)                          :: tng(3,2), duv_ds(2,2,2)
    real(kind=fp)                          :: xyzproj(3)
    integer                                :: stat
    logical                                :: singular
    real(kind=fp), dimension(4)            :: lowerb, upperb
    real(kind=fp)                          :: uvtmp(2,2), dxyz(3), maxdxyz, maxdxyz_rc
    integer, dimension(mesh%nv)            :: idsnew, typnew
    real(kind=fp)                          :: uvnew(2,2,mesh%nv)
    integer                                :: ivert, jvert, ihedg(2)
    integer                                :: iedge, ihype
    integer                                :: iface, jface
    integer                                :: ipass, ivar
    logical                                :: check_change
    real(kind=fp)                          :: frac_conf_ramp
    real(kind=fp)                          :: etot, etot0, etotprev
    real(kind=fp)                          :: detot, detot0
    CHARACTER(2) :: STRNUM
    INTEGER :: FID, I, J, subit
    REAL(KIND=FP) :: DXYZVISU(3,MESH%NV)
    logical :: debugproj

    CALL GET_FREE_UNIT(FID)

    lowerb(:) = -2._fp
    upperb(:) = 2._fp

    ! compute triangle weights
    IF ( .false. ) THEN
       wei(:) = 1._fp
    ELSE
       call compute_triangle_weights( &
            mesh, &
            2._fp, &
            1._fp, &
            wei )
    END IF

    do ipass = 1,passmax
       PRINT *,''
       PRINT *,'OPTIM PASS',IPASS,'/',PASSMAX
       ! compute target edge lengths at vertices
       if ( .false. ) then!ipass < 6 ) then
          hve(:) = 1._fp
       else
          IF ( .false. ) THEN
             hve(:) = huge(1._fp)
             do ivert = 1,mesh%nv
                if ( mesh%typ(ivert) /= 2 ) cycle
                call eval_minimum_curvature_radius( &
                     brep%faces(mesh%ids(ivert))%surface, &
                     mesh%uv(:,1,ivert), &
                     hve(ivert) )
             end do
          ELSE
             call discrete_minimum_curvature_radius( &
                  mesh, &
                  hve )
          END IF
          rc = 0.25_fp * hve**2
          !hve = min(minval(hve)*hmax/hmin, hve)
          !PRINT *,''!MIN(HVE) =', real(minval(hve)), minloc(hve,1)
          hve = hve/minval(hve)
          hve = min(hmax/hmin, hve)
          !hve = min(3._fp, hve)
          ! normalize to [0,1] for plotting
          !hve = (hve - minval(hve)) / (maxval(hve) - minval(hve))

          call smooth_function_mesh( &
               mesh, &
               hve, &
               30 )
       end if

       IF ( .false. ) THEN
          WRITE (STRNUM, '(I2.2)') IPASS
          call write_tecplot_mesh_solv( &
               mesh, &
               '../test/optimmesh/optim_'//strnum//'.dat', &
               'pass_'//strnum, &
               hve )
          call write_vtk_mesh( &
               mesh, &
               '../test/optimmesh/optim_'//strnum//'.vtk' )
       END IF

       
       !! compute triangle weights
       !!wei(:) = 1._fp
       !call compute_triangle_weights( &
       !     mesh, &
       !     5._fp, &
       !     1._fp, &
       !     wei )
       
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
            hminsqr, &
            ener, &
            grad, &
            hess )

       if ( .false. ) then !passmax < 10 ) then
          open(unit=13, file='../tmp/xyz.dat', action='write')
          do i = 1,mesh%nv
             write (13,*) mesh%xyz(:,i)
          end do
          close(13)

          open(unit=13, file='../tmp/ener.dat', action='write')
          do i = 1,mesh%nt
             write (13,*) ener(i)
          end do
          close(13)

          open(unit=13, file='../tmp/grad.dat', action='write')
          do i = 1,mesh%nv
             write (13,*) grad(:,i)
          end do
          close(13)

          open(unit=13, file='../tmp/hess.dat', action='write')
          do i = 1,mesh%nv
             write (13,*) hess(:,:,i)
          end do
          close(13)
       end if
       !
       etot = sum(ener(1:mesh%nt))
       if ( ipass == 1 ) etot0 = etot
       etot = etot/etot0
       PRINT *,'E/E0 =',real(etot)
       !
       ! compute vertex displacements
       !OPEN(UNIT=FID, FILE='Jouke/meshgen/brepmesh/debug_dxyz.dat', ACTION='WRITE')
       maxdxyz = 0._fp
       maxdxyz_rc = 0._fp
       compute_duv : do ivert = 1,mesh%nv
          !PRINT *,'IVERT =', IVERT
          select case ( mesh%typ(ivert) ) ! <-------------------------------+
          case (0) ! -------------------------------------------------------+
             uvnew(:,:,ivert) = mesh%uv(:,:,ivert)                          !
             idsnew(ivert) = mesh%ids(ivert)
             typnew(ivert) = 0
             DXYZ(:) = 0._FP
             DXYZVISU(:,IVERT) = 0._FP
          case (1) ! -------------------------------------------------------+
             call diffgeom_intersection( &                                  !
                  brep%edges(mesh%ids(ivert))%curve%surf, &                 !
                  mesh%uv(:,:,ivert), &                                     !
                  duv_ds, &                                                 !
                  tng, &                                                    !
                  stat )                                                    !
             IF ( .false. ) THEN
                DS = 0._FP
             ELSE
                ds = - dot_product(tng(:,1), grad(:,ivert)) / &                !
                     dot_product(tng(:,1), matmul(hess(:,:,ivert), tng(:,1)))  !
             END IF

             do subit = 1,7
                duv(:,:,ivert) = ds * duv_ds(:,1,:)                            !
                dxyz = ds * tng(:,1)
                maxdxyz = max(maxdxyz, sum(dxyz**2))
                maxdxyz_rc = max(maxdxyz_rc, sum(dxyz**2)/rc(ivert))
                DXYZVISU(:,IVERT) = DXYZ
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
                     .false., &!(ivert == 1666) ,&!
                     stat )
                if ( stat > 0 ) THEN
                   PRINT *,'IVERT =',IVERT
                   PRINT *,'DXYZ =',DXYZ
                   ds = 0.5_fp * ds
                   !PAUSE
                else
                   exit
                END if
             end do
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
             maxdxyz = max(maxdxyz, sum(dxyz**2))
             maxdxyz_rc = max(maxdxyz_rc, sum(dxyz**2)/rc(ivert))
             DXYZVISU(:,IVERT) = DXYZ
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
                !debugproj = (ivert == 23977 .and. passmax<50 .and. &
                !     mesh%uv(1,1,ivert) < -0.738_fp  .and. mesh%uv(2,1,ivert) < -0.985_fp)!.false.!
                debugproj = .false.
                call projection_hyperface( &
                     brep, &
                     iface, &
                     mesh%uv(:,1,ivert), &
                     mesh%xyz(:,ivert), &
                     duv(:,1,ivert), &
                     dxyz, &
                     idsnew(ivert), &
                     uvtmp(:,1), &
                     debugproj, &
                     stat )
                if ( stat > 0 ) then
                   PRINT *,'IVERT =',IVERT
                   PRINT *,'DXYZ =',DXYZ
                   uvtmp(1:2,1) = mesh%uv(1:2,1,ivert)
                   call projection_surface( &
                        brep%faces(mesh%ids(ivert))%surface, &
                        mesh%xyz(1:3,ivert) + dxyz, &
                        uvtmp(1:2,1), &
                        [-2._fp, -2._fp], &
                        [2._fp, 2._fp], &
                        stat, &
                        xyzproj )
                   PRINT *,'UVTMP =',UVTMP(1:2,1)
                   PRINT *,'XYZPROJ =',XYZPROJ
                end if
                if ( stat > 0 ) THEN
                   PRINT *,'IVERT =',IVERT
                   PRINT *,'DXYZ =',DXYZ
                   PAUSE
                END if
                uvnew(:,1,ivert) = uvtmp(:,1)
             end if
             !
          end select ! <----------------------------------------------------+
          !
          IF ( MAXVAL(ABS(UVNEW(:,:,IVERT))) > 1._FP + EPSUV ) PRINT *, IVERT, TYPNEW(IVERT), IDSNEW(IVERT), UVNEW(:,:,IVERT)
          !WRITE (FID,*) DXYZ
       end do compute_duv
       !CLOSE(FID)

       PRINT *,'MAX(DXYZ) =',SQRT(MAXDXYZ)
       PRINT *,'MAX(DXYZ/RC) =',SQRT(MAXDXYZ_RC)
       
       if ( .false. ) then!passmax < 10 ) then
          CALL write_tecplot_mesh_displacement2( &
               mesh, &
               '/d/bandrieu/GitHub/FFTsurf/cases/vortex/output/debug/visu.dat', &
               'visu_disp', &
               dxyzvisu )
       end if

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


       IF ( .false. ) THEN
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


          open( &
               unit = fid, &
               file = 'Jouke/meshgen/brepmesh/xyz_optim_'//strnum//'.dat', &
               action = 'write' )
          do j = 1,mesh%nv
             do i = 1,3
                write (fid, '(e22.15,1x)', advance='no') mesh%xyz(i,j)
             end do
             write (fid,*)
          end do
          close(fid)
       END IF

       if ( .false. ) PAUSE

       if ( maxdxyz/hminsqr < EPSdxyzsqr ) then
          PRINT *,'MAX(DXYZ) << MIN(H)'
          exit
       end if
       
       if ( ipass > 1 ) then
          detot = etotprev - etot
          if ( ipass > 2 ) then
             if ( detot0 < 0._fp ) detot0 = detot
             PRINT *,'DETOT/DETOT0 =',real(detot/detot0)
             if ( detot < 0.01_fp * detot0 ) then ! .or. abs(detot) < 0.005 ) then
                IF ( DETOT < 0._FP ) THEN
                   PRINT *,'NON DECROISSANT'
                ELSE
                   PRINT *,'STATIONNAIRE'
                END IF
                exit
             end if
          else
             detot0 = detot
             PRINT *,'DETOT0 =',real(detot0)
          end if
       end if
       etotprev = etot
       
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
       hminsqr, &
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
    real(kind=fp), intent(out) :: hminsqr        ! smallest edge length (squared)
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
       htri = sum(hve(tri(:,itr)))
       ! ideal triangle area
       area_ide(itr) = htri**2
    end do
    ! normalize ideal triangle areas
    area_ide = area_ide * sum(area_act) / sum(area_ide)
    PRINT *,'MINMAX AREA_ACT',MINVAL(AREA_ACT), MAXVAL(AREA_ACT)
    PRINT *,'MINMAX AREA_IDE',MINVAL(AREA_IDE), MAXVAL(AREA_IDE)


    grad(:,:) = 0._fp
    hess(:,:,:) = 0._fp
    hminsqr = huge(1._fp)
    do itr = 1,ntr ! <-------------------------------------------------------------------------------------+
       ! edge vectors                                                                                      !
       edg_vec = xyz(:,tri([3,1,2],itr)) - xyz(:,tri([2,3,1],itr))                                         !
       ! squared edge lengths                                                                              !
       edg_len = sum(edg_vec**2, 1)                                                                        !
       hminsqr = min(hminsqr, minval(edg_len))                                                             !
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

    PRINT *,'HMIN =',SQRT(HMINsqr)
    
  end subroutine energy_jiao













  subroutine write_tecplot_mesh_displacement( &
       mesh, &
       filename, &
       zonename, &
       dxyz, &
       ener )
    use mod_util
    use mod_mesh
    implicit none
    type(type_surface_mesh), intent(in) :: mesh
    character(*),            intent(in) :: filename
    character(*),            intent(in) :: zonename
    real(kind=fp),           intent(in) :: dxyz(3,mesh%nv)
    real(kind=fp),           intent(in) :: ener(mesh%nt)
    integer                             :: fid, i, j

    call get_free_unit(fid)
    open(unit=fid, file=filename, action='write')

    write (fid,*) 'VARIABLES = "X" "Y" "Z" "U" "V" "W" "E"'

    write (fid,*) 'ZONE T="' // zonename // '"'
    write (fid,'(A2,I0,A3,I0)') 'N=',mesh%nv,' E=',mesh%nt
    write (fid,*) 'ZONETYPE=FETriangle'
    write (fid,*) 'DATAPACKING=BLOCK'

    write (fid,*) 'VARLOCATION = ([1-6]=NODAL, [7]=CELLCENTERED)'

    do i = 1,3
       write (fid,*) ''
       do j = 1,mesh%nv
          write (fid,'(ES15.7)') mesh%xyz(i,j)
       end do
    end do

    do i = 1,3
       write (fid,*) ''
       do j = 1,mesh%nv
          write (fid,'(ES15.7)') dxyz(i,j)
       end do
    end do

    write (fid,*) ''

    do j = 1,mesh%nt
       write (fid,'(ES15.7)') ener(j)
    end do

    write (fid,*) ''

    do j = 1,mesh%nt
       write (fid,'(I0,1x,I0,1x,I0)') mesh%tri(:,j)
    end do

    close(fid)

  end subroutine write_tecplot_mesh_displacement











subroutine write_tecplot_mesh_displacement2( &
       mesh, &
       filename, &
       zonename, &
       dxyz )
    use mod_util
    use mod_mesh
    implicit none
    type(type_surface_mesh), intent(in) :: mesh
    character(*),            intent(in) :: filename
    character(*),            intent(in) :: zonename
    real(kind=fp),           intent(in) :: dxyz(3,mesh%nv)
    integer                             :: fid, j

    call get_free_unit(fid)
    open(unit=fid, file=filename, action='write')

    write (fid,*) 'VARIABLES = "X" "Y" "Z" "U" "V" "W"'

    write (fid,*) 'ZONE T="' // zonename // '"'
    write (fid,*) 'STRANDID=0, SOLUTIONTIME=0'
    write (fid,*) 'Nodes=',mesh%nv,', Elements=',mesh%nt
    write (fid,*) 'ZONETYPE=FETriangle'
    write (fid,*) 'DATAPACKING=POINT'
    
    do j = 1,mesh%nv
       write (fid,*) mesh%xyz(:,j), dxyz(:,j)
    end do

    do j = 1,mesh%nt
       write (fid,'(I0,1x,I0,1x,I0)') mesh%tri(:,j)
    end do

    close(fid)

  end subroutine write_tecplot_mesh_displacement2
  

  







  subroutine compute_triangle_weights( &
       mesh, &
       wclose, &
       wfar, &
       wt )
    use mod_mesh
    use mod_halfedge
    implicit none
    real(kind=fp), parameter             :: EPS = 0.001_fp, gam = -log(EPS), fhavg = 5._fp
    type(type_surface_mesh), intent(in)  :: mesh
    real(kind=fp),           intent(in)  :: wclose
    real(kind=fp),           intent(in)  :: wfar
    real(kind=fp),           intent(out) :: wt(mesh%nt)
    real(kind=fp)                        :: hmin, hmax, alp, bet
    real(kind=fp)                        :: dv(mesh%nv), dij, dt
    logical                              :: visited(mesh%nv)
    integer, dimension(mesh%nv)          :: front, fronttmp
    integer                              :: nfront, nfronttmp
    integer                              :: ivert, jvert, kvert, ihedg(2), iface, jface

    visited(1:mesh%nv) = .false.
    ! initial front
    nfront = 0
    do ivert = 1,mesh%nv
       if ( mesh%typ(ivert) < 2 ) then
          nfront = nfront + 1
          front(nfront) = ivert
       end if
    end do
    if ( nfront < 1 ) then
       wt(1:mesh%nt) = 1._fp
       return
    end if
    dv(:) = huge(1._fp)
    dv(front(1:nfront)) = 0._fp
    visited(front(1:nfront)) = .true.

    hmin = huge(1._fp)
    hmax = 0._fp
    do while ( .not.all(visited) )
       nfronttmp = 0
       do kvert = 1,nfront
          ivert = front(kvert)
          ihedg = mesh%v2h(:,ivert)
          iface = get_face(ihedg)
          adjacent_verts : do
             jvert = get_dest(mesh, ihedg)
             dij = norm2(mesh%xyz(:,ivert) - mesh%xyz(:,jvert))
             hmin = min(hmin, dij)
             hmax = max(hmax, dij)
             dv(jvert) = min(dv(jvert), dv(ivert) + dij)
             if ( .not.visited(jvert) ) then
                nfronttmp = nfronttmp + 1
                fronttmp(nfronttmp) = jvert
                visited(jvert) = .true.
             end if
             ihedg = get_prev(ihedg)
             ihedg = get_twin(mesh, ihedg)
             jface = get_face(ihedg)
             if ( jface < 1 .or. jface == iface ) exit adjacent_verts
          end do adjacent_verts
       end do

       nfront = nfronttmp
       front(1:nfront) = fronttmp(1:nfront)
    end do

    bet = (wclose - wfar)/(1._fp - EPS)
    alp = wclose - bet
    
    dv = dv / maxval(dv)
    do iface = 1,mesh%nt
       dt = sum(dv(mesh%tri(:,iface))) / 3._fp
       !wt(iface) = wclose + (wfar - wclose)*dt
       wt(iface) = alp + bet*exp(-gam*dt*2._fp/(fhavg*(hmin + hmax)))
    end do

    
  end subroutine compute_triangle_weights






  subroutine pre_deformation( &
       brep, &
       mesh, &
       dxyz )
    use mod_types_brep
    use mod_mesh
    use mod_diffgeom
    use mod_linalg
    use mod_halfedge
    use mod_projection
    use mod_tolerances
    implicit none
    real(kind=fp), parameter               :: EPShmin = 1.d-3
    type(type_brep),         intent(in)    :: brep
    type(type_surface_mesh), intent(inout) :: mesh
    real(kind=fp),           intent(inout) :: dxyz(3,mesh%nv)
    integer, dimension(mesh%nv)            :: idsnew
    real(kind=fp)                          :: duv(2), uvnew(2), xyzproj(3)
    real(kind=fp)                          :: dxyz_duv(3,2)
    real(kind=fp)                          :: hmin
    logical                                :: singular
    logical                                :: check_change
    integer                                :: stat
    integer                                :: ivert, jvert, ivar, ihedg(2), iface, jface
    logical :: debugproj

    PRINT *,'PRE_DEFORMATION'
    debugproj = .false.
    
    idsnew(1:mesh%nv) = mesh%ids(1:mesh%nv)
    vertices : do ivert = 1,mesh%nv
       if ( mesh%typ(ivert) < 2 ) cycle
       do ivar = 1,2 ! <-------------------------------+
          call evald1( &                               !
               dxyz_duv(:,ivar), &                     !
               brep%faces(mesh%ids(ivert))%surface, &  !
               mesh%uv(:,1,ivert), &                   !
               ivar )                                  !
       end do ! <--------------------------------------+
       call solve_NxN( &
            duv, &
            matmul(transpose(dxyz_duv), dxyz_duv), &
            matmul(transpose(dxyz_duv), dxyz(:,ivert)), &
            singular )
       !
       dxyz(1:3,ivert) = matmul(dxyz_duv, duv)
       ! handle passing to an adjacent face (crossing of a smooth edge...)
       uvnew = mesh%uv(:,1,ivert) + duv
       ihedg = mesh%v2h(:,ivert) ! mesh halfedge index
       iface = mesh%ids(ivert)   ! brep face index
       check_change = ( maxval(abs(duv)) > 1._fp + EPSuv - maxval(abs(mesh%uv(:,1,ivert))) )
       if ( .not.check_change ) then
          hmin = huge(1._fp)
          adjacent_verts : do ! <------------------------------------------------------+
             jvert = get_dest(mesh, ihedg) ! mesh vertex index                         !
             jface = mesh%ids(jvert)       ! brep face index                           !
             hmin = min(hmin, sum((mesh%xyz(:,jvert) - mesh%xyz(:,ivert))**2))
             ! (quasi) necessary conditions for a change of supporting brep face:      !
             ! - at least one adjacent vertex is supported by a different brep face;   !
             ! - this vertex is in the halfspace pointed by the xyz displacement.      !
             if ( mesh%typ(jvert) /= 2 .or. jface /= iface ) then ! <---------------+  !
                if ( dot_product(dxyz(1:3,ivert), &                                 !  !
                     mesh%xyz(:,jvert) - mesh%xyz(:,ivert)) > 0._fp ) then ! <---+  !  !
                   check_change = .true.                                         !  !  !
                   !exit adjacent_verts                                           !  !  !
                end if ! <-------------------------------------------------------+  !  !
             end if ! <-------------------------------------------------------------+  !
             ! move on to next adjacent vertex                                         !
             ihedg = get_prev(ihedg)       ! outgoing mesh halfedge                    !
             ihedg = get_twin(mesh, ihedg) ! ingoing mesh halfedge                     !
             if ( ihedg(2) < 1 .or. ihedg(2) ==  mesh%v2h(2,ivert) ) exit              !
          end do adjacent_verts ! <----------------------------------------------------+
          if ( sum(dxyz(1:3,ivert)**2) < hmin*EPShmin**2 ) then
             !PRINT *,'IVERT =', IVERT, ', |DXYZ| =',NORM2(dxyz(1:3,ivert)),', HMIN =',SQRT(HMIN)
             dxyz(1:3,ivert) = 0._fp
             cycle vertices
          end if
       end if
       !
       if ( check_change ) then
          !debugproj = (ivert == 21116)
          !debugproj = (ivert == 23977 .and. &
          !           mesh%uv(1,1,ivert) < -0.739_fp  .and. mesh%uv(2,1,ivert) < -0.986_fp)
          call projection_hyperface( &
               brep, &
               iface, &
               mesh%uv(1:2,1,ivert), &
               mesh%xyz(1:3,ivert), &
               duv, &
               dxyz(1:3,ivert), &
               idsnew(ivert), &
               uvnew, &
               debugproj, &
               stat )
          if ( stat > 0 ) then
             PRINT *,'IVERT =',IVERT
             PRINT *,'DXYZ =',DXYZ(:,ivert)
             uvnew = mesh%uv(1:2,1,ivert)
             call projection_surface( &
                  brep%faces(mesh%ids(ivert))%surface, &
                  mesh%xyz(1:3,ivert) + dxyz(1:3,ivert), &
                  uvnew, &
                  [-2._fp, -2._fp], &
                  [2._fp, 2._fp], &
                  stat, &
                  xyzproj )
             idsnew(ivert) = mesh%ids(ivert)
             PRINT *,'UVNEW =',UVnew
             PRINT *,'XYZPROJ =',XYZPROJ
          end if
          if ( stat > 0 ) then
             PRINT *,'IVERT =',IVERT
             PRINT *,'DXYZ =',DXYZ(:,ivert)
             PAUSE
          else
             mesh%uv(1:2,1,ivert) = uvnew
          end if
       else
          mesh%uv(1:2,1,ivert) = uvnew
       end if
    end do vertices

    mesh%ids(1:mesh%nv) = idsnew(1:mesh%nv)

    do ivert = 1,mesh%nv
       if ( mesh%typ(ivert) < 2 ) cycle
       call eval( &
            mesh%xyz(1:3,ivert), &
            brep%faces(idsnew(ivert))%surface, &
            mesh%uv(1:2,1,ivert) )
    end do
    
  end subroutine pre_deformation


















  subroutine spring_displacement_smoothing( &
       mesh, &
       dxyz, &
       npass )
    use mod_mesh
    use mod_halfedge
    use mod_tolerances
    implicit none
    type(type_surface_mesh), intent(in)    :: mesh
    real(kind=fp),           intent(inout) :: dxyz(3,mesh%nv)
    integer,                 intent(in)    :: npass
    real(kind=fp)                          :: maxd0, dxyztmp(3,mesh%nv), w, sumw
    integer                                :: ipass, ivert, jvert, ihedg(2), iface, jface, n

    maxd0 = sqrt(maxval(sum(dxyz**2,1)))
    PRINT *,'MAXD0 =',MAXD0
    
    do ipass = 1,npass
       dxyztmp(1:3,1:mesh%nv) = 0._fp
       do ivert = 1,mesh%nv
          if ( mesh%typ(ivert) == 2 ) then
             ihedg = mesh%v2h(:,ivert)
             iface = get_face(ihedg)
             sumw = 0._fp
             n = 0
             adjacent_verts : do
                jvert = get_dest(mesh, ihedg)
                w = sum((mesh%xyz(:,jvert) - mesh%xyz(:,ivert))**2)
                if ( w > EPSxyzsqr ) then
                   n = n + 1
                   w = 1._fp / sqrt(w)
                   sumw = sumw + w
                   dxyztmp(1:3,ivert) = dxyztmp(1:3,ivert) + w*dxyz(1:3,jvert)
                else
                   PRINT *,'/!\ ADJACENT VERTICES TOO CLOSE, VERTS =',IVERT,JVERT,', DISTANCE =',SQRT(W)
                end if

                ihedg = get_prev(ihedg)
                ihedg = get_twin(mesh, ihedg)
                jface = get_face(ihedg)
                if ( jface < 1 .or. jface == iface ) exit adjacent_verts
             end do adjacent_verts
             if ( sumw < EPSfp ) then
                w = 1._fp
             else
                w = sumw / real(n, kind=fp)
             end if
             dxyztmp(1:3,ivert) = dxyztmp(1:3,ivert) + w * dxyz(1:3,ivert)
             sumw = sumw + w
             dxyztmp(1:3,ivert) = dxyztmp(1:3,ivert) / sumw
          else
             dxyztmp(1:3,ivert) = dxyz(1:3,ivert)
          end if
       end do
       PRINT *,'PASS #',IPASS,', DELTA* =',SQRT(MAXVAL(SUM((DXYZ - DXYZTMP)**2,1))) / maxd0
       
       dxyz(1:3,1:mesh%nv) = dxyztmp(1:3,1:mesh%nv)
    end do

    !PRINT *,'MAX =',SQRT(MAXVAL(SUM(DXYZ(MESH%TYP(1:MESH%NV)
    MAXD0 = 0._FP
    DO IVERT = 1,MESH%NV
       IF ( MESH%TYP(IVERT) == 2 ) MAXD0 = MAX(MAXD0, SUM(DXYZ(:,IVERT)**2))
    END DO
    PRINT *,'MAX DXYZ_FREE =',SQRT(MAXD0)

  end subroutine spring_displacement_smoothing













  subroutine optim_jiao_uv( &
       brep, &
       mesh, &
       frac_conf1, &
       frac_conf2, &
       ipass1, &
       ipass2, &
       passmax, &
       hmin, &
       hmax )
    USE MOD_UTIL
    use mod_diffgeom
    use mod_intersection
    use mod_types_brep
    use mod_brep
    use mod_mesh
    use mod_halfedge
    use mod_linalg
    use mod_projection
    USE MOD_TOLERANCES
    implicit none
    real(kind=fp), parameter               :: EPSdxyz = 5.d-2 ! *min(h)
    real(kind=fp), parameter               :: EPSdxyzsqr = EPSdxyz**2
    type(type_brep),         intent(in)    :: brep
    type(type_surface_mesh), intent(inout) :: mesh
    real(kind=fp),           intent(in)    :: frac_conf1
    real(kind=fp),           intent(in)    :: frac_conf2
    integer,                 intent(in)    :: ipass1
    integer,                 intent(in)    :: ipass2
    integer,                 intent(in)    :: passmax
    real(kind=fp),           intent(in)    :: hmin
    real(kind=fp),           intent(in)    :: hmax
    real(kind=fp)                          :: hve(mesh%nv), rc(mesh%nv)
    real(kind=fp)                          :: wei(mesh%nt)
    real(kind=fp)                          :: hminsqr
    real(kind=fp)                          :: ener(mesh%nt)
    real(kind=fp)                          :: grad(3,mesh%nv)
    real(kind=fp)                          :: hess(3,3,mesh%nv)
    real(kind=fp)                          :: duv(2,mesh%nv), dxyz(3)
    real(kind=fp)                          :: tng(3,2)
    logical                                :: singular
    integer                                :: ivert, iface
    integer                                :: ipass, ivar
    real(kind=fp)                          :: maxdxyz, maxdxyz_rc
    real(kind=fp)                          :: frac_conf_ramp
    real(kind=fp)                          :: etot, etot0, etotprev
    real(kind=fp)                          :: detot, detot0

    ! compute triangle weights
    IF ( .false. ) THEN
       wei(:) = 1._fp
    ELSE
       call compute_triangle_weights( &
            mesh, &
            2._fp, &
            1._fp, &
            wei )
    END IF

    do ipass = 1,passmax
       PRINT *,''
       PRINT *,'OPTIM PASS',IPASS,'/',PASSMAX
       ! compute target edge lengths at vertices
       if ( .false. ) then!ipass < 6 ) then
          hve(:) = 1._fp
       else
          IF ( .false. ) THEN
             hve(:) = huge(1._fp)
             do ivert = 1,mesh%nv
                if ( mesh%typ(ivert) /= 2 ) cycle
                call eval_minimum_curvature_radius( &
                     brep%faces(mesh%ids(ivert))%surface, &
                     mesh%uv(:,1,ivert), &
                     hve(ivert) )
             end do
          ELSE
             call discrete_minimum_curvature_radius( &
                  mesh, &
                  hve )
          END IF
          rc = 0.25_fp * hve**2
          !hve = min(minval(hve)*hmax/hmin, hve)
          !PRINT *,''!MIN(HVE) =', real(minval(hve)), minloc(hve,1)
          hve = hve/minval(hve)
          hve = min(hmax/hmin, hve)
          !hve = min(3._fp, hve)
          ! normalize to [0,1] for plotting
          !hve = (hve - minval(hve)) / (maxval(hve) - minval(hve))

          call smooth_function_mesh( &
               mesh, &
               hve, &
               30 )
       end if

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
            hminsqr, &
            ener, &
            grad, &
            hess )

       etot = sum(ener(1:mesh%nt))
       if ( ipass == 1 ) etot0 = etot
       etot = etot/etot0
       PRINT *,'E/E0 =',real(etot)
       
       ! compute vertex displacements
       maxdxyz = 0._fp
       maxdxyz_rc = 0._fp
       duv(1:2,1:mesh%nv) = 0._fp
       compute_duv : do ivert = 1,mesh%nv
          if ( mesh%typ(ivert) /= 2 ) cycle
          do ivar = 1,2 ! <-------------------------------+
             call evald1( &                               !
                  tng(:,ivar), &                          !
                  brep%faces(mesh%ids(ivert))%surface, &  !
                  mesh%uv(:,1,ivert), &                   !
                  ivar )                                  !
          end do ! <--------------------------------------+
          call solve_NxN( &
               duv(1:2,ivert), &
               matmul(matmul(transpose(tng), hess(1:3,1:3,ivert)), tng), &
               -matmul(transpose(tng), grad(1:3,ivert)), &
               singular )
          dxyz = matmul(tng, duv(1:2,ivert))
          maxdxyz = max(maxdxyz, sum(dxyz**2))
          maxdxyz_rc = max(maxdxyz_rc, sum(dxyz**2)/rc(ivert))
       end do compute_duv
       PRINT *,'MAX(DXYZ) =',SQRT(MAXDXYZ)
       PRINT *,'MAX(DXYZ/RC) =',SQRT(MAXDXYZ_RC)

       ! update vertices coordinates
       mesh%uv(1:2,1,1:mesh%nv) = mesh%uv(1:2,1,1:mesh%nv) + duv(1:2,1:mesh%nv)
       update_uvxyz : do ivert = 1,mesh%nv
          if ( mesh%typ(ivert) /= 2 ) cycle
          iface = mesh%ids(ivert)
          call eval( &
               mesh%xyz(:,ivert), &
               brep%faces(iface)%surface, &
               mesh%uv(:,1,ivert) )
       end do update_uvxyz

       if ( maxdxyz/hminsqr < EPSdxyzsqr ) then
          PRINT *,'MAX(DXYZ) << MIN(H)'
          exit
       end if

       if ( ipass > 1 ) then
          detot = etotprev - etot
          if ( ipass > 2 ) then
             if ( detot0 < 0._fp ) detot0 = detot
             PRINT *,'DETOT/DETOT0 =',real(detot/detot0)
             if ( detot < 0.01_fp * detot0 ) then ! .or. abs(detot) < 0.005 ) then
                IF ( DETOT < 0._FP ) THEN
                   PRINT *,'NON DECROISSANT'
                ELSE
                   PRINT *,'STATIONNAIRE'
                END IF
                exit
             end if
          else
             detot0 = detot
             PRINT *,'DETOT0 =',real(detot0)
          end if
       end if
       etotprev = etot
       
    end do
    
  end subroutine optim_jiao_uv









  subroutine swap_edge( &
       mesh, &
       ihedg )
    use mod_mesh
    use mod_halfedge
    implicit none
    type(type_surface_mesh), intent(inout) :: mesh
    integer,                 intent(in)    :: ihedg(2)
    integer                                :: hedgpair(2,2)
    integer, dimension(3,2)                :: ijk, tri
    integer                                :: twinj(2,2)
    integer                                :: ih, jh, iface, jface, ivert

    hedgpair(1:2,1) = ihedg
    hedgpair(1:2,2) = get_twin(mesh, ihedg)

    do ih = 1,2
       iface = hedgpair(2,ih)
       ivert = hedgpair(1,ih)
       ! local vertex index in face
       ijk(1,ih) = ivert
       ijk(2,ih) = 1 + mod(ivert,3)
       ijk(3,ih) = 1 + mod(ivert+1,3)
       do ivert = 1,3
          tri(ivert,ih) = mesh%tri(ijk(ivert,ih),iface)
       end do

       twinj(1:2,ih) = mesh%twin(1:2,ijk(2,ih),iface)
    end do

    do ih = 1,2
       jh = 1 + mod(ih,2)
       iface = hedgpair(2,ih)
       jface = hedgpair(2,jh)
       
       ! change v2h
       ivert = tri(1,ih)
       if ( mesh%v2h(1,ivert) == hedgpair(1,ih) .and. &
            mesh%v2h(2,ivert) == hedgpair(2,ih) ) then
          mesh%v2h(1,ivert) = ijk(2,jh)
          mesh%v2h(2,ivert) = jface
       end if

       ! change tri (f2v)
       mesh%tri(ijk(2,ih),iface) = tri(3,jh)

       ! change twins
       mesh%twin(1,ijk(2,ih),iface) = ijk(2,jh)
       mesh%twin(2,ijk(2,ih),iface) = jface

       mesh%twin(1:2,ijk(1,ih),iface) = twinj(1:2,jh)
       if ( all(twinj(1:2,jh) > 0) ) then
          mesh%twin(1,twinj(1,jh),twinj(2,jh)) = ijk(1,ih)
          mesh%twin(2,twinj(1,jh),twinj(2,jh)) = iface
       end if
    end do
    
  end subroutine swap_edge
  
end module mod_optimmesh
