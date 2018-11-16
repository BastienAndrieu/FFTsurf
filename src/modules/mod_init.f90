module mod_init

contains





  subroutine init_from_surfaces( &
       fileoptions, &
       options, &
       surf, &
       nsurf, &
       interdata, &
       brep, &
       hypergraph, &
       mesh )
    use mod_util
    use mod_options
    use mod_diffgeom
    use mod_types_intersection
    use mod_types_brep
    use mod_intersection
    use mod_brep
    use mod_hypergraph
    use mod_mesh
    use mod_polynomial
    use mod_optimmesh
    use mod_tolerances
    ! mode = 0 : propagation des surfaces
    !        1 : regeneration BREP + hypergraphe
    !        2 : isomorphisme hypergrap
    implicit none
    character(*),                    intent(in)    :: fileoptions
    type(type_options),              intent(inout) :: options
    type(type_surface), allocatable, intent(inout) :: surf(:)
    integer,                         intent(out)   :: nsurf
    type(type_intersection_data),    intent(inout) :: interdata
    type(type_brep),                 intent(inout) :: brep
    type(type_hypergraph),           intent(inout) :: hypergraph
    type(type_surface_mesh),         intent(out)   :: mesh
    character(99)                                  :: dir
    character(3)                                   :: strnum3
    integer                                        :: fid
    logical, allocatable                           :: feat_edge(:), feat_vert(:)
    integer                                        :: isurf, icurv, i, j

    ! Read options file
    call read_options( &
         fileoptions, &
         options )
    call print_options( &
         options)

    dir = trim(options%directory)
    if ( options%reprise ) then
       dir = trim(dir) // 'checkpoint/'
    else
       dir = trim(dir) // 'init/'
    end if
    print *,'dir =',trim(dir)

    call get_free_unit(fid)

    ! Import surfaces
    open( &
         unit = fid, &
         file = trim(dir) // 'surftag.dat', &
         action = 'read')
    read (fid,*) nsurf
    allocate(surf(nsurf))
    do isurf = 1,nsurf
       read (fid,*) surf(isurf)%tag
    end do
    close(fid)

    do isurf = 1,nsurf
       write (strnum3,'(i3.3)') isurf
       call read_polynomial( &
            surf(isurf)%x, &
            trim(dir) // 'coef/c_' // strnum3 // '.cheb', &
            nvar=2, &
            base=1 )

       call economize2( &
            surf(isurf)%x, &
            EPSmath )

       call compute_deriv1(surf(isurf))
       call compute_deriv2(surf(isurf))
       call compute_pseudonormal(surf(isurf))
       call economize2( &
            surf(isurf)%pn, &
            EPSmath )
    end do

    if ( options%mode > 0 ) then
       ! Read initial tangential intersection curves
       call read_intersection_curves( &
            trim(dir) // 'tangent_curves.dat', &
            surf, &
            interdata )
       do icurv = 1,interdata%nc
          interdata%curves(icurv)%smooth = .true.
       end do

       ! Generate initial BREP model
       call make_brep( &
            surf, &
            nsurf, &
            interdata, &
            brep, &
            options%chord_err, &
            options%hmin, &
            options%hmax )

       ! debugging >> ..................
       call write_brep_files( &
            brep, &
            trim(dir) // 'brep/verts.dat', &
            trim(dir) // 'brep/edges.dat', &
            trim(dir) // 'brep/faces.dat' )

       call write_intersection_data( &
            interdata, &
            trim(dir) // 'brep/intersection_points.dat', &
            trim(dir) // 'brep/intersection_curves.dat' )
       ! .................. <<

       ! Make hypergraph
       call make_hypergraph( &
            brep, &
            hypergraph, &
            feat_edge, &
            feat_vert )

       ! debugging >> ..................
       call get_free_unit( fid )
       open(unit=fid, file=trim(dir) // 'brep/hyperfaces.dat', action='write')
       write (fid,*) hypergraph%nhf
       do i = 1,hypergraph%nhf
          write (fid,*) hypergraph%hyperfaces(i)%nf
          write (fid,*) hypergraph%hyperfaces(i)%faces(1:hypergraph%hyperfaces(i)%nf)
       end do
       close(fid)

       open(unit=fid, file=trim(dir) // 'brep/hyperedges.dat', action='write')
       write (fid,*) hypergraph%nhe
       do i = 1,hypergraph%nhe
          write (fid,*) hypergraph%hyperedges(i)%ne
          write (fid,*) hypergraph%hyperedges(i)%verts
          write (fid,*) hypergraph%hyperedges(i)%hyperfaces
          do j = 1,hypergraph%hyperedges(i)%ne
             write (fid,*) hypergraph%hyperedges(i)%halfedges(1:2,j)
          end do
       end do
       close(fid)
       ! .................. <<

    end if

    if ( options%mode > 1 ) then      
       ! Generate a first mesh that conforms to the BREP
       call generate_brep_conforming_mesh( &
            brep, &
            options%hmin, &
            options%hmax, &
            options%chord_err, &
            feat_edge(1:brep%ne), &
            feat_vert(1:brep%nv), &
            hypergraph%hyperedges(1:hypergraph%nhe), &
            hypergraph%nhe, &
            mesh )

       if ( allocated(feat_edge) ) deallocate(feat_edge)
       if ( allocated(feat_vert) ) deallocate(feat_vert)
       
       ! debugging >> ..................
       call write_connectivity( &
            '../debug/', &
            mesh, &
            0 )
       call write_xyz_positions( &
            '../debug/', &
            mesh, &
            0 )
       ! .................. <<

       ! Eliminate edges that are too short
       call contract_small_edges( &
            brep, &
            hypergraph%hyperedges(1:hypergraph%nhe), &
            hypergraph%nhe, &
            mesh, &
            0.25_fp*options%hmin )

       ! debugging >> ..................
       call write_connectivity( &
            '../debug/', &
            mesh, &
            1 )
       call write_xyz_positions( &
            '../debug/', &
            mesh, &
            1 )
       ! .................. <<

       ! Make mesh halfedges
       call make_halfedges( &
            mesh )
       !END IF

       ! debugging >> ..................
       call write_mesh_files( &
            mesh, &
            trim(dir) // 'mesh/tri.dat', &
            trim(dir) // 'mesh/xyz.dat', &
            trim(dir) // 'mesh/uv.dat', &
            trim(dir) // 'mesh/idstyp.dat', &
            trim(dir) // 'mesh/paths.dat' )

       call get_free_unit(fid)
       open(unit=fid, file=trim(dir) // 'mesh/mv2h.dat', action='write')
       do i = 1,mesh%nv
          write (fid,*) mesh%v2h(:,i)
       end do
       close(fid)
       open(unit=fid, file=trim(dir) // 'mesh/mtwin.dat', action='write')
       do i = 1,mesh%nt
          write (fid,*) mesh%twin(:,:,i)
       end do
       close(fid)
       ! .................. <<

       ! CHECK UVs
       IF ( .true. ) call check_uvs(brep, mesh)

       PAUSE

       IF ( .TRUE. ) THEN
          do i = 1,interdata%nc
             if ( interdata%curves(i)%smooth ) then
                PRINT *,'REFINE POLYLINE #',I
                call refine_intersection_polyline( &
                     interdata%curves(i), &
                     1.d-4, &
                     options%chord_err )
             end if
          end do
       END IF

       ! Mesh smoothing
       call optim_jiao( &
            brep, &
            hypergraph%hyperedges(1:hypergraph%nhe), &
            hypergraph%nhe, &
            mesh, &
            1._fp, &!PARAM_frac_conf1, &
            0.7_fp, &!PARAM_frac_conf2, &
            5, &!PARAM_ipass1, &
            15, &!PARAM_ipass2, &
            30, &
            options%hmin, &
            options%hmax )


       ! CHECK UVs
       IF ( .true. ) call check_uvs(brep, mesh)

    end if

  end subroutine init_from_surfaces










  subroutine generate_brep_conforming_mesh( &
       brep, &
       hmin, &
       hmax, &
       tol, &
       feat_edge, &
       feat_vert, &
       hyperedges, &
       nhe, &
       mesh )
    use mod_brep
    use mod_hypergraph
    use mod_mesh
    use mod_util
    use mod_polynomial
    use mod_halfedge
    implicit none
    type(type_brep),         intent(in)    :: brep
    real(kind=fp),           intent(in)    :: hmin
    real(kind=fp),           intent(in)    :: hmax
    real(kind=fp),           intent(in)    :: tol
    logical,                 intent(in)    :: feat_edge(brep%ne)
    logical,                 intent(in)    :: feat_vert(brep%nv)
    integer,                 intent(in)    :: nhe
    type(type_hyperedge),    intent(inout) :: hyperedges(nhe)
    type(type_surface_mesh), intent(inout) :: mesh
    integer                                :: finf, fpts, fedg
    integer                                :: bv2mv(brep%nv), be2mv(2,brep%ne)
    integer, allocatable                   :: loc2glob(:), idsf(:), typf(:)
    real(kind=fp), allocatable             :: uvf(:,:,:), xyzf(:,:), uvi(:,:)
    integer, allocatable                   :: trif(:,:)
    integer                                :: npf, np, ntrif
    integer                                :: ifirstedge, ifirstpoint, sens, ifirst, ilast
    type(type_path)                        :: pathtmp
    integer                                :: iface, iwire, ihedg(2), ivert, ihype, iedge, i
    character(3)                           :: strnum

    ! write 'info' file (common to all faces)
    call get_free_unit(finf)
    open(unit=finf, file='../tmp/info.dat', action='write')
    write (finf,*) hmin
    write (finf,*) hmax
    write (finf,*) tol
    close(finf)

    allocate(loc2glob(100), idsf(100), typf(100))
    mesh%nv = 0
    mesh%nt = 0
    bv2mv(:) = 0
    be2mv(:,:) = 0

    faces : do iface = 1,brep%nf ! <-----------------------------------------------------------+
       ! write polynomial parametric surface                                                   !
       call write_polynomial(brep%faces(iface)%surface%x, '../tmp/c.cheb')                     !
       !                                                                                       !
       ! open boundary points & edges files                                                    !
       call get_free_unit(fpts)                                                                !
       open(unit=fpts, file='../tmp/bpts.dat', action='write')                                 !
       call get_free_unit(fedg)                                                                !
       open(unit=fedg, file='../tmp/bedg.dat', action='write')                                 !
       !                                                                                       !
       npf = 0                                                                                 !
       wires : do iwire = 0,brep%faces(iface)%ninner ! <-----------------------------------+   !
          ! first halfedge of the wire                                                     !   !
          if ( iwire == 0 ) then ! <-------------------+                                   !   !
             ihedg = brep%faces(iface)%outer           !                                   !   !
          else ! --------------------------------------+                                   !   !
             ihedg = brep%faces(iface)%inner(:,iwire)  !                                   !   !
          end if ! <-----------------------------------+                                   !   !
          ifirstedge = ihedg(1)                                                            !   !
          !                                                                                !   !
          ifirstpoint = npf + 1                                                            !   !
          ! traverse the wire                                                              !   !
          halfedges : do ! <------------------------------------------------------------+  !   !
             ! get polyline points                                                      !  !   !
             call get_polyline_endpoints( &                                             !  !   !
                  brep, &                                                               !  !   !
                  ihedg, &                                                              !  !   !
                  ifirst, &                                                             !  !   !
                  ilast, &                                                              !  !   !
                  sens, &                                                               !  !   !
                  np )                                                                  !  !   !
             !                                                                          !  !   !
             if ( .not.allocated(loc2glob) .or. &                                       !  !   !
                  size(loc2glob) < npf + np - 1 ) then ! <------+                       !  !   !
                call reallocate_list(loc2glob, npf + np + 100)  !                       !  !   !
             end if ! <-----------------------------------------+                       !  !   !
             !                                                                          !  !   !
             if ( .not.allocated(idsf) .or. &                                           !  !   !
                  size(idsf) < npf + np - 1 ) then ! <----------+                       !  !   !
                call reallocate_list(idsf, npf + np + 100)      !                       !  !   !
             end if ! <-----------------------------------------+                       !  !   !
             !                                                                          !  !   !
             if ( .not.allocated(typf) .or. &                                           !  !   !
                  size(typf) < npf + np - 1 ) then ! <----------+                       !  !   !
                call reallocate_list(typf, npf + np + 100)      !                       !  !   !
             end if ! <-----------------------------------------+                       !  !   !
             !                                                                          !  !   !
             if ( allocated(uvf) .and. size(uvf,3) < np ) deallocate(uvf)               !  !   !
             if ( .not.allocated(uvf) ) allocate(uvf(2,2,np))                           !  !   !
             ! polyligne uv dans le sens du contour de la face                          !  !   !
             if ( sens == 2 ) then ! <----------------------------------------+         !  !   !
                uvf(1:2,1:2,1:np) = brep%edges(ihedg(1))%curve%polyline%uv&   !         !  !   !
                     (1:2,[2,1],ifirst:ilast)                                 !         !  !   !
             else ! ----------------------------------------------------------+         !  !   !
                uvf(1:2,1:2,1:np) = brep%edges(ihedg(1))%curve%polyline%uv&   !         !  !   !
                     (1:2,1:2,ifirst:ilast:-1)                                !         !  !   !
             end if ! <-------------------------------------------------------+         !  !   !
             !                                                                          !  !   !
             if ( allocated(xyzf) .and. size(xyzf,2) < np ) deallocate(xyzf)            !  !   !
             if ( .not.allocated(xyzf) ) allocate(xyzf(3,np))                           !  !   !
             xyzf(1:3,1:np) = brep%edges(ihedg(1))%curve%polyline%xyz&                  !  !   !
                  (1:3,ifirst:ilast:(-1)**sens)                                         !  !   !
             !                                                                          !  !   !
             ! get brep vertex index of first point                                     !  !   !
             ivert = get_orig(brep, ihedg)                                              !  !   !
             if ( bv2mv(ivert) == 0 ) then ! <-----------------------------------+      !  !   !
                ! the brep vertex is not yet present in the mesh                 !      !  !   !
                bv2mv(ivert) = mesh%nv + 1                                       !      !  !   !
                if ( feat_vert(ivert) ) then ! <------------------------------+  !      !  !   !
                   ! the mesh vertex is locked on a brep vertex               !  !      !  !   !
                   idsf(npf+1) = ivert                                        !  !      !  !   !
                   typf(npf+1) = 0                                            !  !      !  !   !
                else ! -------------------------------------------------------+  !      !  !   !
                   if ( feat_edge(ihedg(1)) ) then ! <---------------------+  !  !      !  !   !
                      ! the mesh vertex is locked on a hyperedge           !  !  !      !  !   !
                      idsf(npf+1) = ihedg(1)                               !  !  !      !  !   !
                      typf(npf+1) = 1                                      !  !  !      !  !   !
                   else ! -------------------------------------------------+  !  !      !  !   !
                      ! the mesh vertex is free to move along a hyperface  !  !  !      !  !   !
                      idsf(npf+1) = iface                                  !  !  !      !  !   !
                      typf(npf+1) = 2                                      !  !  !      !  !   !
                   end if ! <----------------------------------------------+  !  !      !  !   !
                end if ! <----------------------------------------------------+  !      !  !   !
                !                                                                !      !  !   !
                ! append new vertex to the mesh                                  !      !  !   !
                call append_vertices( &                                          !      !  !   !
                     mesh, &                                                     !      !  !   !
                     xyzf(1:3,1), &                                              !      !  !   !
                     uvf(1:2,1:2,1), &                                           !      !  !   !
                     idsf(npf+1), &                                              !      !  !   !
                     typf(npf+1), &                                              !      !  !   !
                     1 )                                                         !      !  !   !
             else ! -------------------------------------------------------------+      !  !   !
                if ( feat_edge(ihedg(1)) ) then ! <------------------------+     !      !  !   !
                   if ( mesh%typ(bv2mv(ivert)) == 2 ) then ! <----------+  !     !      !  !   !
                      mesh%ids(bv2mv(ivert)) = ihedg(1)                 !  !     !      !  !   !
                      mesh%typ(bv2mv(ivert)) = 1                        !  !     !      !  !   !
                      mesh%uv(1:2,1:2,bv2mv(ivert)) = uvf(1:2,1:2,1)    !  !     !      !  !   !
                   end if ! <-------------------------------------------+  !     !      !  !   !
                end if ! <-------------------------------------------------+     !      !  !   !
             end if ! <----------------------------------------------------------+      !  !   !
             loc2glob(npf+1) = bv2mv(ivert)                                             !  !   !
             !                                                                          !  !   !
             if ( be2mv(1,ihedg(1)) == 0 ) then ! <------------------------------+      !  !   !
                ! the brep edge is not yet present in the mesh                   !      !  !   !
                be2mv(1,ihedg(1)) = mesh%nv + 1                                  !      !  !   !
                be2mv(2,ihedg(1)) = mesh%nv + np - 2                             !      !  !   !
                !                                                                !      !  !   !
                if ( feat_edge(ihedg(1)) ) then ! <-----------------------+      !      !  !   !
                   ! the mesh vertices are locked on a hyperedge          !      !      !  !   !
                   idsf(npf+2:npf+np-1) = ihedg(1)                        !      !      !  !   !
                   typf(npf+2:npf+np-1) = 1                               !      !      !  !   !
                else ! ---------------------------------------------------+      !      !  !   !
                   ! the mesh vertices are free to move along a hyperface !      !      !  !   !
                   idsf(npf+2:npf+np-1) = iface                           !      !      !  !   !
                   typf(npf+2:npf+np-1) = 2                               !      !      !  !   !
                end if ! <------------------------------------------------+      !      !  !   !
                !                                                                !      !  !   !
                if ( sens == 1 ) then ! <-----------------------------------+    !      !  !   !
                   loc2glob(npf+2:npf+np-1) = mesh%nv + [(i, i=np-2,1,-1)]  !    !      !  !   !
                   call append_vertices( &                                  !    !      !  !   !
                        mesh, &                                             !    !      !  !   !
                        xyzf(1:3,np-1:2:-1), &                              !    !      !  !   !
                        uvf(1:2,1:2,np-1:2:-1), &                           !    !      !  !   !
                        idsf(npf+2:npf+np-1), &                             !    !      !  !   !
                        typf(npf+2:npf+np-1), &                             !    !      !  !   !
                        np-2 )                                              !    !      !  !   !
                else ! -----------------------------------------------------+    !      !  !   !
                   loc2glob(npf+2:npf+np-1) = mesh%nv + [(i, i=1,np-2)]     !    !      !  !   !
                   call append_vertices( &                                  !    !      !  !   !
                        mesh, &                                             !    !      !  !   !
                        xyzf(1:3,2:np-1), &                                 !    !      !  !   !
                        uvf(1:2,1:2,2:np-1), &                              !    !      !  !   !
                        idsf(npf+2:npf+np-1), &                             !    !      !  !   !
                        typf(npf+2:npf+np-1), &                             !    !      !  !   !
                        np-2 )                                              !    !      !  !   !
                end if ! <--------------------------------------------------+    !      !  !   !
             else ! -------------------------------------------------------------+      !  !   !
                if ( sens == 1 ) then ! <-----------------------------------+    !      !  !   !
                   loc2glob(npf+2:npf+np-1) = &                             !    !      !  !   !
                        [(i, i=be2mv(2,ihedg(1)),be2mv(1,ihedg(1)),-1)]     !    !      !  !   !
                else ! -----------------------------------------------------+    !      !  !   !
                   loc2glob(npf+2:npf+np-1) = &                             !    !      !  !   !
                        [(i, i=be2mv(1,ihedg(1)),be2mv(2,ihedg(1)))]        !    !      !  !   !
                end if ! <--------------------------------------------------+    !      !  !   !
             end if ! <----------------------------------------------------------+      !  !   !
             !                                                                          !  !   !
             ! write points...                                                          !  !   !
             do i = 1,np-1 ! <---------------------------+                              !  !   !
                write (fpts,*) uvf(:,1,i)                !                              !  !   !
             end do ! <----------------------------------+                              !  !   !
             ! ...and edges                                                             !  !   !
             do i = 1,np-2 ! <---------------------------+                              !  !   !
                write (fedg,*) npf + i, npf + i + 1      !                              !  !   !
             end do ! <----------------------------------+                              !  !   !
             !                                                                          !  !   !
             npf = npf + np - 1                                                         !  !   !
             !                                                                          !  !   !
             ! move on to the next halfedge on the wire                                 !  !   !
             ihedg = get_next(brep, ihedg)                                              !  !   !
             !                                                                          !  !   !
             if ( ihedg(1) == ifirstedge ) then ! <------+                              !  !   !
                ! the wire is complete                   !                              !  !   !
                write (fedg,*) npf, ifirstpoint          !                              !  !   !
                exit                                     !                              !  !   !
             else ! -------------------------------------+                              !  !   !
                write (fedg,*) npf, npf+1                !                              !  !   !
             end if ! <----------------------------------+                              !  !   !
             !                                                                          !  !   !
          end do halfedges ! <----------------------------------------------------------+  !   !
       end do wires ! <--------------------------------------------------------------------+   !
       !                                                                                       !
       ! close files                                                                           !
       close(fpts)                                                                             !
       close(fedg)                                                                             !
       !                                                                                       !
       ! run mesher                                                                            !
       PRINT *,'MESH FACE #',IFACE
       call system( &!'/home/bastien/MeshGen/./meshgen.out &
            '/stck/bandrieu/Bureau/MeshGen/./meshgen.out &
            & ../tmp/c.cheb &
            & ../tmp/bpts.dat &
            & ../tmp/bedg.dat &
            & ../tmp/info.dat &
            & ../tmp/tri.dat &
            & ../tmp/uv.dat &
            & ../tmp/xyz.dat' )
       write (strnum,'(i3.3)') iface                                                           !
       ! read face submesh                                                                     !
       call read_triangles( &                                                                  !
            '../tmp/tri.dat', &                                                                !
            trif, &                                                                            !
            ntrif )                                                                            !
       !                                                                                       !
       call read_points( &                                                                     !
            '../tmp/uv.dat', &                                                                 !
            uvi, &                                                                             !
            2, &                                                                               !
            np )                                                                               !
       if ( allocated(uvf) ) deallocate(uvf)                                                   !
       allocate(uvf(2,2,np-npf))                                                               !
       uvf(1:2,1,1:np-npf) = uvi(1:2,npf+1:np)                                                 !
       !                                                                                       !
       call read_points( &                                                                     !
            '../tmp/xyz.dat', &                                                                !
            xyzf, &                                                                            !
            3, &                                                                               !
            np )                                                                               !
       !                                                                                       !
       if ( size(loc2glob) < np ) call reallocate_list(loc2glob, np)                           !
       loc2glob(npf+1:np) = mesh%nv + [(i, i=1,np-npf)]                                        !
       !                                                                                       !
       if ( size(idsf) < np ) call reallocate_list(idsf, np)                                   !
       idsf(npf+1:np) = iface                                                                  !
       !                                                                                       !
       if ( size(typf) < np ) call reallocate_list(typf, np)                                   !
       typf(npf+1:np) = 2                                                                      !
       !                                                                                       !
       do i = 1,3 ! <----------------------------------+                                       !
          trif(i,1:ntrif) = loc2glob(trif(i,1:ntrif))  !                                       !
       end do ! <--------------------------------------+                                       !
       ! add new triangles                                                                     !
       call append_triangles( &                                                                !
            mesh, &                                                                            !
            trif(1:3,1:ntrif), &                                                               !
            [(brep%faces(iface)%hyperface, i=1,ntrif)], &                                      !
            ntrif )                                                                            !
       !                                                                                       !
       ! add new mesh vertices                                                                 !
       call append_vertices( &                                                                 !
            mesh, &                                                                            !
            xyzf(1:3,npf+1:np), &                                                              !
            uvf(1:2,1:2,1:np-npf), &                                                           !
            idsf(npf+1:np), &                                                                  !
            typf(npf+1:np), &                                                                  !
            np-npf )                                                                           !
       !                                                                                       !
    end do faces ! <---------------------------------------------------------------------------+

    ! create a path for each hyperedge
    do ihype = 1,nhe
       pathtmp%nv = 0
       pathtmp%hyperedge = ihype
       !
       if ( hyperedges(ihype)%verts(1) < 1 ) then ! <---+
          ihedg = hyperedges(ihype)%halfedges(:,1)      !
          ivert = get_orig(brep, ihedg)                 !
          hyperedges(ihype)%verts(1:2) = ivert          !
       end if ! <---------------------------------------+
       !
       do iedge = 1,hyperedges(ihype)%ne ! <---------------+
          ihedg = hyperedges(ihype)%halfedges(:,iedge)     !
          sens = 1 + mod(ihedg(2),2)                       !
          !                                                !
          call insert_after( &                             !
               pathtmp%verts, &                            !
               pathtmp%nv, &                               !
               bv2mv(get_orig(brep, ihedg)), &             !
               pathtmp%nv )                                !
          !                                                !
          np = be2mv(2,ihedg(1)) - be2mv(1,ihedg(1)) + 1   !
          if ( sens == 1 ) then ! <------+                 !
             ifirst = be2mv(2,ihedg(1))  !                 !
             ilast = be2mv(1,ihedg(1))   !                 !
          else ! ------------------------+                 !
             ifirst = be2mv(1,ihedg(1))  !                 !
             ilast = be2mv(2,ihedg(1))   !                 !
          end if ! <---------------------+                 !
          !                                                !
          call insert_n_after( &                           !
               pathtmp%verts, &                            !
               pathtmp%nv, &                               !
               np, &                                       !
               [(i, i=ifirst,ilast,(-1)**sens)], &         !
               pathtmp%nv )                                !
       end do ! <------------------------------------------+
       !
       if ( hyperedges(ihype)%verts(2) > 0 ) then ! <---+
          call insert_after( &                          !
               pathtmp%verts, &                         !
               pathtmp%nv, &                            !
               bv2mv(hyperedges(ihype)%verts(2)), &     !
               pathtmp%nv )                             !
       end if ! <---------------------------------------+
       !
       call append_path( &
            mesh, &
            pathtmp )
       !
       ! compute discrete curvilinear abscissa along the path
       allocate(mesh%paths(mesh%npaths)%s(mesh%paths(mesh%npaths)%nv))
       mesh%paths(mesh%npaths)%s(1) = 0._fp
       do i = 2,mesh%paths(mesh%npaths)%nv ! <-------------------------------+
          mesh%paths(mesh%npaths)%s(i) = mesh%paths(mesh%npaths)%s(i-1) + &  !
               norm2( &                                                      !
               mesh%xyz(:,mesh%paths(mesh%npaths)%verts(i)) - &              !
               mesh%xyz(:,mesh%paths(mesh%npaths)%verts(i-1)) )              !
       end do ! <------------------------------------------------------------+
       !
    end do

  end subroutine generate_brep_conforming_mesh








  subroutine read_msh( &
       filename, &
       xyz, &
       nv, &
       tri, &
       ref, &
       nt )
    use mod_math
    use mod_util
    implicit none
    character(*),               intent(in)    :: filename
    real(kind=fp), allocatable, intent(inout) :: xyz(:,:)
    integer,                    intent(out)   :: nv
    integer, allocatable,       intent(inout) :: tri(:,:)
    integer, allocatable,       intent(inout) :: ref(:)
    integer,                    intent(out)   :: nt
    integer                                   :: fid, err
    character(100)                            :: line
    real(kind=fp)                             :: vert(4)
    integer                                   :: ne, ntags, elem(10)
    logical, allocatable                      :: keepv(:)
    integer, allocatable                      :: old2new(:)
    integer                                   :: ivert, ielem

    call get_free_unit(fid)
    open( &
         unit = fid, &
         file = filename, &
         action = 'read' )
    do ! <-------------------------------------------------------+
       read (fid, *, iostat=err) line                            !
       if ( err /= 0 ) exit                                      !
       line = trim(line)                                         !
       !                                                         !
       if ( line(1:6) == '$Nodes' ) then ! <-----------------+   !
          read (fid,*) nv                                    !   !
          if ( allocated(xyz) ) then ! <---------------+     !   !
             if ( size(xyz,2) < nv ) deallocate(xyz)   !     !   !
          end if ! <-----------------------------------+     !   !
          if ( .not.allocated(xyz) ) allocate(xyz(3,nv))     !   !
          !                                                  !   !
          do ivert = 1,nv ! <-------------+                  !   !
             read (fid,*) vert            !                  !   !
             xyz(1:3,ivert) = vert(2:4)   !                  !   !
          end do ! <----------------------+                  !   !
          !                                                  !   !
       elseif ( line(1:9) == '$Elements' ) then ! -----------+   !
          read (fid,*) ne                                    !   !
          if ( allocated(tri) ) then ! <-----------------+   !   !
             if ( size(tri,2) < ne ) deallocate(tri)     !   !   !
          end if ! <-------------------------------------+   !   !
          if ( .not.allocated(tri) ) allocate(tri(3,ne))     !   !
          if ( allocated(ref) ) then ! <----------------+   !   !
             if ( size(ref) < ne ) deallocate(ref)       !   !   !
          end if ! <-------------------------------------+   !   !
          if ( .not.allocated(ref) ) allocate(ref(ne))       !   !
          !                                                  !   !
          nt = 0                                             !   !
          do ielem = 1,ne ! <----------------------------+   !   !
             read (fid,*) elem(1:3)                      !   !   !
             if ( elem(2) == 2 ) then ! <------------+   !   !   !
                nt = nt + 1                          !   !   !   !
                ntags = elem(3)                      !   !   !   !
                backspace(fid)                       !   !   !   !
                read (fid,*) elem(1:ntags+6)         !   !   !   !
                tri(1:3,nt) = elem(ntags+4:ntags+6)  !   !   !   !
                ref(nt) = elem(ntags+3)              !   !   !   !
             end if ! <------------------------------+   !   !   !
          end do ! <-------------------------------------+   !   !
          !                                                  !   !
       end if ! <--------------------------------------------+   !
    end do ! <---------------------------------------------------+
    close(fid)


    ! fix
    allocate(keepv(nv))
    keepv(1:nv) = .false.
    do ielem = 1,nt
       keepv(tri(1:3,ielem)) = .true.
    end do

    allocate(old2new(nv))
    old2new(1:nv) = 0
    nv = 0
    do ivert = 1,size(old2new)
       if ( keepv(ivert) ) then
          nv = nv + 1
          old2new(ivert) = nv
       end if
    end do

    do ielem = 1,nt
       tri(1:3,ielem) = old2new(tri(1:3,ielem))
    end do

    xyz(1:3,1:nv) = reshape( &
         pack(xyz(1:3,1:size(keepv)), spread(keepv, dim=1, ncopies=3)), &
         [3,nv] )

    deallocate(keepv, old2new)

  end subroutine read_msh










  subroutine check_uvs( &
       brep, &
       mesh )
    use mod_types_brep
    use mod_mesh
    use mod_diffgeom
    use mod_tolerances
    implicit none
    type(type_brep),         intent(in) :: brep
    type(type_surface_mesh), intent(in) :: mesh
    real(kind=fp)                       :: xyzverif(3,2)
    integer                             :: i, j, iface

    PRINT *,'CHECK UVs'
    do i = 1,mesh%nv
       select case(mesh%typ(i))
       case (0)
       case (1)
          do j = 1,2
             iface = brep%edges(mesh%ids(i))%halfedges(1+mod(j,2))%face
             call eval( &
                  xyzverif(:,j), &
                  brep%faces(iface)%surface, &
                  mesh%uv(:,j,i) )
          end do
          IF ( MAX(norm2(mesh%xyz(:,i) - xyzverif(:,1)), norm2(mesh%xyz(:,i) - xyzverif(:,2))) > EPSxyz ) THEN
             PRINT *,'I =', I,', ERR =', norm2(mesh%xyz(:,i) - xyzverif(:,1)), norm2(mesh%xyz(:,i) - xyzverif(:,2))
             do j = 1,2
                iface = brep%edges(mesh%ids(i))%halfedges(j)%face
                call eval( &
                     xyzverif(:,j), &
                     brep%faces(iface)%surface, &
                     mesh%uv(:,j,i) )
             end do
             PRINT *,'    I =', I,', ERR =', norm2(mesh%xyz(:,i) - xyzverif(:,1)), norm2(mesh%xyz(:,i) - xyzverif(:,2))
          END IF
       case (2)
          call eval( &
               xyzverif(:,1), &
               brep%faces(mesh%ids(i))%surface, &
               mesh%uv(:,1,i) )
          IF ( norm2(mesh%xyz(:,i) - xyzverif(:,1)) > 1.D-13 ) THEN
             PRINT *,'I =', I,', ERR =', norm2(mesh%xyz(:,i) - xyzverif(:,1))
          END IF
       end select
    end do

  end subroutine check_uvs























  subroutine init_from_mesh( &
       fileoptions, &
       options, &
       surf, &
       nsurf, &
       interdata, &
       brep, &
       hypergraph, &
       mesh )
    use mod_import
    use mod_util
    use mod_options
    use mod_diffgeom
    use mod_types_intersection
    use mod_types_brep
    use mod_intersection
    use mod_brep
    use mod_hypergraph
    use mod_mesh
    use mod_polynomial
    use mod_optimmesh
    use mod_tolerances
    use mod_geometry
    use mod_projection
    implicit none
    real(kind=fp), parameter                               :: TOLang = 30._fp ! [degrees]
    real(kind=fp), parameter                               :: TOLcosang = cos(CSTpi*TOLang/180._fp)
    character(*),                            intent(in)    :: fileoptions
    type(type_options),                      intent(inout) :: options
    type(type_surface), allocatable, target, intent(inout) :: surf(:)
    integer,                                 intent(out)   :: nsurf
    type(type_intersection_data), target,    intent(inout) :: interdata
    type(type_brep),                         intent(inout) :: brep
    type(type_hypergraph),                   intent(inout) :: hypergraph
    type(type_surface_mesh),                 intent(out)   :: mesh
    character(100)                                         :: dir, filemesh
    character(3)                                           :: strnum3
    integer                                                :: fid
    integer, allocatable                                   :: remap(:)
    integer                                                :: bv2bf(3), nbv2bf
    real(kind=fp)                                          :: uvpoint(2,0)
    type(ptr_surface)                                      :: surfpoint(0)
    integer                                                :: idpoint
    logical, dimension(:,:), allocatable                   :: visited, bnd_hedg
    real(kind=fp)                                          :: nor(3,2)
    integer                                                :: surfpair(2)
    integer, allocatable                                   :: polyverts(:), curve2verts(:,:)
    integer                                                :: npolyverts, nc
    type(type_intersection_curve), pointer                 :: curve => null()
    type(type_point_on_surface), pointer                   :: pos => null()
    integer                                                :: statproj
    real(kind=fp)                                          :: xyzproj(3)
    logical, dimension(:), allocatable                     :: feat_edge, feat_vert
    type(type_path)                                        :: path
    integer                                                :: head, tail, stride
    integer                                                :: isurf
    integer                                                :: iface, jface
    integer                                                :: iref, jref
    integer, dimension(2)                                  :: ihedg, jhedg
    integer                                                :: iedge
    integer                                                :: iendpoint, ipos, imv, iip, ipv, jpv
    integer                                                :: ivert, jvert, kvert, lvert
    integer                                                :: ihype, icurv
    integer                                                :: nfeatedg
    INTEGER :: I, J
    real(kind=fp), dimension(3,2) :: xyzverif

    ! Read options file
    call read_options( &
         fileoptions, &
         options )
    call print_options( &
         options)
    
    dir = trim(options%directory)
    if ( options%reprise ) then
       dir = trim(dir) // 'checkpoint/'
    else
       dir = trim(dir) // 'init/'
    end if
    
    ! Import surfaces
    call get_free_unit(fid)
    open( &
         unit = fid, &
         !file = 'import_geometry/surftag.dat', & !*****
         file = trim(dir) // 'surftag.dat', &
         action = 'read')
    read (fid,*) nsurf
    allocate(surf(nsurf))
    do isurf = 1,nsurf
       read (fid,*) surf(isurf)%tag
    end do
    close(fid)

    PRINT *,'NSURF =',NSURF
    do isurf = 1,nsurf
       write (strnum3,'(i3.3)') isurf
       call read_polynomial( &
            surf(isurf)%x, &
            !'import_geometry/coef/c_' // strnum3 // '.cheb', & !*****
            trim(dir) // 'coef/c_' // strnum3 // '.cheb', &
            nvar=2, &
            base=1 )

       call economize2( &
            surf(isurf)%x, &
            EPSmath )

       call compute_deriv1(surf(isurf))
       call compute_deriv2(surf(isurf))
       call compute_pseudonormal(surf(isurf))
       call economize2( &
            surf(isurf)%pn, &
            EPSmath )
    end do
    
    
    ! Read xyz, tri and face ref
    !filemesh = 'import_geometry/base_sym2.msh' !*****
    filemesh = trim(dir) // 'mesh/initmesh.msh'
    call read_msh( &
         trim(filemesh), &
         mesh%xyz, &
         mesh%nv, &
         mesh%tri, &
         mesh%ihf, &
         mesh%nt )

    ! fix mesh face ref -> contiguous BREP face indices
    allocate(remap(maxval(mesh%ihf(1:mesh%nt))))
    remap(:) = -1
    brep%nf = 0
    do iface = 1,mesh%nt
       iref = mesh%ihf(iface)
       if ( remap(iref) < 0 ) then
          brep%nf = brep%nf + 1
          remap(iref) = brep%nf
       end if
    end do
    mesh%ihf(1:mesh%nt) = remap(mesh%ihf(1:mesh%nt))
    deallocate(remap)

    allocate(mesh%uv(2,2,mesh%nv), mesh%ids(mesh%nv), mesh%typ(mesh%nv))
    mesh%typ(1:mesh%nv) = 2
    ! read uv coordinates
    call get_free_unit(fid)
    open( &
         unit = fid, &
         !file = 'import_geometry/uv0.dat', & !*****
         file = trim(dir) // 'mesh/uv0.dat', &
         action = 'read' )
    do ivert = 1,mesh%nv
       read (fid,*) mesh%uv(1:2,1:2,ivert)
    end do
    close(fid)
    
    ! read ids
    call get_free_unit(fid)
    open( &
         unit = fid, &
         !file = 'import_geometry/ids0.dat', & !*****
         file = trim(dir) // 'mesh/ids0.dat', &
         action = 'read' )
    do ivert = 1,mesh%nv
       read (fid,*) mesh%ids(ivert)
    end do
    close(fid)

    ! make mesh halfedges
    call make_halfedges( &
         mesh )

    ! make intersection points
    do ivert = 1,mesh%nv ! <------------------------------------+
       nbv2bf = 0                                               !
       ihedg = mesh%v2h(:,ivert) ! outgoing                     !
       iface = get_face(ihedg)                                  !
       jface = iface                                            !
       do ! <-----------------------------------------------+   !
          if ( jface < 1 ) then ! <---+                     !   !
             jref = -1                !                     !   !
          else ! ---------------------+                     !   !
             jref = mesh%ihf(jface)   !                     !   !
          end if ! -------------------+                     !   !
          !                                                 !   !
          do iref = 1,nbv2bf ! <------------------------+   !   !
             if ( jref == bv2bf(iref) ) exit            !   !   !
          end do ! <------------------------------------+   !   !
          !                                                 !   !
          if ( nbv2bf < 1 .or. iref > nbv2bf ) then ! <-+   !   !
             nbv2bf = nbv2bf + 1                        !   !   !
             bv2bf(nbv2bf) = jref                       !   !   !
          end if ! <------------------------------------+   !   !
          if ( nbv2bf > 2 ) exit                            !   !
          !                                                 !   !
          if ( jface < 1 ) exit                             !   !
          !                                                 !   !
          ihedg = get_prev(ihedg)       ! ingoing           !   !
          ihedg = get_twin(mesh, ihedg) ! outgoing          !   !
          jface = get_face(ihedg)                           !   !
          if ( jface == iface ) exit                        !   !
       end do ! <-------------------------------------------+   !
       !                                                        !
       if ( nbv2bf > 2 ) then ! <---------------------------+   !
          call add_intersection_point( &                    !   !
               uvpoint, &                                   !   !
               mesh%xyz(1:3,ivert), &                       !   !
               surfpoint, &                                 !   !
               0, &                                         !   !
               interdata, &                                 !   !
               idpoint )                                    !   !
          interdata%points(idpoint)%ivert = ivert           !   !
          mesh%typ(ivert) = 0                               !   !
          mesh%ids(ivert) = idpoint                         !   !
       end if ! <-------------------------------------------+   !
    end do ! <--------------------------------------------------+

    ! mark halfedges on BREP faces boundaries
    allocate(visited(3,mesh%nt), bnd_hedg(3,mesh%nt))
    visited(1:3,1:mesh%nt) = .false.
    do iface = 1,mesh%nt ! <--------------------------------------+
       iref = mesh%ihf(iface)                                     !
       do iedge = 1,3 ! <--------------------------------------+  !
          if ( visited(iedge,iface) ) cycle                    !  !
          ihedg = mesh%twin(1:2,iedge,iface)                   !  !
          if ( ihedg(2) < 1 ) then ! <--+                      !  !
             jref = -1                  !                      !  !
          else ! -----------------------+                      !  !
             jref = mesh%ihf(ihedg(2))  !                      !  !
          end if ! <--------------------+                      !  !
          bnd_hedg(iedge,iface) = ( jref /= iref )             !  !
          visited(iedge,iface) = .true.                        !  !
          bnd_hedg(ihedg(1),ihedg(2)) = bnd_hedg(iedge,iface)  !  !
          visited(ihedg(1),ihedg(2)) = .true.                  !  !
       end do ! <----------------------------------------------+  !
    end do ! <----------------------------------------------------+
    
    ! make intersection curves
    nc = 0
    visited(1:3,1:mesh%nt) = .false.
    allocate(polyverts(100))
    do jvert = 1,interdata%np
       ivert = interdata%points(jvert)%ivert
       ihedg = mesh%v2h(:,ivert)
       iface = get_face(ihedg)
       jface = iface
       do
          if ( .not.visited(ihedg(1),ihedg(2)) .and. &
               bnd_hedg(ihedg(1),ihedg(2)) ) then ! <-----------------------------------------------+
             ! add a new intersection curve                                                         !
             call extract_intersection_polyline_vertices( &                                         !
                  mesh, &                                                                           !
                  ihedg, &                                                                          !
                  visited, &                                                                        !
                  bnd_hedg, &                                                                       !
                  polyverts, &                                                                      !
                  npolyverts )                                                                      !
             !                                                                                      !
             call add_intersection_curve( &                                                         !
                  interdata, &                                                                      !
                  [0._fp, 0._fp, 0._fp], &                                                          !
                  mesh%ids(polyverts([1,npolyverts])), &                                            !
                  spread(spread([-1._fp, 1._fp], 2, 2), 3, 2) )                                     !
             ! is it a smooth (tangential) intersection?                                            !
             curve => interdata%curves(interdata%nc)                                                !
             jhedg = get_twin(mesh, ihedg)                                                          !
             jface = get_face(jhedg)                                                                !
             if ( jface < 1 ) then ! <------------------------------------------+                   !
                curve%smooth = .false.                                          !                   !
                surfpair(1) = 0                                                 !                   !
             else ! ------------------------------------------------------------+                   !
                nor(:,1) = triangle_normal( &                                   !                   !
                     mesh%xyz(:,mesh%tri(:,ihedg(2))), &                        !                   !
                     .true. )                                                   !                   !
                nor(:,2) = triangle_normal( &                                   !                   !
                     mesh%xyz(:,mesh%tri(:,jface)), &                           !                   !
                     .true. )                                                   !                   !
                curve%smooth = ( dot_product(nor(:,1), nor(:,2)) > TOLcosang )  !                   !
                surfpair(1) = mesh%ihf(jface)                                   !                   !
             end if ! <---------------------------------------------------------+                   !
             !                                                                                      !
             surfpair(2) = mesh%ihf(ihedg(2))                                                       !
             if ( surfpair(2) < surfpair(1) ) then ! <--------------------------+                   !
                mesh%uv(1:2,1:2,polyverts(2:npolyverts-1)) = &                  !                   !
                     mesh%uv(1:2,[2,1],polyverts(2:npolyverts-1))               !                   !
             end if ! <---------------------------------------------------------+                   !
             do isurf = 1,2 ! <-------------------------------------------------+                   !
                if ( surfpair(isurf) > 0 ) curve%surf(isurf)%ptr => &           !                   !
                     surf(surfpair(isurf))                                      !                   !
             end do ! <---------------------------------------------------------+                   !
             !                                                                                      !
             curve%isplit(2,1:2) = [1,npolyverts]                                                   !
             allocate(curve%iedge(1))                                                               !
             curve%iedge = 0                                                                        !
             allocate(curve%polyline)                                                               !
             curve%polyline%np = npolyverts                                                         !
             allocate(curve%polyline%xyz(3,npolyverts), &                                           !
                  curve%polyline%uv(2,2,npolyverts))                                                !
             curve%polyline%xyz(1:3,1:npolyverts) = &                                               !
                  mesh%xyz(1:3,polyverts(1:npolyverts))                                             !
             curve%polyline%uv(1:2,1:2,2:npolyverts-1) = &                                          !
                  mesh%uv(1:2,1:2,polyverts(2:npolyverts-1))                                        !
             !                                                                                      !
             if ( curve%smooth ) then ! <---------------------------+                               !
                mesh%typ(polyverts(2:npolyverts-1)) = 2             !                               !
                mesh%ids(polyverts(2:npolyverts-1)) = surfpair(1)   !                               !
             else ! ------------------------------------------------+                               !
                mesh%typ(polyverts(2:npolyverts-1)) = 1             !                               !
                mesh%ids(polyverts(2:npolyverts-1)) = interdata%nc  !                               !
             end if ! <---------------------------------------------+                               !
             !                                                                                      !
             call insert_column_after( &                                                            !
                  curve2verts, &                                                                    !
                  npolyverts, &                                                                     !
                  nc, &                                                                             !
                  polyverts(1:npolyverts), &                                                        !
                  nc )                                                                              !
             !                                                                                      !
             do iendpoint = 1,2 ! <--------------------------------------------------------------+  !
                ipv = 1 + (npolyverts - 1)*(iendpoint - 1) ! polyline vertex                     !  !
                jpv = ipv + (-1)**(iendpoint + 1) ! nearest interior poly. vertex                !  !
                imv = polyverts(ipv) ! mesh vertex                                               !  !
                iip = mesh%ids(imv) ! intersection point                                         !  !
                !                                                                                !  !
                incident_surfaces : do isurf = 1,2 ! <----------------------------------------+  !  !
                   if ( interdata%points(iip)%npos > 0 ) then ! <--------------------------+  !  !  !
                      pos => interdata%points(iip)%pos                                     !  !  !  !
                      do ipos = 1,interdata%points(iip)%npos ! <------------------------+  !  !  !  !
                         if ( associated(pos%surf,curve%surf(isurf)%ptr) ) then ! <--+  !  !  !  !  !
                            curve%polyline%uv(1:2,isurf,ipv) = pos%uv                !  !  !  !  !  !
                            cycle incident_surfaces                                  !  !  !  !  !  !
                         end if ! <--------------------------------------------------+  !  !  !  !  !
                         if ( ipos < interdata%points(iip)%npos ) pos => pos%next       !  !  !  !  !
                      end do ! <--------------------------------------------------------+  !  !  !  !
                      interdata%points(iip)%npos = interdata%points(iip)%npos + 1          !  !  !  !
                      allocate(pos%next)                                                   !  !  !  !
                      pos => pos%next                                                      !  !  !  !
                   else ! -----------------------------------------------------------------+  !  !  !
                      interdata%points(iip)%npos = 1                                       !  !  !  !
                      allocate(interdata%points(iip)%pos)                                  !  !  !  !
                      pos => interdata%points(iip)%pos                                     !  !  !  !
                   end if ! <--------------------------------------------------------------+  !  !  !
                   pos%surf => curve%surf(isurf)%ptr                                          !  !  !
                   pos%uv = curve%polyline%uv(1:2,isurf,jpv)                                  !  !  !
                   call projection_surface( &                                                 !  !  !
                        pos%surf, &                                                           !  !  !
                        mesh%xyz(1:3,imv), &                                                  !  !  !
                        pos%uv, &                                                             !  !  !
                        -(1._fp + EPSuv)*[1._fp, 1._fp], &                                    !  !  !
                        (1._fp + EPSuv)*[1._fp, 1._fp], &                                     !  !  !
                        statproj, &                                                           !  !  !
                        xyzproj )                                                             !  !  !
                   IF ( STATPROJ > 0 ) THEN
                      PRINT *,'FAILED TO PROJECT VERTEX #', imv, ' ON SURFACE #', SURFPAIR(ISURF)
                      PAUSE
                   END IF
                   curve%polyline%uv(1:2,isurf,ipv) = pos%uv                                  !  !  !
                end do incident_surfaces ! <--------------------------------------------------+  !  !
             end do ! <--------------------------------------------------------------------------+  !
             !                                                                                      !
             ! CHECK UVs --------->>>
             do lvert = 2,npolyverts-1
                kvert = polyverts(lvert)
                do isurf = 1,2
                   call eval( &
                        xyzverif(:,isurf), &
                        curve%surf(isurf)%ptr, &
                        mesh%uv(:,isurf,kvert) )
                end do
                IF ( MAX(norm2(mesh%xyz(:,kvert) - xyzverif(:,1)), &
                     norm2(mesh%xyz(:,kvert) - xyzverif(:,2))) > EPSxyz ) THEN
                   PRINT *,'I =', kVERT,', ERR =', &
                        norm2(mesh%xyz(:,kvert) - xyzverif(:,1)), &
                        norm2(mesh%xyz(:,kvert) - xyzverif(:,2))
                END IF
             end do
             !
             !OPEN(UNIT=FID, FILE='../debug/polyline_uv.dat', ACTION='WRITE')
             !DO IPV = 1,CURVE%POLYLINE%NP
             !   WRITE (FID,*) CURVE%POLYLINE%UV(1:2,1:2,IPV)
             !END DO
             !CLOSE(FID)
             !PAUSE
             ! -----------<<<
             !
             !                                                                                      !
          end if ! <--------------------------------------------------------------------------------+
          !
          ! cyle halfedges
          ihedg = get_prev(ihedg)       ! ingoing
          ihedg = get_twin(mesh, ihedg) ! outgoing
          jface = get_face(ihedg)
          if ( jface < 1 .or. jface == iface ) exit
       end do
    end do
    deallocate(visited, bnd_hedg) 

    !IF ( .FALSE. ) THEN
    !   DO IVERT = 1,INTERDATA%NP
    !      POS => INTERDATA%POINTS(IVERT)%POS
    !      PRINT *,'POINT #',IVERT
    !      DO IPOS = 1,INTERDATA%POINTS(IVERT)%NPOS
    !         IF ( .NOT.ASSOCIATED(POS%SURF) ) PRINT *,'N/A'
    !         DO ISURF = 1,NSURF
    !            IF ( ASSOCIATED(POS%SURF, SURF(ISURF)) ) PRINT *,'   SURF #',ISURF,', UV =',POS%UV
    !         END DO
    !         POS => POS%NEXT
    !      END DO
    !   END DO
    !END IF
    
    ! make BREP
    brep%nf = 0
    brep%ne = 0
    brep%nv = 0
    do ivert = 1,interdata%np
       interdata%points(ivert)%ivert = 0
    end do
    call make_brep_from_intersection_data( &
         surf, &
         nsurf, &
         interdata, &
         brep )

    ! debugging >> ..................
       call write_brep_files( &
            brep, &
            trim(dir) // 'brep/verts.dat', &
            trim(dir) // 'brep/edges.dat', &
            trim(dir) // 'brep/faces.dat' )

       call write_intersection_data( &
            interdata, &
            trim(dir) // 'brep/intersection_points.dat', &
            trim(dir) // 'brep/intersection_curves.dat' )
       ! .................. <<
    
    do ivert = 1,mesh%nv ! <------------------------------------------------------+
       select case ( mesh%typ(ivert) ) ! <-----------------------------+          !
       case (0) ! -----------------------------------------------------+          !
          mesh%ids(ivert) = interdata%points(mesh%ids(ivert))%ivert    !          !
       case (1) ! -----------------------------------------------------+          !
          mesh%ids(ivert) = interdata%curves(mesh%ids(ivert))%iedge(1) !          !
       end select ! <--------------------------------------------------+          !
       !                                                                          !
       if ( mesh%typ(ivert) == 0 ) then ! <------------------------------------+  !
          nfeatedg = 0                                                         !  !
          ihedg = brep%verts(mesh%ids(ivert))%halfedge ! outgoing              !  !
          iface = get_face(brep, ihedg)                                        !  !
          do ! <--------------------------------------------------+            !  !
             if ( .not.is_smooth(brep, ihedg(1)) .or. &           !            !  !
                  is_boundary_edge(brep, ihedg(1)) ) then ! <--+  !            !  !
                nfeatedg = nfeatedg + 1                        !  !            !  !
                if ( nfeatedg == 1 ) then ! <--+               !  !            !  !
                   iedge = ihedg(1)            !               !  !            !  !
                end if ! <---------------------+               !  !            !  !
             end if ! <----------------------------------------+  !            !  !
             ihedg = get_prev(brep, ihedg) ! ingoing              !            !  !
             ihedg = get_twin(ihedg)       ! outgoing             !            !  !
             jface = get_face(brep, ihedg)                        !            !  !
             if ( jface == iface .or. jface < 1 ) exit            !            !  !
          end do ! <----------------------------------------------+            !  !
          !                                                                    !  !
          if ( nfeatedg == 0 ) then ! <-------------------------------------+  !  !
             mesh%typ(ivert) = 2                                            !  !  !
             jvert = mesh%ids(ivert)                                        !  !  !
             pos => brep%verts(jvert)%point%pos                             !  !  !
             mesh%ids(ivert) = iface                                        !  !  !
             do ipos = 1,brep%verts(jvert)%point%npos ! <------+            !  !  !
                if ( associated(pos%surf, &                    !            !  !  !
                     brep%faces(iface)%surface) ) then ! <--+  !            !  !  !
                   mesh%uv(1:2,1,ivert) = pos%uv            !  !            !  !  !
                   exit                                     !  !            !  !  !
                end if ! <----------------------------------+  !            !  !  !
                pos => pos%next                                !            !  !  !
             end do ! <----------------------------------------+            !  !  !
             !                                                              !  !  !
          elseif ( nfeatedg == 2 ) then ! ----------------------------------+  !  !
             mesh%typ(ivert) = 1                                            !  !  !
             jvert = mesh%ids(ivert)                                        !  !  !
             mesh%ids(ivert) = iedge                                        !  !  !
             do jface = 1,2 ! <------------------------------------------+  !  !  !
                iface = brep%edges(iedge)%halfedges(1+mod(jface,2))%face !  !  !  !
                pos => brep%verts(jvert)%point%pos                       !  !  !  !
                do ipos = 1,brep%verts(jvert)%point%npos ! <------+      !  !  !  !
                   if ( associated(pos%surf, &                    !      !  !  !  !
                        brep%faces(iface)%surface) ) then ! <--+  !      !  !  !  !
                      mesh%uv(1:2,jface,ivert) = pos%uv        !  !      !  !  !  !
                      exit                                     !  !      !  !  !  !
                   end if ! <----------------------------------+  !      !  !  !  !
                   pos => pos%next                                !      !  !  !  !
                end do ! <----------------------------------------+      !  !  !  !
             end do ! <--------------------------------------------------+  !  !  !
          end if ! <--------------------------------------------------------+  !  !
          !                                                                    !  !
       end if ! <--------------------------------------------------------------+  !
    end do ! <--------------------------------------------------------------------+

    IF ( .true. ) call check_uvs(brep, mesh)
    
    ! make hypergraph
    call make_hypergraph( &
         brep, &
         hypergraph, &
         feat_edge, &
         feat_vert )

    ! debugging >> ..................
    call get_free_unit( fid )
    open(unit=fid, file=trim(dir) // 'brep/hyperfaces.dat', action='write')
    write (fid,*) hypergraph%nhf
    do i = 1,hypergraph%nhf
       write (fid,*) hypergraph%hyperfaces(i)%nf
       write (fid,*) hypergraph%hyperfaces(i)%faces(1:hypergraph%hyperfaces(i)%nf)
    end do
    close(fid)

    open(unit=fid, file=trim(dir) // 'brep/hyperedges.dat', action='write')
    write (fid,*) hypergraph%nhe
    do i = 1,hypergraph%nhe
       write (fid,*) hypergraph%hyperedges(i)%ne
       write (fid,*) hypergraph%hyperedges(i)%verts
       write (fid,*) hypergraph%hyperedges(i)%hyperfaces
       do j = 1,hypergraph%hyperedges(i)%ne
          write (fid,*) hypergraph%hyperedges(i)%halfedges(1:2,j)
       end do
    end do
    close(fid)
    ! .................. <<

    do iface = 1,mesh%nt ! <------------------------------------+
       mesh%ihf(iface) = brep%faces(mesh%ihf(iface))%hyperface  !
    end do ! <--------------------------------------------------+

    do ihype = 1,hypergraph%nhe ! <--------------------------------------------------+
       path%nv = 0                                                                   !
       do iedge = 1,hypergraph%hyperedges(ihype)%ne ! <---------------------------+  !
          ihedg = hypergraph%hyperedges(ihype)%halfedges(1:2,iedge)               !  !
          do icurv = 1,interdata%nc ! <----------------------------------------+  !  !
             if ( interdata%curves(icurv)%iedge(1) == ihedg(1) ) then ! <---+  !  !  !
                if ( ihedg(2) == 2 ) then ! <----------------------------+  !  !  !  !
                   head = interdata%curves(icurv)%polyline%np            !  !  !  !  !
                   tail = 2                                              !  !  !  !  !
                   stride = -1                                           !  !  !  !  !
                else ! --------------------------------------------------+  !  !  !  !
                   head = 1                                              !  !  !  !  !
                   tail = interdata%curves(icurv)%polyline%np - 1        !  !  !  !  !
                   stride = 1                                            !  !  !  !  !
                end if ! <-----------------------------------------------+  !  !  !  !
                if ( iedge == hypergraph%hyperedges(ihype)%ne ) then ! <-+  !  !  !  !
                   tail = tail + stride                                  !  !  !  !  !
                end if ! <-----------------------------------------------+  !  !  !  !
                call insert_n_after( &                                      !  !  !  !
                     path%verts, &                                          !  !  !  !
                     path%nv, &                                             !  !  !  !
                     (tail - head)/stride + 1, &                            !  !  !  !
                     curve2verts(head:tail:stride,icurv), &                 !  !  !  !
                     path%nv )                                              !  !  !  !
                exit                                                        !  !  !  !
             end if ! <-----------------------------------------------------+  !  !  !
          end do ! <-----------------------------------------------------------+  !  !
       end do ! <-----------------------------------------------------------------+  !
       !                                                                             !
       path%hyperedge = ihype                                                        !
       call append_path( &                                                           !
            mesh, &                                                                  !
            path )                                                                   !
    end do ! <-----------------------------------------------------------------------+
    
    deallocate(curve2verts)


    call write_mesh_files( &
         mesh, &
         trim(dir) // 'mesh/tri.dat', &
         trim(dir) // 'mesh/xyz.dat', &
         trim(dir) // 'mesh/uv.dat', &
         trim(dir) // 'mesh/idstyp.dat', &
         trim(dir) // 'mesh/paths.dat' )


    ! padd tangent intersection polylines with new points
    IF ( .TRUE. ) THEN
       do i = 1,interdata%nc
          if ( interdata%curves(i)%smooth ) then
             PRINT *,'REFINE POLYLINE #',I
             call refine_intersection_polyline( &
                  interdata%curves(i), &
                  1.d-4, &
                  options%chord_err )
          end if
       end do
    END IF


  end subroutine init_from_mesh





  


  subroutine extract_intersection_polyline_vertices( &
       mesh, &
       ihedg_first, &
       visited, &
       bnd_hedg, &
       verts, &
       n )
    use mod_util
    use mod_mesh
    use mod_halfedge
    use mod_types_intersection
    implicit none
    type(type_surface_mesh),    intent(in)    :: mesh
    integer,                    intent(in)    :: ihedg_first(2)
    logical,                    intent(inout) :: visited(3,mesh%nt)
    logical,                    intent(in)    :: bnd_hedg(3,mesh%nt)
    integer, allocatable,       intent(inout) :: verts(:)
    integer,                    intent(out)   :: n
    integer                                   :: orig, dest
    integer, dimension(2)                     :: ihedg, jhedg

    ihedg = ihedg_first
    n = 0
    do ! <-------------------------------------------------+
       if ( .not.visited(ihedg(1),ihedg(2)) .and. &        !
            bnd_hedg(ihedg(1),ihedg(2)) ) then ! <-----+   !
          visited(ihedg(1),ihedg(2)) = .true.          !   !
          jhedg = get_twin(mesh, ihedg)                !   !
          if ( all(jhedg > 0) ) then ! <----------+    !   !
             visited(jhedg(1),jhedg(2)) = .true.  !    !   !
          end if ! <------------------------------+    !   !
          !                                            !   !
          orig = get_orig(mesh, ihedg)                 !   !
          call insert_after( &                         !   !
               verts, &                                !   !
               n, &                                    !   !
               orig, &                                 !   !
               n )                                     !   !
          !                                            !   !
          dest = get_orig(mesh, jhedg)                 !   !
          if ( mesh%typ(dest) == 0 ) then ! <---+      !   !
             call insert_after( &               !      !   !
                  verts, &                      !      !   !
                  n, &                          !      !   !
                  dest, &                       !      !   !
                  n )                           !      !   !
             return                             !      !   !
          end if ! <----------------------------+      !   !
          !                                            !   !
          ihedg = get_next(ihedg) ! ingoing            !   !
          cycle                                        !   !
       end if ! <--------------------------------------+   !
       ! cyle halfedges around current vertex              !
       ihedg = get_twin(mesh, ihedg) ! ingoing             !
       ihedg = get_next(ihedg)       ! outgoing            !
    end do ! <---------------------------------------------+
    
  end subroutine extract_intersection_polyline_vertices


  
end module mod_init
