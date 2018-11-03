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

    faces : do iface = 1,brep%nf
       ! write polynomial parametric surface 
       call write_polynomial(brep%faces(iface)%surface%x, '../tmp/c.cheb')
       !
       ! open boundary points & edges files
       call get_free_unit(fpts)
       open(unit=fpts, file='../tmp/bpts.dat', action='write')
       call get_free_unit(fedg)
       open(unit=fedg, file='../tmp/bedg.dat', action='write')
       !
       npf = 0
       wires : do iwire = 0,brep%faces(iface)%ninner ! <-----------------------------------+
          ! first halfedge of the wire                                                     !
          if ( iwire == 0 ) then ! <-------------------+                                   !
             ihedg = brep%faces(iface)%outer           !                                   !
          else ! --------------------------------------+                                   !
             ihedg = brep%faces(iface)%inner(:,iwire)  !                                   !
          end if ! <-----------------------------------+                                   !
          ifirstedge = ihedg(1)                                                            !
          !                                                                                !
          ifirstpoint = npf + 1                                                            !
          ! traverse the wire                                                              !
          halfedges : do ! <------------------------------------------------------------+  !
             ! get polyline points                                                      !  !
             call get_polyline_endpoints( &
                  brep, &
                  ihedg, &
                  ifirst, &
                  ilast, &
                  sens, &
                  np )
             !
             if ( .not.allocated(loc2glob) .or. &
                  size(loc2glob) < npf + np - 1 ) then ! <------+
                call reallocate_list(loc2glob, npf + np + 100)  !
             end if ! <-----------------------------------------+
             !
             if ( .not.allocated(idsf) .or. &
                  size(idsf) < npf + np - 1 ) then ! <----------+
                call reallocate_list(idsf, npf + np + 100)      !
             end if ! <-----------------------------------------+
             !
             if ( .not.allocated(typf) .or. &
                  size(typf) < npf + np - 1 ) then ! <----------+
                call reallocate_list(typf, npf + np + 100)      !
             end if ! <-----------------------------------------+
             !
             if ( allocated(uvf) .and. size(uvf,3) < np ) deallocate(uvf)
             if ( .not.allocated(uvf) ) allocate(uvf(2,2,np))
             ! polyligne uv dans le sens du contour de la face
             if ( sens == 2 ) then
                uvf(1:2,1:2,1:np) = brep%edges(ihedg(1))%curve%polyline%uv&
                     (1:2,[2,1],ifirst:ilast)
             else
                uvf(1:2,1:2,1:np) = brep%edges(ihedg(1))%curve%polyline%uv&
                     (1:2,1:2,ifirst:ilast:-1)
             end if
             !
             if ( allocated(xyzf) .and. size(xyzf,2) < np ) deallocate(xyzf)
             if ( .not.allocated(xyzf) ) allocate(xyzf(3,np))
             xyzf(1:3,1:np) = brep%edges(ihedg(1))%curve%polyline%xyz&
                  (1:3,ifirst:ilast:(-1)**sens)
             !                                                                          !  !
             ! get brep vertex index of first point
             ivert = get_orig(brep, ihedg)
             if ( bv2mv(ivert) == 0 ) then ! <-----------------------------------+
                ! the brep vertex is not yet present in the mesh                 !
                bv2mv(ivert) = mesh%nv + 1                                       !
                if ( feat_vert(ivert) ) then ! <------------------------------+  !
                   ! the mesh vertex is locked on a brep vertex               !  !
                   idsf(npf+1) = ivert                                        !  !
                   typf(npf+1) = 0                                            !  !
                else ! -------------------------------------------------------+  !
                   if ( feat_edge(ihedg(1)) ) then ! <---------------------+  !  !
                      ! the mesh vertex is locked on a hyperedge           !  !  !
                      idsf(npf+1) = ihedg(1)                               !  !  !
                      typf(npf+1) = 1                                      !  !  !
                   else ! -------------------------------------------------+  !  !
                      ! the mesh vertex is free to move along a hyperface  !  !  !
                      idsf(npf+1) = iface                                  !  !  !
                      typf(npf+1) = 2                                      !  !  !
                   end if ! <----------------------------------------------+  !  !
                end if ! <----------------------------------------------------+  !
                !                                                                !
                ! append new vertex to the mesh                                  !
                call append_vertices( &                                          !
                     mesh, &                                                     !
                     xyzf(1:3,1), &                                              !
                     uvf(1:2,1:2,1), &                                           !
                     idsf(npf+1), &                                              !
                     typf(npf+1), &                                              !
                     1 )                                                         !
             else ! -------------------------------------------------------------+
                if ( feat_edge(ihedg(1)) ) then ! <------------------------+     !
                   if ( mesh%typ(bv2mv(ivert)) == 2 ) then ! <----------+  !     !
                      mesh%ids(bv2mv(ivert)) = ihedg(1)                 !  !     !
                      mesh%typ(bv2mv(ivert)) = 1                        !  !     !
                      mesh%uv(1:2,1:2,bv2mv(ivert)) = uvf(1:2,[2,1],1)  !  !     !
                   end if ! <-------------------------------------------+  !     !
                end if ! <-------------------------------------------------+     !
             end if ! <----------------------------------------------------------+
             loc2glob(npf+1) = bv2mv(ivert)
             !
             if ( be2mv(1,ihedg(1)) == 0 ) then ! <-----------------------------------------+
                ! the brep edge is not yet present in the mesh                              !
                be2mv(1,ihedg(1)) = mesh%nv + 1                                             !
                be2mv(2,ihedg(1)) = mesh%nv + np - 2                                        !
                !                                                                           !
                if ( feat_edge(ihedg(1)) ) then ! <-----------------------+                 !
                   ! the mesh vertices are locked on a hyperedge          !                 !
                   idsf(npf+2:npf+np-1) = ihedg(1)                        !                 !
                   typf(npf+2:npf+np-1) = 1                               !                 !
                else ! ---------------------------------------------------+                 !
                   ! the mesh vertices are free to move along a hyperface !                 !
                   idsf(npf+2:npf+np-1) = iface                           !                 !
                   typf(npf+2:npf+np-1) = 2                               !                 !
                end if ! <------------------------------------------------+                 !
                !                                                                           !
                if ( sens == 1 ) then ! <-----------------------------------+               !
                   loc2glob(npf+2:npf+np-1) = mesh%nv + [(i, i=np-2,1,-1)]  !               !
                   call append_vertices( &                                  !               !
                        mesh, &                                             !               !
                        xyzf(1:3,np-1:2:-1), &                              !               !
                        uvf(1:2,1:2,np-1:2:-1), &                           !               !
                        idsf(npf+2:npf+np-1), &                             !               !
                        typf(npf+2:npf+np-1), &                             !               !
                        np-2 )                                              !               !
                else ! -----------------------------------------------------+               !
                   loc2glob(npf+2:npf+np-1) = mesh%nv + [(i, i=1,np-2)]     !               !
                   call append_vertices( &                                  !               !
                        mesh, &                                             !               !
                        xyzf(1:3,2:np-1), &                                 !               !
                        uvf(1:2,1:2,2:np-1), &                              !               !
                        idsf(npf+2:npf+np-1), &                             !               !
                        typf(npf+2:npf+np-1), &                             !               !
                        np-2 )                                              !               !
                end if ! <--------------------------------------------------+               !
             else ! ------------------------------------------------------------------------+
                if ( sens == 1 ) then ! <-----------------------------------+               !
                   loc2glob(npf+2:npf+np-1) = &                             !               !
                        [(i, i=be2mv(2,ihedg(1)),be2mv(1,ihedg(1)),-1)]     !               !
                else ! -----------------------------------------------------+               !
                   loc2glob(npf+2:npf+np-1) = &                             !               !
                        [(i, i=be2mv(1,ihedg(1)),be2mv(2,ihedg(1)))]        !               !
                end if ! <--------------------------------------------------+               !
             end if ! <---------------------------------------------------------------------+
             !
             ! write points...                            
             do i = 1,np-1 ! <---------------------------+
                write (fpts,*) uvf(:,1,i)                !
             end do ! <----------------------------------+
             ! ...and edges
             do i = 1,np-2 ! <---------------------------+
                write (fedg,*) npf + i, npf + i + 1      !
             end do ! <----------------------------------+
             !
             npf = npf + np - 1
             !
             ! move on to the next halfedge on the wire
             ihedg = get_next(brep, ihedg)
             !
             if ( ihedg(1) == ifirstedge ) then ! <------+
                ! the wire is complete                   !
                write (fedg,*) npf, ifirstpoint          !
                exit                                     !
             else ! -------------------------------------+
                write (fedg,*) npf, npf+1                !
             end if ! <----------------------------------+
             !           
          end do halfedges
       end do wires
       !
       ! close files
       close(fpts)
       close(fedg)
       !
       ! run mesher
       PRINT *,'MESH FACE #',IFACE
       call system( &!'/stck/bandrieu/Bureau/MeshGen/./meshgen.out &
            '/home/bastien/MeshGen/./meshgen.out &
            & ../tmp/c.cheb &
            & ../tmp/bpts.dat &
            & ../tmp/bedg.dat &
            & ../tmp/info.dat &
            & ../tmp/tri.dat &
            & ../tmp/uv.dat &
            & ../tmp/xyz.dat' )
       write (strnum,'(i3.3)') iface
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
          hyperedges(ihype)%verts(1) = ivert            !
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



  
end module mod_init
