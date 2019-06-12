program demo_EoS_brep

  use mod_util
  use mod_options
  use mod_math
  use mod_types_intersection
  use mod_types_brep
  use mod_brep
  use mod_hypergraph
  use mod_mesh
  use mod_intersection
  use mod_tolerances
  use mod_init
  use mod_optimmesh
  use mod_halfedge
  use mod_polynomial
  use mod_propagation
  use mod_eos

  !-------------------------------------------------
  implicit none
  character(10)                           :: argstr
  integer                                 :: argnum
  logical                                 :: MAKE_EOS
  
  type(type_options)                      :: options
  type(type_surface), allocatable, target :: surf(:)
  integer                                 :: nsurf, isurf, jsurf
  type(type_intersection_data), target    :: interdata
  integer                                 :: icurv
  integer, allocatable                    :: curvetype(:), edgetype(:), curvetypetmp(:)
  type(type_brep)                         :: brep
  integer                                 :: iface, jface, ivert, iedge
  type(type_hypergraph)                   :: hypergraph
  logical, allocatable                    :: feat_edge(:), feat_vert(:)
  type(type_surface_mesh)                 :: mesh

  type(type_surface), allocatable, target :: surf_new(:)
  integer                                 :: nsurf_new
  type(type_intersection_data), target    :: interdata_new
  type(type_brep)                         :: brep_new

  integer                                 :: fid
  character(3)                            :: strnum3

  integer                                 :: ihedg(2)

  real(kind=fp)                           :: xyzverif(3,2)

  type(ptr_surface), allocatable          :: surf_eos(:), edge_eos(:), vert_eos(:)
  integer                                 :: m, n, i, j, head, tail, sens, np
  real(kind=fp), allocatable              :: xyz(:,:,:), speed(:,:), s2e(:,:,:), tp(:)
  type(ptr_surface)                       :: enve_rl(2), surfpair(2)
  integer                                 :: jcurv, ipoint, jpoint

  real(kind=fp), allocatable              :: xyzll(:,:), uvlle(:,:), uvll(:,:), tmp(:,:), toto(:,:)
  integer                                 :: stat, nptot, nconvex, iborder
  integer, allocatable                    :: iconvex(:,:) ![iedge;iborder]
  type(type_surface)                      :: surfll
  real(kind=fp)                           :: longlat_range(2), solidangle
  real(kind=fp)                           :: uvendpoint(2,2), xyzendpoint(3)
  integer                                 :: idendpoint(2)
  type(type_intersection_curve), pointer  :: curve => null(), curve_new => null()

  real(kind=fp), allocatable              :: spheres(:,:)

  integer                                 :: parameterization_mode = 2
  !-------------------------------------------------


  !! Read argument
  if ( command_argument_count() < 1 ) then
     MAKE_EOS = .false.
  else
     call get_command_argument(1, argstr)
     read (argstr,*) argnum
     MAKE_EOS = (argnum > 0)
  end if



  

  call get_free_unit(fid)




  ! >>>>>-----------------
  IF ( .false. )  THEN
     open(unit=fid, file='demo_EoS_brep/debug/longlatpatch_xyz.dat', action='read')
     read (fid,*) n
     allocate(xyzll(3,n), uvll(2,n))
     do i = 1,n
        read (fid,*) xyzll(1:3,i)
     end do
     close(fid)

     call long_lat_patch_from_points( &
          [0._fp, 1._fp, 1._fp], &
          xyzll, &
          n, &
          16, &
          stat, &
          surfll, &
          uvll, &
          longlat_range, &
          solidangle )
     print *, longlat_range*180._fp/CSTpi, solidangle

     call write_polynomial( &
          surfll%x, &
          'demo_EoS_brep/debug/longlatpatch_surf.cheb' )

     open(unit=fid, file='demo_EoS_brep/debug/longlatpatch_uv.dat', action='write')
     do i = 1,n
        write (fid,'(es22.15,1x,es22.15)') uvll(1,i), uvll(2,i)
     end do
     close(fid)
     
     STOP
  END IF
  ! -----------------<<<<<


  !! read options file
  call read_options( &
       'demo_EoS_brep/demo_brep.opt', &
       options )

  !! import surfaces
  open( &
       unit = fid, &
       file = 'demo_EoS_brep/init/surftag.dat', &
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
          'demo_EoS_brep/init/coef/c_' // strnum3 // '.cheb', &
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

  !! read intersection data
  call read_intersection_curves_new( &
       'demo_EoS_brep/init/curves.dat', &
       surf, &
       interdata, &
       curvetype )
  do icurv = 1,interdata%nc
     interdata%curves(icurv)%smooth = (curvetype(icurv) == 0)
     allocate(interdata%curves(icurv)%iedge(interdata%curves(icurv)%nsplit-1))
     interdata%curves(icurv)%iedge(:) = 0
  end do

  call intersect_all_surfaces( &
          surf(1:nsurf), &
          nsurf, &
          interdata, &
          [((i == 2 .or. i == 8 .or. i == 9), i=1,nsurf)], &
          options%chord_err, &
          options%hmin, &
          options%hmax )

  call move_alloc(from=curvetype, to=curvetypetmp)
  allocate(curvetype(interdata%nc))
  curvetype(1:size(curvetypetmp)) = curvetypetmp(1:size(curvetypetmp))
  curvetype(size(curvetypetmp)+1:interdata%nc) = 1 ! concave intersection
  deallocate(curvetypetmp)
  PRINT *,'OK'

  ! --->>> DEBUG 
  call write_intersection_data( &
       interdata, &
       'demo_EoS_brep/debug/intersection_points.dat', &
       'demo_EoS_brep/debug/intersection_curves.dat' )
  ! <<<---


  !! make brep
  brep%nv = 0
  brep%ne = 0
  brep%nf = 0
  call make_brep_from_intersection_data( &
       surf, &
       nsurf, &
       interdata, &
       brep )

  ! --->>> DEBUG 
  call write_brep_files( &
       brep, &
       'demo_EoS_brep/debug/verts.dat', &
       'demo_EoS_brep/debug/edges.dat', &
       'demo_EoS_brep/debug/faces.dat' )
  call write_brep_edges_geometry( &
       brep, &
       'demo_EoS_brep/debug/edges_xyz.dat', &
       'demo_EoS_brep/debug/edges_uv.dat' )
  ! <<<---



  ! >>> ---------- Mark edges as smooth, concave or convex
  allocate(edgetype(brep%ne))
  do iedge = 1,brep%ne
     do icurv = 1,interdata%nc
        if ( associated(brep%edges(iedge)%curve, interdata%curves(icurv)) ) then
           edgetype(iedge) = curvetype(icurv)
           exit
        end if
     end do

     !call edge_convexity( &
     !     brep, &
     !     iedge, &
     !     convexity )
     !PRINT *, IEDGE, convexity, edgetype(iedge)
  end do
  !stop
  ! ----------<<<
  


IF ( .NOT.MAKE_EOS ) THEN

  IF ( .true. ) THEN
     do iface = 1,brep%nf
        PRINT *,'brepmesh face #',iface
        write (strnum3,'(i3.3)') iface
        call generate_face_mesh( &       
             brep, &
             iface, &
             'demo_EoS_brep/brepmesh/c_'//strnum3//'.cheb', &
             'demo_EoS_brep/brepmesh/bpts_'//strnum3//'.dat', &
             'demo_EoS_brep/brepmesh/bedg_'//strnum3//'.dat', &
             'demo_EoS_brep/brepmesh/info.dat', &
             'demo_EoS_brep/brepmesh/tri_'//strnum3//'.dat', &
             'demo_EoS_brep/brepmesh/uv_'//strnum3//'.dat', &
             'demo_EoS_brep/brepmesh/xyz_'//strnum3//'.dat' )
     end do
  END IF


  IF ( .true. ) THEN
     !! make hypergraph
     call make_hypergraph( &
          brep, &
          hypergraph, &
          feat_edge, &
          feat_vert )

     call write_hypergraph( &
          hypergraph, &
          'demo_EoS_brep/debug/hyperfaces.dat', &
          'demo_EoS_brep/debug/hyperedges.dat' )

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

     !do iface = 1,mesh%nt
     !   mesh%ihf(iface) = brep%faces(mesh%ihf(iface))%hyperface
     !end do


     call make_halfedges( &
          mesh )

     

     ! >>> ----- A DEBUGGER -----
     do ivert = 1,mesh%nv
        if ( mesh%typ(ivert) == 1 ) then
           do jface = 1,2
              iface = brep%edges(mesh%ids(ivert))%halfedges(1+mod(jface,2))%face
              call eval( &
                   xyzverif(:,jface), &
                   brep%faces(iface)%surface, &
                   mesh%uv(:,jface,ivert) )
           end do
           IF ( MAX(norm2(mesh%xyz(:,ivert) - xyzverif(:,1)), norm2(mesh%xyz(:,ivert) - xyzverif(:,2))) > EPSxyz ) THEN
              !PRINT *,mesh%ids(ivert)
              mesh%uv(1:2,1:2,ivert) = mesh%uv(1:2,[2,1],ivert) 
           END IF
        end if
     end do
     ! ----------<<<

     call check_uvs(brep, mesh)

     call write_connectivity( &
          'demo_EoS_brep/mesh/', &
          mesh, &
          0 )
     call write_xyz_positions( &
          'demo_EoS_brep/mesh/', &
          mesh, &
          0 )

     call write_inria_mesh( &
          mesh, &
          'demo_EoS_brep/mesh/mesh.mesh' )
     call write_obj_mesh( &
          mesh, &
          'demo_EoS_brep/mesh/mesh.obj' )
     call write_vtk_mesh( &
       mesh, &
       'demo_EoS_brep/mesh/mesh.vtk' )


     do iface = 1,mesh%nt
        mesh%ihf(iface) = brep%faces(mesh%ihf(iface))%hyperface
     end do
     call write_obj_mesh( &
          mesh, &
          'demo_EoS_brep/mesh/mesh_hyp.obj' )
     call write_vtk_mesh( &
       mesh, &
       'demo_EoS_brep/mesh/mesh_hyp.vtk' )

     if ( .true. ) then
        call optim_jiao_uv( &
             brep, &
             mesh, &
             1._fp, &!PARAM_frac_conf1, &
             0.7_fp, &!PARAM_frac_conf2, &
             0, &!PARAM_ipass1, &
             0, &!PARAM_ipass2, &
             20, &
             options%hmin, &
             options%hmax )
     else
        call optim_jiao( &
             brep, &
             hypergraph%hyperedges(1:hypergraph%nhe), &
             hypergraph%nhe, &
             mesh, &
             1._fp, &!PARAM_frac_conf1, &
             0.7_fp, &!PARAM_frac_conf2, &
             0, &!PARAM_ipass1, &
             0, &!PARAM_ipass2, &
             20, &
             options%hmin, &
             options%hmax )
     end if

     call write_mesh_files( &
          mesh, &
          'demo_EoS_brep/mesh/tri.dat', &
          'demo_EoS_brep/mesh/xyz.dat', &
          'demo_EoS_brep/mesh/uv.dat', &
          'demo_EoS_brep/mesh/idstyp.dat', &
          'demo_EoS_brep/mesh/paths.dat' )
     call get_free_unit(fid)
       open(unit=fid, file='demo_EoS_brep/mesh/mv2h.dat', action='write')
       do i = 1,mesh%nv
          write (fid,*) mesh%v2h(:,i)
       end do
       close(fid)
       open(unit=fid, file='demo_EoS_brep/mesh/mtwin.dat', action='write')
       do i = 1,mesh%nt
          write (fid,*) mesh%twin(:,:,i)
       end do
       close(fid)
       

     call write_xyz_positions( &
          'demo_EoS_brep/mesh/', &
          mesh, &
          1 )

     call write_inria_mesh( &
          mesh, &
          'demo_EoS_brep/mesh/mesh_optim.mesh' )
     call write_obj_mesh( &
          mesh, &
          'demo_EoS_brep/mesh/mesh_optim.obj' )
     call write_vtk_mesh( &
       mesh, &
       'demo_EoS_brep/mesh/mesh_optim.vtk' )


     IF ( .FALSE. ) THEN
        open(unit=fid, file='demo_EoS_brep/spheres.dat', action='write')
        do i = 1,mesh%nv
           write (fid,*) mesh%xyz(1:3,i), normal_speed(mesh%xyz(1,i),mesh%xyz(2,i),mesh%xyz(3,i),1,1)
        end do
        close(fid)
     ELSE
        call sample_spheres_brute_force( &
             mesh, &
             spheres )
        open(unit=fid, file='demo_EoS_brep/spheres.dat', action='write')
        do i = 1,size(spheres,2)
           write (fid,*) spheres(1:4,i)
        end do
        close(fid)
     END IF
     
     STOP
  END IF

END IF



  !call check_intersection_points( &
  !     interdata )
  !call check_intersection_polylines( &
  !     interdata )
  !PAUSE
  do isurf = 1,nsurf
     surf(isurf)%x%degr(1) = size(surf(isurf)%x%coef,1) - 1
     surf(isurf)%x%degr(2) = size(surf(isurf)%x%coef,2) - 1 
     call compute_deriv1(surf(isurf))
     call compute_deriv2(surf(isurf))
     call compute_pseudonormal(surf(isurf))
  end do
  

  allocate(surf_new(nsurf+brep%ne+brep%nv))
  ! >>> ---------- EoS from faces 
  allocate(surf_eos(nsurf))
  nsurf_new = 0
  do isurf = 1,nsurf
     PRINT *,'ISURF =',ISURF
     write (strnum3,'(i3.3)') isurf
     nsurf_new = nsurf_new + 1
     PRINT *,'NSURF_NEW =',NSURF_NEW
     
     m = 2*size(surf(isurf)%x%coef,1)
     n = 2*size(surf(isurf)%x%coef,2)
     !print *, m, n, size(surf(isurf)%x%coef,3)
     
     allocate(xyz(m,n,3), speed(m,n), s2e(m,n,3))
     call ifcht2( &
         surf(isurf)%x%coef, &
         xyz, &
         m, &
         n, &
         3 )
     
     speed = normal_speed( &
          xyz(1:m,1:n,1), &
          xyz(1:m,1:n,2), &
          xyz(1:m,1:n,3), &
          m, &
          n )
     
     call eos_surface2( &
          xyz, &
          speed, &
          m, &
          n, &
          1._fp, &
          1, &
          s2e )

     do i = 1,3
        xyz(1:m,1:n,i) = xyz(1:m,1:n,i) + 1._fp*speed(1:m,1:n)*s2e(1:m,1:n,i)
     end do

     call reset_polynomial(surf_new(nsurf_new)%x, 2, 1, [m-1,n-1], 3)
     call fcht2( &
          xyz, &
          surf_new(nsurf_new)%x%coef(1:m,1:n,1:3), &
          m, &
          n, &
          3, &
          EPSmath )
     surf_eos(isurf)%ptr => surf_new(nsurf_new)
     
     deallocate(xyz, speed, s2e)

     call economize2( &
          surf_eos(isurf)%ptr%x, &
          EPSmath )
     call compute_deriv1(surf_eos(isurf)%ptr)
     call compute_deriv2(surf_eos(isurf)%ptr)
     call compute_pseudonormal(surf_eos(isurf)%ptr)
     
     call write_polynomial( &
          surf_eos(isurf)%ptr%x, &
          'demo_EoS_brep/debug/eos_c_'//strnum3 //'.cheb' )
  end do
  ! ----------<<<

  

  ! >>> ---------- EoS from convex edges
  allocate(edge_eos(brep%ne))
  do iedge = 1,brep%ne
     if ( edgetype(iedge) == 1 ) cycle

     curve => brep%edges(iedge)%curve
     do jsurf = 1,2
        do isurf = 1,nsurf
           if ( associated(curve%surf(jsurf)%ptr, surf(isurf)) ) then
              enve_rl(jsurf)%ptr => surf_eos(isurf)%ptr
              exit
           end if
        end do
     end do

     if ( edgetype(iedge) == 0 ) then
        ! preserve smooth intersection curve
        np = curve%polyline%np
        do jpoint = 1,2
           ipoint = 1 + (jpoint-1)*(np-1)
           uvendpoint(1:2,1:2) = curve%polyline%uv(1:2,1:2,ipoint)
           call eval( &
                xyzendpoint(1:3), &
                enve_rl(1)%ptr, &
                uvendpoint(1:2,1) )
           call add_intersection_point( &
                uvendpoint, &
                xyzendpoint, &
                enve_rl, &
                2, &
                interdata_new, &
                idendpoint(jpoint) )
        end do

        call add_intersection_curve( &
               interdata_new, &
               [0._fp, 0._fp, 0._fp], &
               idendpoint, &
               spread(spread([-1._fp-EPSuv, 1._fp+EPSuv], 2, 2), 3, 2) )
        curve_new => interdata_new%curves(interdata_new%nc)
        curve_new%smooth = .true.
        curve_new%surf = enve_rl
        curve_new%isplit(2,1:2) = [1,np]
        allocate(curve_new%iedge(1))
        curve_new%iedge(:) = 0
        allocate(curve_new%polyline)
        curve_new%polyline%np = np
        allocate(curve_new%polyline%uv(2,2,np))
        curve_new%polyline%uv(1:2,1:2,1:np) = curve%polyline%uv(1:2,1:2,1:np)
        ! update xyz
        allocate(curve_new%polyline%xyz(3,np))
        do ipoint = 1,np
           call eval( &
                curve_new%polyline%xyz(1:3,ipoint), &
                curve_new%surf(1)%ptr, &
                curve_new%polyline%uv(1:2,1,ipoint) )
        end do
        cycle
     end if
     
     PRINT *,'IEDGE =',IEDGE
     write (strnum3,'(i3.3)') iedge
     nsurf_new = nsurf_new + 1
     PRINT *,'NSURF_NEW =',NSURF_NEW

     call get_polyline_endpoints( &
         brep, &
         [iedge,1], &
         head, &
         tail, &
         sens, &
         np )


     allocate(tp(np))
     call eos_from_curve( &
          parameterization_mode, &
          curve, &
          head, &
          tail, &
          enve_rl, &
          tp, &
          xyz, &
          m, &
          n )
     PRINT *,'M, N =', M, N

     call reset_polynomial(surf_new(nsurf_new)%x, 2, 1, [m-1,n-1], 3)
     call fcht2( &
          xyz, &
          surf_new(nsurf_new)%x%coef(1:m,1:n,1:3), &
          m, &
          n, &
          3, &
          EPSmath )

     edge_eos(iedge)%ptr => surf_new(nsurf_new)

     call economize2( &
          edge_eos(iedge)%ptr%x, &
          EPSmath )
     call compute_deriv1(edge_eos(iedge)%ptr)
     call compute_deriv2(edge_eos(iedge)%ptr)
     call compute_pseudonormal(edge_eos(iedge)%ptr)
     
     call write_polynomial( &
          edge_eos(iedge)%ptr%x, &
          'demo_EoS_brep/debug/eos_edge_c_'//strnum3//'.cheb' )

     ! 1st curve (R-E), 2nd curve (E-L)
     do jcurv = 1,2
        surfpair(jcurv)%ptr => enve_rl(jcurv)%ptr
        surfpair(1+mod(jcurv,2))%ptr => edge_eos(iedge)%ptr
        do jpoint = 1,2
           ipoint = 1 + (jpoint-1)*(np-1)
           uvendpoint(1:2,jcurv) = curve%polyline%uv(1:2,jcurv,ipoint)
           uvendpoint(1,1+mod(jcurv,2)) = real((-1)**jpoint, kind=fp)
           uvendpoint(2,1+mod(jcurv,2)) = real((-1)**jcurv, kind=fp)
           call eval( &
                xyzendpoint(1:3), &
                surfpair(1)%ptr, &
                uvendpoint(1:2,1) )
           call add_intersection_point( &
                uvendpoint, &
                xyzendpoint, &
                surfpair, &
                2, &
                interdata_new, &
                idendpoint(jpoint) )
        end do

        call add_intersection_curve( &
               interdata_new, &
               [0._fp, 0._fp, 0._fp], &
               idendpoint, &
               spread(spread([-1._fp-EPSuv, 1._fp+EPSuv], 2, 2), 3, 2) )
        curve_new => interdata_new%curves(interdata_new%nc)
        curve_new%smooth = .true.
        curve_new%surf = surfpair
        curve_new%isplit(2,1:2) = [1,np]
        allocate(curve_new%iedge(1))
        curve_new%iedge(:) = 0
        allocate(curve_new%polyline)
        curve_new%polyline%np = np
        allocate(curve_new%polyline%uv(2,2,np))
        curve_new%polyline%uv(1:2,jcurv,1:np) = curve%polyline%uv(1:2,jcurv,1:np)
        curve_new%polyline%uv(1,1+mod(jcurv,2),1:np) = tp(1:np)
        curve_new%polyline%uv(2,1+mod(jcurv,2),1:np) = real((-1)**jcurv, kind=fp)
        ! update xyz
        allocate(curve_new%polyline%xyz(3,np))
        do ipoint = 1,np
           call eval( &
                curve_new%polyline%xyz(1:3,ipoint), &
                curve_new%surf(1)%ptr, &
                curve_new%polyline%uv(1:2,1,ipoint) )
        end do
     end do
     
     deallocate(tp, xyz)
  end do
  ! ----------<<<

  
  
  
  
  ! >>> ---------- EoS from convex vertices
  IF ( .true. ) THEN
     allocate(vert_eos(brep%nv))
     np = 30
     allocate(xyzll(3,np), uvlle(2,np))
     do ivert = 1,brep%nv
        PRINT *,'VERTEX #',ivert
        ihedg = brep%verts(ivert)%halfedge
        iface = get_face(brep, ihedg)
        nconvex = 0
        nptot = 0
        do
           if ( edgetype(ihedg(1)) == 2 ) then
              if ( ihedg(2) == 1 ) then
                 iborder = 4
              else
                 iborder = 2
              end if

              if ( .not.allocated(xyzll) ) allocate(xyzll(3,nptot+np))
              if ( .not.allocated(uvlle) ) allocate(uvlle(2,nptot+np))

              if ( size(xyzll,2) < nptot+np ) then
                 call move_alloc(from=xyzll, to=tmp)
                 allocate(xyzll(3,nptot+np))
                 xyzll(1:3,1:nptot) = tmp(1:3,1:nptot)
                 deallocate(tmp)
              end if
              if ( size(uvlle,2) < nptot+np ) then
                 call move_alloc(from=uvlle, to=tmp)
                 allocate(uvlle(2,nptot+np))
                 uvlle(1:2,1:nptot) = tmp(1:2,1:nptot)
                 deallocate(tmp)
              end if

              allocate(tmp(5,np))
              call trace_border_polyline( &
                   edge_eos(ihedg(1))%ptr, &
                   iborder, &
                   np, &
                   tmp(1:2,1:np), &
                   tmp(3:5,1:np) )
              xyzll(1:3,nptot+1:nptot+np) = tmp(3:5,1:np)
              uvlle(1:2,nptot+1:nptot+np) = tmp(1:2,1:np)

              call append_vec_i( &
                   [ihedg(1),iborder,nptot+1,nptot+np], &
                   4, &
                   iconvex, &
                   nconvex )

              nptot = nptot + np
              deallocate(tmp)
           end if

           ! traverse halfedges counter-clockwise
           ihedg = get_prev(brep, ihedg) ! previous halfedge
           ihedg = get_twin(ihedg) ! twin halfedge
           jface = get_face(brep, ihedg)

           if ( jface == iface ) exit
        end do

        if ( nconvex == 2 ) then
           do iedge = 1,2
              surfpair(iedge)%ptr => edge_eos(iconvex(1,1+mod(iedge,2)))%ptr
           end do

           ! for now we assume np is equal for both incident curves' polylines
           ! (=> in particular ntot = 2*np)
           head = iconvex(3,1)
           tail = iconvex(4,1)
           np = tail - head + 1
           do jpoint = 1,2
              ipoint = 1 + (jpoint-1)*(np-1)
              uvendpoint(1:2,1) = uvlle(1:2,2*np-1+ipoint)
              uvendpoint(1:2,2) = uvlle(1:2,ipoint)
              call add_intersection_point( &
                   uvendpoint, &
                   xyzll(:,ipoint), &
                   surfpair, &
                   2, &
                   interdata_new, &
                   idendpoint(jpoint) )
           end do

           call add_intersection_curve( &
                interdata_new, &
                [0._fp, 0._fp, 0._fp], &
                idendpoint, &
                spread(spread([-1._fp-EPSuv, 1._fp+EPSuv], 2, 2), 3, 2) )
           curve_new => interdata_new%curves(interdata_new%nc)
           curve_new%smooth = .true.
           curve_new%surf = surfpair
           curve_new%isplit(2,1:2) = [1,np]
           allocate(curve_new%iedge(1))
           curve_new%iedge(:) = 0
           allocate(curve_new%polyline)
           curve_new%polyline%np = np
           allocate(curve_new%polyline%uv(2,2,np), curve_new%polyline%xyz(3,np))
           curve_new%polyline%uv(1:2,1,1:np) = uvlle(1:2,2*np:np+1:-1)
           curve_new%polyline%uv(1:2,2,1:np) = uvlle(1:2,1:np)
           curve_new%polyline%xyz(1:3,1:np) = &
                0.5_fp*(xyzll(1:3,1:np) + xyzll(1:3,2*np:np+1:-1))

        elseif ( nconvex > 2 ) then

           write (strnum3,'(i3.3)') ivert
           nsurf_new = nsurf_new + 1
           PRINT *,'NSURF_NEW =',NSURF_NEW

           allocate(uvll(2,nptot))
           ALLOCATE(TOTO(3,NPTOT))
           TOTO(1:3,1:NPTOT) = XYZLL(1:3,1:NPTOT)
           call long_lat_patch_from_points( &
                brep%verts(ivert)%point%xyz, &
                xyzll(1:3,1:nptot), &
                nptot, &
                16, &
                stat, &
                surf_new(nsurf_new), &
                uvll, &
                longlat_range, &
                solidangle )

           vert_eos(ivert)%ptr => surf_new(nsurf_new)
           call economize2( &
                vert_eos(ivert)%ptr%x, &
                EPSmath )
           call compute_deriv1(vert_eos(ivert)%ptr)
           call compute_deriv2(vert_eos(ivert)%ptr)
           call compute_pseudonormal(vert_eos(ivert)%ptr)

           call write_polynomial( &
                vert_eos(ivert)%ptr%x, &
                'demo_EoS_brep/debug/eos_vert_c_'//strnum3//'.cheb' )

           ALLOCATE(TMP(3,NPTOT))
           DO I = 1,NPTOT
              CALL EVAL( &
                   TMP(1:3,I), &
                   VERT_EOS(IVERT)%PTR, &
                   UVLL(1:2,I) )
           end do
           PRINT *,SQRT(MAXVAL(SUM((TMP(1:3,1:NP) - TOTO(1:3,1:NP))**2,1)))
           DEALLOCATE(TMP, TOTO)

           surfpair(1)%ptr => vert_eos(ivert)%ptr
           do iedge = 1,nconvex
              surfpair(2)%ptr => edge_eos(iconvex(1,iedge))%ptr
              do i = 1,2
                 j = iconvex(2+i,iedge)
                 uvendpoint(1:2,1) = uvll(1:2,j)
                 uvendpoint(1:2,2) = uvlle(1:2,j)
                 call add_intersection_point( &
                      uvendpoint, &
                      xyzll(1:3,j), &
                      surfpair, &
                      2, &
                      interdata_new, &
                      idendpoint(i) ) 
              end do

              call add_intersection_curve( &
                   interdata_new, &
                   [0._fp, 0._fp, 0._fp], &
                   idendpoint, &
                   spread(spread([-1._fp-EPSuv, 1._fp+EPSuv], 2, 2), 3, 2) )

              curve_new => interdata_new%curves(interdata_new%nc)
              curve_new%smooth = .true.
              curve_new%surf = surfpair
              head = iconvex(3,iedge)
              tail = iconvex(4,iedge)
              np = tail - head + 1
              curve_new%isplit(2,1:2) = [1,np]
              allocate(curve_new%iedge(1))
              curve_new%iedge(:) = 0
              allocate(curve_new%polyline)
              curve_new%polyline%np = np
              allocate(curve_new%polyline%uv(2,2,np), curve_new%polyline%xyz(3,np))
              curve_new%polyline%uv(1:2,1,1:np) = uvll(1:2,head:tail)
              curve_new%polyline%uv(1:2,2,1:np) = uvlle(1:2,head:tail)
              curve_new%polyline%xyz(1:3,1:np) = xyzll(1:3,head:tail)
           end do
        end if

        if ( allocated(uvll) ) deallocate(uvll)

     end do
     if ( allocated(xyzll) ) deallocate(xyzll)
  END IF
  ! ----------<<<

  call check_intersection_points( &
       interdata_new )
  !STOP
  call check_intersection_polylines( &
       interdata_new )
  !STOP

  
  IF ( .true. ) THEN
     call intersect_all_surfaces( &
          surf_new(1:nsurf_new), &
          nsurf_new, &
          interdata_new, &
          [((i == 2 .or. i == 8 .or. i == 9), i=1,nsurf_new)], &
          options%chord_err, &
          options%hmin, &
          options%hmax )
  END IF

  
  ! --->>> DEBUG 
  call write_intersection_data( &
       interdata_new, &
       'demo_EoS_brep/debug/intersection_points_new.dat', &
       'demo_EoS_brep/debug/intersection_curves_new.dat' )
  ! <<<---

  !STOP
  
  !! make new brep
  brep_new%nv = 0
  brep_new%ne = 0
  brep_new%nf = 0
  call make_brep_from_intersection_data( &
       surf_new, &
       nsurf_new, &
       interdata_new, &
       brep_new )


  ! --->>> DEBUG 
  call write_brep_files( &
       brep_new, &
       'demo_EoS_brep/debug/verts_new.dat', &
       'demo_EoS_brep/debug/edges_new.dat', &
       'demo_EoS_brep/debug/faces_new.dat' )
  open(unit=fid, file='demo_EoS_brep/debug/edges_xyz_new.dat', action='write')
  write (fid,*) brep_new%ne
  do iedge = 1,brep_new%ne
     call get_polyline_endpoints( &
          brep_new, &
          [iedge,1], &
          head, &
          tail, &
          sens, &
          np )
     write (fid,*) np
     do ivert = head,tail,(-1)**sens
        write (fid,*) brep_new%edges(iedge)%curve%polyline%xyz(1:3,ivert)
     end do
  end do
  close(fid)
  ! <<<---


  IF ( .true. ) THEN
     do iface = 1,brep_new%nf
        PRINT *,'brepmesh face #',iface
        write (strnum3,'(i3.3)') iface
        call generate_face_mesh( &       
             brep_new, &
             iface, &
             'demo_EoS_brep/brepmesh_eos/c_'//strnum3//'.cheb', &
             'demo_EoS_brep/brepmesh_eos/bpts_'//strnum3//'.dat', &
             'demo_EoS_brep/brepmesh_eos/bedg_'//strnum3//'.dat', &
             'demo_EoS_brep/brepmesh_eos/info.dat', &
             'demo_EoS_brep/brepmesh_eos/tri_'//strnum3//'.dat', &
             'demo_EoS_brep/brepmesh_eos/uv_'//strnum3//'.dat', &
             'demo_EoS_brep/brepmesh_eos/xyz_'//strnum3//'.dat' )
     end do
  END IF




  IF ( .true. ) THEN
     !! make hypergraph
     call make_hypergraph( &
          brep_new, &
          hypergraph, &
          feat_edge, &
          feat_vert )

     call generate_brep_conforming_mesh( &
          brep_new, &
          options%hmin, &
          options%hmax, &
          options%chord_err, &
          feat_edge(1:brep%ne), &
          feat_vert(1:brep%nv), &
          hypergraph%hyperedges(1:hypergraph%nhe), &
          hypergraph%nhe, &
          mesh )

     !do iface = 1,mesh%nt
     !   mesh%ihf(iface) = brep%faces(mesh%ihf(iface))%hyperface
     !end do

     call write_inria_mesh( &
          mesh, &
          'demo_EoS_brep/mesh_eos/mesh_eos.mesh' )
     call write_obj_mesh( &
          mesh, &
          'demo_EoS_brep/mesh_eos/mesh_eos.obj' )
     call write_vtk_mesh( &
       mesh, &
       'demo_EoS_brep/mesh_eos/mesh_eos.vtk' )

     call make_halfedges( &
          mesh )

     ! >>> ----- A DEBUGGER -----
     do ivert = 1,mesh%nv
        if ( mesh%typ(ivert) == 1 ) then
           do jface = 1,2
              iface = brep_new%edges(mesh%ids(ivert))%halfedges(1+mod(jface,2))%face
              call eval( &
                   xyzverif(:,jface), &
                   brep_new%faces(iface)%surface, &
                   mesh%uv(:,jface,ivert) )
           end do
           IF ( MAX(norm2(mesh%xyz(:,ivert) - xyzverif(:,1)), norm2(mesh%xyz(:,ivert) - xyzverif(:,2))) > EPSxyz ) THEN
              !PRINT *,mesh%ids(ivert)
              mesh%uv(1:2,1:2,ivert) = mesh%uv(1:2,[2,1],ivert) 
           END IF
        end if
     end do
     ! ----------<<<

     call check_uvs(brep_new, mesh)

     call write_connectivity( &
          'demo_EoS_brep/mesh_eos/', &
          mesh, &
          0 )
     call write_xyz_positions( &
          'demo_EoS_brep/mesh_eos/', &
          mesh, &
          0 )

     do iface = 1,mesh%nt
        mesh%ihf(iface) = brep_new%faces(mesh%ihf(iface))%hyperface
     end do
     
     call write_obj_mesh( &
          mesh, &
          'demo_EoS_brep/mesh_eos/mesh_eos_hyp.obj' )
     call write_vtk_mesh( &
       mesh, &
       'demo_EoS_brep/mesh_eos/mesh_eos_hyp.vtk' )

     if ( .true. ) then
        call optim_jiao_uv( &
             brep_new, &
             mesh, &
             1._fp, &!PARAM_frac_conf1, &
             0.7_fp, &!PARAM_frac_conf2, &
             0, &!PARAM_ipass1, &
             0, &!PARAM_ipass2, &
             20, &
             options%hmin, &
             options%hmax )
     else
        call optim_jiao( &
             brep_new, &
             hypergraph%hyperedges(1:hypergraph%nhe), &
             hypergraph%nhe, &
             mesh, &
             1._fp, &!PARAM_frac_conf1, &
             0.7_fp, &!PARAM_frac_conf2, &
             0, &!PARAM_ipass1, &
             0, &!PARAM_ipass2, &
             20, &
             options%hmin, &
             options%hmax )
     end if
  
     call write_xyz_positions( &
          'demo_EoS_brep/mesh_eos/', &
          mesh, &
          1 )


     call write_inria_mesh( &
          mesh, &
          'demo_EoS_brep/mesh_eos/mesh_eos_optim.mesh' )
     call write_obj_mesh( &
          mesh, &
          'demo_EoS_brep/mesh_eos/mesh_eos_optim.obj' )
     call write_vtk_mesh( &
          mesh, &
          'demo_EoS_brep/mesh_eos/mesh_eos_optim.vtk' )
     
     STOP
  END IF

  

  
  STOP

  if ( .true. ) then
     call optim_jiao_uv( &
          brep, &
          mesh, &
          1._fp, &!PARAM_frac_conf1, &
          0.7_fp, &!PARAM_frac_conf2, &
          0, &!PARAM_ipass1, &
          0, &!PARAM_ipass2, &
          20, &
          options%hmin, &
          options%hmax )
  else
     call optim_jiao( &
          brep, &
          hypergraph%hyperedges(1:hypergraph%nhe), &
          hypergraph%nhe, &
          mesh, &
          1._fp, &!PARAM_frac_conf1, &
          0.7_fp, &!PARAM_frac_conf2, &
          0, &!PARAM_ipass1, &
          0, &!PARAM_ipass2, &
          20, &
          options%hmin, &
          options%hmax )
  end if

  call write_xyz_positions( &
       'demo_EoS_brep/mesh/', &
       mesh, &
       1 )

contains

  subroutine read_intersection_curves_new( &
       filename, &
       surf, &
       interdata, &
       curvetype )
    use mod_util
    use mod_diffgeom
    use mod_types_intersection
    use mod_tolerances
    implicit none
    character(*),                 intent(in)    :: filename
    type(type_surface), target,   intent(in)    :: surf(:)
    type(type_intersection_data), intent(inout) :: interdata
    integer, allocatable,         intent(inout) :: curvetype(:)
    integer                                     :: fid, ncurves, np
    real(kind=fp), allocatable                  :: xyz(:,:), uv(:,:,:)
    type(ptr_surface)                           :: surfpair(2)
    integer                                     :: ic, ip, pair(2), endpt(2)

    call get_free_unit(fid)
    open( &
         unit=fid, &
         file=filename, &
         action='read' )
    read (fid,*) ncurves
    allocate(curvetype(ncurves))
    do ic = 1,ncurves
       read (fid,*) pair
       surfpair(1)%ptr => surf(pair(1))
       surfpair(2)%ptr => surf(pair(2))
       read (fid,*) curvetype(ic)
       read (fid,*) np

       allocate(xyz(3,np), uv(2,2,np))
       do ip = 1,np
          read (fid,*) xyz(:,ip)
       end do
       do ip = 1,np
          read (fid,*) uv(:,1,ip), uv(:,2,ip)
       end do

       call add_intersection_point( &
            uv(:,:,1), &
            xyz(:,1), &
            surfpair, &
            2, &
            interdata, &
            endpt(1) ) 
       call add_intersection_point( &
            uv(:,:,np), &
            xyz(:,np), &
            surfpair, &
            2, &
            interdata, &
            endpt(2) )

       call add_intersection_curve( &
            interdata, &
            [0._fp, 0._fp, 0._fp], &
            endpt, &
            spread(spread([-1._fp-EPSuv, 1._fp+EPSuv], 2, 2), 3, 2) )
       interdata%curves(interdata%nc)%surf(1)%ptr => surf(pair(1))
       interdata%curves(interdata%nc)%surf(2)%ptr => surf(pair(2))
       interdata%curves(interdata%nc)%isplit(2,1:2) = [1,np]

       allocate(interdata%curves(interdata%nc)%polyline)
       interdata%curves(interdata%nc)%polyline%np = np
       call move_alloc(from=xyz, to=interdata%curves(interdata%nc)%polyline%xyz)
       call move_alloc(from=uv , to=interdata%curves(interdata%nc)%polyline%uv )

    end do
    close(fid)

  end subroutine read_intersection_curves_new




  subroutine check_intersection_points( &
       interdata )
    implicit none
    type(type_intersection_data), intent(in) :: interdata
    type(type_point_on_surface), pointer     :: pos => null()
    real(kind=fp)                            :: xyz(3), err
    integer                                  :: ip, is

    do ip = 1,interdata%np
       pos => interdata%points(ip)%pos
       do is = 1,interdata%points(ip)%npos
          call eval( &
               xyz, &
               pos%surf, &
               pos%uv )
          err = norm2(xyz - interdata%points(ip)%xyz)
          if ( err > 1.d-7 ) then
             print *, 'POINT #', ip, is, err!, 'U=',pos%uv,', X=',xyz!, &
             !'/', interdata%points(ip)%xyz
             !print *, pos%surf%x%coef(1,2,1:3)
          end if
          pos => pos%next
       end do
    end do
  end subroutine check_intersection_points




  subroutine check_intersection_polylines( &
       interdata )
    implicit none
    type(type_intersection_data), intent(in) :: interdata
    real(kind=fp)                            :: xyz(3,2), err(2)
    integer                                  :: ic, ip, is

    do ic = 1,interdata%nc
       do ip = 1,interdata%curves(ic)%polyline%np
          do is = 1,2
             call eval( &
                  xyz(1:3,is), &
                  interdata%curves(ic)%surf(is)%ptr, &
                  interdata%curves(ic)%polyline%uv(1:2,is,ip) )
             err(is) = norm2(xyz(1:3,is) - &
                  interdata%curves(ic)%polyline%xyz(1:3,ip))
          end do
          if ( any(err > 1.d-7) ) then
             PRINT *, 'CURVE #',ic, ip, err
             do is = 1,2
                call eval( &
                     xyz(1:3,is), &
                     interdata%curves(ic)%surf(1+mod(is,2))%ptr, &
                     interdata%curves(ic)%polyline%uv(1:2,is,ip) )
                err(is) = norm2(xyz(1:3,is) - &
                     interdata%curves(ic)%polyline%xyz(1:3,ip))
             end do
             PRINT *, err
          end if
       end do
    end do
    
  end subroutine check_intersection_polylines
       





  subroutine sample_spheres( &
       mesh, &
       spheres )
    implicit none
    type(type_surface_mesh),    intent(in)    :: mesh
    real(kind=fp), allocatable, intent(inout) :: spheres(:,:)
    real(kind=fp)                             :: r(1,mesh%nv)
    integer                                   :: seeds(mesh%nv), nseeds    
    logical, dimension(mesh%nv)               :: active, visited
    integer                                   :: ivert
    
    r = normal_speed( &
         mesh%xyz(1,1:mesh%nv), &
         mesh%xyz(2,1:mesh%nv), &
         mesh%xyz(3,1:mesh%nv), &
         1, &
         mesh%nv )
    
    nseeds = 0
    active(1:mesh%nv) = .true.
    do while ( any(active) )
       do ivert = 1,mesh%nv
          if ( active(ivert) ) then
             print *,count(active),'/',mesh%nv
             print *,'+1 seed:',ivert
             visited = .not.active
             nseeds = nseeds + 1
             seeds(nseeds) = ivert
             call grow_seed( &
                  ivert, &
                  mesh%xyz(1:3,ivert), &
                  r(1,ivert), &
                  mesh, &
                  r, &
                  active, &
                  visited )
             !exit
          end if
       end do
    end do

    allocate(spheres(4,nseeds))
    do ivert = 1,nseeds
       spheres(1:3,ivert) = mesh%xyz(1:3,seeds(ivert))
       spheres(4,ivert) = r(1,seeds(ivert))
    end do
    
  end subroutine sample_spheres




  recursive subroutine grow_seed( &
       ivert, &
       xyz_seed, &
       r_seed, &
       mesh, &
       r, &
       active, &
       visited )
    implicit none
    integer,                 intent(in)    :: ivert
    real(kind=fp),           intent(in)    :: xyz_seed(3)
    real(kind=fp),           intent(in)    :: r_seed
    type(type_surface_mesh), intent(in)    :: mesh
    real(kind=fp),           intent(in)    :: r(mesh%nv)
    logical,                 intent(inout) :: active(mesh%nv)
    logical,                 intent(inout) :: visited(mesh%nv)
    integer, dimension(2)                  :: ihedg, jhedg
    integer                                :: iface, jface, jvert
    !INTEGER :: K

    visited(ivert) = .true.
    !PRINT *,'   IVERT =',IVERT,norm2(xyz_seed - mesh%xyz(1:3,ivert)) - (r_seed + r(ivert))
    !PRINT *,'   ',count(active),'/',mesh%nv
    
    if ( norm2(xyz_seed - mesh%xyz(1:3,ivert)) > r_seed + r(ivert) ) return

    active(ivert) = .false.
    
    ihedg = mesh%v2h(1:2,ivert) ! outgoing
    iface = get_face(ihedg)
    !PRINT *,'   IFACE =',IFACE
    !K = 0
    do
       jvert = get_dest(mesh, ihedg)
       !PRINT *,'AVANT :',IVERT,JVERT,active(jvert),visited(jvert)
       if ( active(jvert) .and. .not.visited(jvert) ) then
          call grow_seed( &
               jvert, &
               xyz_seed, &
               r_seed, &
               mesh, &
               r, &
               active, &
               visited )
       end if
       !PRINT *,'APRES :',IVERT,JVERT,active(jvert),visited(jvert)
       jhedg = get_prev(ihedg) ! ingoing
       ihedg = get_twin(mesh, jhedg) ! outgoing
       jface = get_face(ihedg)
       !PRINT *,JVERT,JFACE,ihedg
       if ( jface == iface .or. jface < 1 ) exit
       !K = K + 1
       !IF ( K > 50 ) STOP
    end do
    
  end subroutine grow_seed
  









  subroutine get_v2v( &
       mesh, &
       ivert, &
       adjv, &
       nadjv, &
       nadjvmax )
    use mod_halfedge
    use mod_mesh
    implicit none
    type(type_surface_mesh), intent(in)  :: mesh
    integer,                 intent(in)  :: ivert
    integer,                 intent(in)  :: nadjvmax
    integer,                 intent(out) :: adjv(nadjvmax)
    integer,                 intent(out) :: nadjv
    integer, dimension(2)                :: ihedg, jhedg
    integer                              :: iface, jface, jvert

    ihedg = mesh%v2h(1:2,ivert)
    iface = get_face(ihedg)
    nadjv = 0
    do
       jvert = get_dest(mesh, ihedg)
       if ( nadjv + 1 > nadjvmax ) STOP 'nadjv + 1 > nadjvmax'
       nadjv = nadjv + 1
       adjv(nadjv) = jvert
       
       jhedg = get_prev(ihedg)       ! ingoing
       ihedg = get_twin(mesh, jhedg) ! outgoing
       jface = get_face(ihedg)
       
       if ( jface == iface .or. jface < 1 ) exit
    end do
    
  end subroutine get_v2v







  subroutine sample_spheres2( &
       mesh, &
       spheres )
    implicit none
    type(type_surface_mesh),    intent(in)    :: mesh
    real(kind=fp), allocatable, intent(inout) :: spheres(:,:)
    real(kind=fp)                             :: r(1,mesh%nv)
    integer                                   :: seeds(mesh%nv), nseeds
    integer                                   :: front(mesh%nv), nfront
    logical, dimension(mesh%nv)               :: active, visited
    integer                                   :: ivert, jvert
    
    r = normal_speed( &
         mesh%xyz(1,1:mesh%nv), &
         mesh%xyz(2,1:mesh%nv), &
         mesh%xyz(3,1:mesh%nv), &
         1, &
         mesh%nv )
    
    nseeds = 0
    active(1:mesh%nv) = .true.
    do while ( any(active) )
       print *,count(active),'/',mesh%nv
       do ivert = 1,mesh%nv
          if ( .not.active(ivert) ) cycle

          print *,'+1 seed:',ivert
          visited = .not.active
          nseeds = nseeds + 1
          seeds(nseeds) = ivert

          call get_v2v( &
               mesh, &
               ivert, &
               front, &
               nfront, &
               mesh%nv )
          visited(front(1:nfront)) = .true.
          do jvert = 1,nfront
             if ( norm2(mesh%xyz(1:3,jvert) - mesh%xyz(1:3,ivert)) > r(1,jvert) + r(1,ivert) ) then
                !candidate = jvert
             end if
          end do
          

          
          call grow_seed( &
               ivert, &
               mesh%xyz(1:3,ivert), &
               r(1,ivert), &
               mesh, &
               r, &
               active, &
               visited )
       end do
    end do

    allocate(spheres(4,nseeds))
    do ivert = 1,nseeds
       spheres(1:3,ivert) = mesh%xyz(1:3,seeds(ivert))
       spheres(4,ivert) = r(1,seeds(ivert))
    end do
    
  end subroutine sample_spheres2




  subroutine sample_spheres_brute_force( &
       mesh, &
       spheres )
    implicit none
    type(type_surface_mesh),    intent(in)    :: mesh
    real(kind=fp), allocatable, intent(inout) :: spheres(:,:)
    real(kind=fp)                             :: r(1,mesh%nv)
    integer                                   :: perm(mesh%nv)
    integer                                   :: seeds(mesh%nv), nseeds
    logical, dimension(mesh%nv)               :: active
    integer                                   :: i, j, ivert, jvert

    call randperm(perm, mesh%nv)

    r = 0.25_fp*normal_speed( &
         mesh%xyz(1,1:mesh%nv), &
         mesh%xyz(2,1:mesh%nv), &
         mesh%xyz(3,1:mesh%nv), &
         1, &
         mesh%nv )

    nseeds = 0
    active(1:mesh%nv) = .true.
    do i = 1,mesh%nv
       ivert = perm(i)
       if ( .not.active(ivert) ) cycle
       nseeds = nseeds + 1
       seeds(nseeds) = ivert
       
       do j = i+1,mesh%nv
          jvert = perm(j)
          if ( .not.active(jvert) ) cycle
          active(jvert) = norm2(mesh%xyz(1:3,jvert) - mesh%xyz(1:3,ivert)) > r(1,jvert) + r(1,ivert)
       end do
    end do

    allocate(spheres(4,nseeds))
    do ivert = 1,nseeds
       spheres(1:3,ivert) = mesh%xyz(1:3,seeds(ivert))
       spheres(4,ivert) = r(1,seeds(ivert))
    end do
    
  end subroutine sample_spheres_brute_force













  subroutine edge_convexity( &
       brep, &
       iedge, &
       convexity )
    use mod_tolerances
    implicit none
    type(type_brep), intent(in)            :: brep
    integer,         intent(in)            :: iedge
    integer,         intent(out)           :: convexity
    type(type_intersection_curve), pointer :: curve => null()
    real(kind=fp)                          :: uv(2)
    real(kind=fp)                          :: normals(3,2), jacobian(3,2)
    real(kind=fp)                          :: ndotn, tangent(3), angle
    integer                                :: isurf, ivar, ipoint

    curve => brep%edges(iedge)%curve
    ipoint = curve%polyline%np/2
    !print *,'   IPOINT =',ipoint

    do isurf = 1,2
       uv = curve%polyline%uv(1:2,isurf,ipoint)
       do ivar = 1,2
          call evald1( &
               jacobian(1:3,ivar), &
               curve%surf(isurf)%ptr, &
               uv, &
               ivar )
       end do
       normals(1:3,isurf) = cross(jacobian(1:3,1), jacobian(1:3,2))
       normals(1:3,isurf) = normals(1:3,isurf)/norm2(normals(1:3,isurf))
    end do

    ndotn = dot_product(normals(1:3,1),normals(1:3,2))
    !print *,'   NDOTN =',ndotn
    
    if ( abs(ndotn) > 1._fp - 1.d-9 ) then
       convexity = 0 ! smooth
       return
    end if

    tangent = cross(normals(1:3,1),normals(1:3,2))
    tangent = tangent/norm2(tangent)

    angle = atan2(dot_product(normals(1:3,2), cross(normals(1:3,1),tangent)), ndotn)
    PRINT *,'   ANGLE =',angle

    if ( angle < 0._fp ) then
       convexity = 2 ! convex
    else
       convexity = 1 ! concave
    end if
    
  end subroutine edge_convexity


  
end program demo_EoS_brep
