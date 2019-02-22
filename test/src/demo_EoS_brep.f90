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
  type(type_options)                      :: options
  type(type_surface), allocatable, target :: surf(:)
  integer                                 :: nsurf, isurf, jsurf
  type(type_intersection_data), target    :: interdata
  integer                                 :: icurv
  integer, allocatable                    :: curvetype(:), edgetype(:)
  type(type_brep)                         :: brep
  integer                                 :: iface, jface, ivert, iedge
  type(type_hypergraph)                   :: hypergraph
  logical, allocatable                    :: feat_edge(:), feat_vert(:)
  type(type_surface_mesh)                 :: mesh

  integer                                 :: fid
  character(3)                            :: strnum3

  integer                                 :: ihedg(2)

  real(kind=fp)                           :: xyzverif(3,2)

  type(type_surface), allocatable, target :: surf_eos(:), edge_eos(:)
  integer                                 :: m, n, i, head, tail, sens, np
  real(kind=fp), allocatable              :: xyz(:,:,:), speed(:,:), s2e(:,:,:)
  type(ptr_surface)                       :: enve_rl(2)

  real(kind=fp), allocatable              :: xyzll(:,:), uvll(:,:)
  integer                                 :: stat
  type(type_surface)                      :: surfll
  real(kind=fp)                           :: longlat_range(2), solidangle
  !-------------------------------------------------

  call get_free_unit(fid)




  ! >>>>>-----------------
  IF ( .TRUE. )  THEN
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
     !call economize2( &
     !     surf(isurf)%x, &
     !     EPSmath )

     call compute_deriv1(surf(isurf))
     call compute_deriv2(surf(isurf))
     call compute_pseudonormal(surf(isurf))
     !call economize2( &
     !     surf(isurf)%pn, &
     !     EPSmath )
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
  open(unit=fid, file='demo_EoS_brep/debug/edges_xyz.dat', action='write')
  write (fid,*) brep%ne
  do iedge = 1,brep%ne
     call get_polyline_endpoints( &
          brep, &
          [iedge,1], &
          head, &
          tail, &
          sens, &
          np )
     write (fid,*) np
     do ivert = head,tail,(-1)**sens
        write (fid,*) brep%edges(iedge)%curve%polyline%xyz(1:3,ivert)
     end do
  end do
  close(fid)
  ! <<<---



  ! >>> ---------- Mark edges as smooth, concave or convex
  allocate(edgetype(brep%ne))
  do iedge = 1,brep%ne
     do icurv = 1,interdata%nc
        if ( associated(brep%edges(iedge)%curve, interdata%curves(icurv)) ) then
           edgetype(iedge) = curvetype(icurv)
        end if
     end do
  end do
  ! ----------<<<
  





  IF ( .false. ) THEN
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


  IF ( .FALSE. ) THEN
     !! make hypergraph
     call make_hypergraph( &
          brep, &
          hypergraph, &
          feat_edge, &
          feat_vert )

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

     do iface = 1,mesh%nt
        mesh%ihf(iface) = brep%faces(mesh%ihf(iface))%hyperface
     end do


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
     STOP
  END IF

  ! >>> ---------- EoS from faces 
  allocate(surf_eos(nsurf))
  do isurf = 1,nsurf
     write (strnum3,'(i3.3)') isurf
     
     m = size(surf(isurf)%x%coef,1)
     n = size(surf(isurf)%x%coef,2)
     !print *, m, n, size(surf(isurf)%x%coef,3)
     
     allocate(xyz(m,n,3), speed(m,n), s2e(m,n,3))
     call ifcht2( &
         surf(isurf)%x%coef(1:m,1:n,1:3), &
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
     
     call eos_surface( &
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

     call reset_polynomial(surf_eos(isurf)%x, 2, 1, [m-1,n-1], 3)
     call fcht2( &
          xyz, &
          surf_eos(isurf)%x%coef(1:m,1:n,1:3), &
          m, &
          n, &
          3, &
          EPSmath )
     
     deallocate(xyz, speed, s2e)

     call economize2( &
          surf_eos(isurf)%x, &
          EPSmath )
     call write_polynomial( &
          surf_eos(isurf)%x, &
          'demo_EoS_brep/debug/eos_c_'//strnum3 //'.cheb' )
  end do
  ! ----------<<<

  

  ! >>> ---------- EoS from convex edges
  allocate(edge_eos(brep%ne))
  do iedge = 1,brep%ne
     if ( edgetype(iedge) /= 2 ) cycle
     PRINT *,'IEDGE =',IEDGE
     write (strnum3,'(i3.3)') iedge

     call get_polyline_endpoints( &
         brep, &
         [iedge,1], &
         head, &
         tail, &
         sens, &
         np )

     do jsurf = 1,2
        do isurf = 1,nsurf
           if ( associated(brep%edges(iedge)%curve%surf(jsurf)%ptr, surf(isurf)) ) then
              enve_rl(jsurf)%ptr => surf_eos(isurf)
              exit
           end if
        end do
     end do
     
     call eos_from_curve( &
          brep%edges(iedge)%curve%polyline, &
          head, &
          tail, &
          enve_rl, &
          xyz, &
          m, &
          n )
     PRINT *,'M, N =', M, N

     call reset_polynomial(edge_eos(iedge)%x, 2, 1, [m-1,n-1], 3)
     call fcht2( &
          xyz, &
          edge_eos(iedge)%x%coef(1:m,1:n,1:3), &
          m, &
          n, &
          3, &
          EPSmath )

     call economize2( &
          edge_eos(iedge)%x, &
          EPSmath )
     call write_polynomial( &
          edge_eos(iedge)%x, &
          'demo_EoS_brep/debug/eos_edge_c_'//strnum3//'.cheb' )
     
     deallocate(xyz)
  end do
  ! ----------<<<



  

  ! >>> ---------- EoS from convex vertices
  do ivert = 1,brep%nv
     PRINT *,'VERTEX #',ivert
     ihedg = brep%verts(ivert)%halfedge
     iface = get_face(brep, ihedg)
     do
        PRINT *,'   ->',get_dest(brep, ihedg), ' :', edgetype(ihedg(1))
        ! traverse halfedges counter-clockwise
        ihedg = get_prev(brep, ihedg) ! previous halfedge
        ihedg = get_twin(ihedg) ! twin halfedge
        jface = get_face(brep, ihedg)

        if ( jface == iface ) exit
     end do
     PRINT *,''
  end do
  ! ----------<<<
  

  
  




  
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

  


end program demo_EoS_brep
