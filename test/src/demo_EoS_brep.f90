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
  
  !-------------------------------------------------
  implicit none
  type(type_options)                      :: options
  type(type_surface), allocatable, target :: surf(:)
  integer                                 :: nsurf, isurf
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
  !-------------------------------------------------

  call get_free_unit(fid)

  
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
  ! <<<---


  !=============================================
  allocate(edgetype(brep%ne))
  do iedge = 1,brep%ne
     do icurv = 1,interdata%nc
        if ( associated(brep%edges(iedge)%curve, interdata%curves(icurv)) ) then
           edgetype(iedge) = curvetype(icurv)
        end if
     end do
  end do
  
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
  !=============================================




  

  do iface = 1,brep%nf
     PRINT *,'brepmesh face#',iface
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
           mesh%uv(1:2,1:2,ivert) = mesh%uv(1:2,[2,1],ivert) 
        END IF
     end if
  end do
  ! -----<<<
  
  call check_uvs(brep, mesh)
  
  call write_connectivity( &
       'demo_EoS_brep/mesh/', &
       mesh, &
       0 )
  call write_xyz_positions( &
       'demo_EoS_brep/mesh/', &
       mesh, &
       0 )

  !STOP

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
