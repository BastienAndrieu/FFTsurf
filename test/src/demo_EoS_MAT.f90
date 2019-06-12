program demo_EoS_MAT

  use mod_options
  use mod_util
  use mod_polynomial
  use mod_diffgeom
  use mod_types_intersection
  use mod_types_brep
  use mod_tolerances
  use mod_eos
  use mod_intersection
  use mod_brep
  use mod_hypergraph
  use mod_mesh
  use mod_init
  use mod_optimmesh

  implicit none

  type(type_options)                      :: options

  character(10)                           :: argstr
  real(kind=fp)                           :: timestep

  integer                                 :: parameterization_mode = 2

  type(type_surface), target              :: surf(2)
  integer                                 :: npb
  type(type_intersection_data), target    :: interdata
  type(type_brep)                         :: brep
  type(type_hypergraph)                   :: hypergraph
  logical, allocatable                    :: feat_edge(:), feat_vert(:)

  integer                                 :: iborder, ipoint, jpoint
  type(type_intersection_curve), pointer  :: curve => null()
  integer                                 :: curvetype(4)
  integer                                 :: idendpoint(2)
  real(kind=fp)                           :: uvendpoint(2,2)
  type(ptr_surface)                       :: surfpair(2)
  real(kind=fp), allocatable              :: uvxyzb(:,:)

  type(type_surface), allocatable, target :: surf_new(:)
  integer                                 :: nsurf_new
  type(type_intersection_data), target    :: interdata_new
  type(type_brep)                         :: brep_new
  type(type_hypergraph)                   :: hypergraph_new

  type(type_surface_mesh)                 :: mesh

  integer                                 :: isurf, iface, j

  options%chord_err = 5.d-4
  options%hmin      = 1.d-3
  options%hmax      = 1.d-2


  !! Read argument
  if ( command_argument_count() < 1 ) then
     timestep = 1._fp
  else
     call get_command_argument(1, argstr)
     read (argstr,*) timestep
  end if

  !! read skeleton surface parameterization
  do isurf = 1,2
     call read_polynomial( &
          surf(isurf)%x, &
          'demo_EoS_MAT/c_skeleton.cheb', &
          nvar=2, &
          base=1 )

     if ( isurf == 2 ) then
        do j = 1,size(surf(isurf)%x%coef,2)
           if ( mod(j,2) == 0 ) then
              surf(isurf)%x%coef(:,j,1:3) = -surf(isurf)%x%coef(:,j,1:3)
           end if
        end do
     end if

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

  do isurf = 1,2
     surfpair(isurf)%ptr => surf(1+mod(isurf,2))
  end do


  !! border --> interdata
  npb = 100
  allocate(uvxyzb(5,npb))
  do iborder = 1,4
     call trace_border_polyline( &
          surf(1), &
          iborder, &
          npb, &
          uvxyzb(1:2,1:npb), &
          uvxyzb(3:5,1:npb) )

     ! add endpoints
     do jpoint = 1,2
        ipoint = 1 + (jpoint-1)*(npb-1)
        uvendpoint(1:2,2) = uvxyzb(1:2,ipoint)
        uvendpoint(1,1) = uvxyzb(1,ipoint)
        uvendpoint(2,1) = -uvxyzb(2,ipoint)
        call add_intersection_point( &
             uvendpoint, &
             uvxyzb(3:5,ipoint), &
             surfpair, &
             2, &
             interdata, &
             idendpoint(jpoint) )
     end do

     ! add curve
     call add_intersection_curve( &
          interdata, &
          [0._fp, 0._fp, 0._fp], &
          idendpoint, &
          spread(spread([-1._fp-EPSuv, 1._fp+EPSuv], 2, 2), 3, 2) )
     curve => interdata%curves(interdata%nc)
     curvetype(interdata%nc) = 2
     curve%smooth = .false.
     curve%param_vector = uvxyzb(3:5,npb) - uvxyzb(3:5,1)
     curve%param_vector = curve%param_vector/norm2(curve%param_vector)
     curve%surf = surfpair
     curve%isplit(2,1:2) = [1,npb]
     allocate(curve%iedge(1))
     curve%iedge(:) = 0
     allocate(curve%polyline)
     curve%polyline%np = npb
     allocate(curve%polyline%uv(2,2,npb), curve%polyline%xyz(3,npb))
     curve%polyline%uv(1:2,2,1:npb) = uvxyzb(1:2,1:npb)
     curve%polyline%uv(  1,1,1:npb) = uvxyzb(1,1:npb)
     curve%polyline%uv(  2,1,1:npb) = -uvxyzb(2,1:npb)
     curve%polyline%xyz(1:3,1:npb) = uvxyzb(3:5,1:npb)
  end do


  ! --->>> DEBUG 
  call write_intersection_data( &
       interdata, &
       'demo_EoS_MAT/intersection_points.dat', &
       'demo_EoS_MAT/intersection_curves.dat' )
  ! <<<---

  !! make brep
  brep%nv = 0
  brep%ne = 0
  brep%nf = 0
  call make_brep_from_intersection_data( &
       surf, &
       2, &
       interdata, &
       brep )


  ! --->>> DEBUG
  call write_brep_files( &
       brep, &
       'demo_EoS_MAT/verts.dat', &
       'demo_EoS_MAT/edges.dat', &
       'demo_EoS_MAT/faces.dat' )
  call write_brep_edges_geometry( &
       brep, &
       'demo_EoS_MAT/edges_xyz.dat', &
       'demo_EoS_MAT/edges_uv.dat' )
  ! <<<---

  !! make hypergraph
  call make_hypergraph( &
       brep, &
       hypergraph, &
       feat_edge, &
       feat_vert )
  call write_hypergraph( &
       hypergraph, &
       'demo_EoS_MAT/hyperfaces.dat', &
       'demo_EoS_MAT/hyperedges.dat' )


  !!...
  call make_partial_eos( &
       parameterization_mode, &
       surf, &
       2, &
       interdata, &
       curvetype, &
       brep, &
       timestep, &
       surf_new, &
       nsurf_new, &
       interdata_new )


  !! ...
  do isurf = 1,nsurf_new
     call compute_pseudonormal(surf_new(isurf))
     call economize2( &
          surf_new(isurf)%pn, &
          EPSmath )
  end do

  !! make new brep
  brep_new%nv = 0
  brep_new%ne = 0
  brep_new%nf = 0
  call make_brep_from_intersection_data( &
       surf_new, &
       nsurf_new, &
       interdata_new, &
       brep_new )


  call write_brep_files( &
       brep_new, &
       'demo_EoS_MAT/verts_eos.dat', &
       'demo_EoS_MAT/edges_eos.dat', &
       'demo_EoS_MAT/faces_eos.dat' )
  call write_brep_edges_geometry( &
       brep, &
       'demo_EoS_MAT/edges_eos_xyz.dat', &
       'demo_EoS_MAT/edges_eos_uv.dat' )


  if ( allocated(feat_edge) ) deallocate(feat_edge)
  if ( allocated(feat_vert) ) deallocate(feat_vert)

  !! make hypergraph
  call make_hypergraph( &
       brep_new, &
       hypergraph_new, &
       feat_edge, &
       feat_vert )


  !! make mesh
  call generate_brep_conforming_mesh( &
       brep_new, &
       options%hmin, &
       options%hmax, &
       options%chord_err, &
       feat_edge(1:brep_new%ne), &
       feat_vert(1:brep_new%nv), &
       hypergraph_new%hyperedges(1:hypergraph_new%nhe), &
       hypergraph_new%nhe, &
       mesh )
  PRINT *, MESH%NV, ' VERTICES,', MESH%NT, ' TRIANGLES'
  call make_halfedges( &
          mesh )

  call write_inria_mesh( &
       mesh, &
       'demo_EoS_MAT/mesh_eos.mesh' )
  call write_obj_mesh( &
       mesh, &
       'demo_EoS_MAT/mesh_eos.obj' )
  call write_vtk_mesh( &
       mesh, &
       'demo_EoS_MAT/mesh_eos.vtk' )
  
  do iface = 1,mesh%nt
     mesh%ihf(iface) = brep_new%faces(mesh%ihf(iface))%hyperface
  end do

  call optim_jiao( &
       brep_new, &
       hypergraph_new%hyperedges(1:hypergraph_new%nhe), &
       hypergraph_new%nhe, &
       mesh, &
       1._fp, &!PARAM_frac_conf1, &
       0.7_fp, &!PARAM_frac_conf2, &
       0, &!PARAM_ipass1, &
       0, &!PARAM_ipass2, &
       20, &
       options%hmin, &
       options%hmax )

  call write_inria_mesh( &
       mesh, &
       'demo_EoS_MAT/mesh_eos_optim.mesh' )
  call write_obj_mesh( &
       mesh, &
       'demo_EoS_MAT/mesh_eos_optim.obj' )
  call write_vtk_mesh( &
       mesh, &
       'demo_EoS_MAT/mesh_eos_optim.vtk' )

end program demo_EoS_MAT
