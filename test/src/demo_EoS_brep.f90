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
  integer                                 :: fid
  character(3)                            :: strnum3

  type(type_surface), allocatable, target :: surf(:)
  integer                                 :: nsurf, isurf
  type(type_intersection_data), target    :: interdata
  integer                                 :: icurv
  integer, allocatable                    :: curvetype(:), curvetypetmp(:)

  type(type_brep)                         :: brep
  integer                                 :: iface, jface, ivert
  type(type_hypergraph)                   :: hypergraph
  logical, allocatable                    :: feat_edge(:), feat_vert(:)
  type(type_surface_mesh)                 :: mesh

  real(kind=fp)                           :: timestep
  type(type_surface), allocatable, target :: surf_new(:)
  integer                                 :: nsurf_new
  type(type_intersection_data), target    :: interdata_new
  type(type_brep)                         :: brep_new

  integer                                 :: parameterization_mode = -1


  real(kind=fp)                           :: xyzverif(3,2)
  integer                                 :: i
  !-------------------------------------------------

  !! Read arguments
  if ( command_argument_count() < 1 ) then
     MAKE_EOS = .false.
  else
     call get_command_argument(1, argstr)
     read (argstr,*) argnum
     MAKE_EOS = (argnum > 0)

     if ( command_argument_count() < 2 ) then
          timestep = 1._fp
     else
          call get_command_argument(2, argstr)
          read (argstr,*) timestep
     end if
  end if

  !! read options file
  call read_options( &
       'demo_EoS_brep/demo_brep.opt', &
       options )

  !! import surfaces
  call get_free_unit(fid)
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
  END IF
  
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

  IF ( .NOT.MAKE_EOS ) THEN
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

     call write_mesh_files( &
          mesh, &
          'demo_EoS_brep/mesh/tri.dat', &
          'demo_EoS_brep/mesh/xyz.dat', &
          'demo_EoS_brep/mesh/uv.dat', &
          'demo_EoS_brep/mesh/idstyp.dat', &
          'demo_EoS_brep/mesh/paths.dat' )
     call write_halfedges( &
          mesh, &
          'demo_EoS_brep/mesh/mv2h.dat', &
          'demo_EoS_brep/mesh/mtwin.dat' )

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

     STOP
  END IF







  ! ************ MAKE EoB ************

  !timestep = 1._fp
  ! de-economize
  do isurf = 1,nsurf
     surf(isurf)%x%degr(1) = size(surf(isurf)%x%coef,1) - 1
     surf(isurf)%x%degr(2) = size(surf(isurf)%x%coef,2) - 1 
     call compute_deriv1(surf(isurf))
     call compute_deriv2(surf(isurf))
     call compute_pseudonormal(surf(isurf))
  end do


  call make_partial_eos( &
       parameterization_mode, &
       surf(1:nsurf), &
       nsurf, &
       interdata, &
       curvetype, &
       brep, &
       timestep, &
       surf_new, &
       nsurf_new, &
       interdata_new )
  
  call check_intersection_points( &
       interdata_new )
  call check_intersection_polylines( &
       interdata_new )

  do isurf = 1,nsurf_new
     call economize2( &
          surf_new(isurf)%x, &
          EPSmath )
     call compute_deriv1(surf_new(isurf))
     call compute_deriv2(surf_new(isurf))
     call compute_pseudonormal(surf_new(isurf))
  end do


  call intersect_all_surfaces( &
       surf_new(1:nsurf_new), &
       nsurf_new, &
       interdata_new, &
       [((i == 2 .or. i == 8 .or. i == 9), i=1,nsurf_new)], &
       options%chord_err, &
       options%hmin, &
       options%hmax )

  ! --->>> DEBUG 
  call write_intersection_data( &
       interdata_new, &
       'demo_EoS_brep/debug/intersection_points_new.dat', &
       'demo_EoS_brep/debug/intersection_curves_new.dat' )
  ! <<<---

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
  call write_brep_edges_geometry( &
       brep_new, &
       'demo_EoS_brep/debug/edges_xyz_new.dat', &
       'demo_EoS_brep/debug/edges_uv_new.dat' )
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

  call write_inria_mesh( &
       mesh, &
       'demo_EoS_brep/mesh_eos/mesh_eos_optim.mesh' )
  call write_obj_mesh( &
       mesh, &
       'demo_EoS_brep/mesh_eos/mesh_eos_optim.obj' )
  call write_vtk_mesh( &
       mesh, &
       'demo_EoS_brep/mesh_eos/mesh_eos_optim.vtk' )


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

       ! ***********
       interdata%curves(interdata%nc)%param_vector = xyz(1:3,np) - xyz(1:3,1)
       interdata%curves(interdata%nc)%param_vector = interdata%curves(interdata%nc)%param_vector/&
            norm2(interdata%curves(interdata%nc)%param_vector)
       ! ***********

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
    use mod_intersection
    use mod_projection
    use mod_tolerances
    implicit none
    type(type_intersection_data), intent(in) :: interdata
    real(kind=fp)                            :: xyz(3,2), err(2)
    real(kind=fp), dimension(4)              :: upperb, lowerb
    real(kind=fp)                            :: uvc(2,2), xyzc(3)
    integer                                  :: stat
    integer                                  :: ic, ip, is

    upperb(1:4) = 1._fp + EPSuv
    lowerb(1:4) = -upperb(1:4)
    
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
          if ( any(err > 100._fp*EPSxyz) ) then
             PRINT *, 'CURVE #',ic, ip, err
             cycle
             !   do is = 1,2
             !      call eval( &
             !           xyz(1:3,is), &
             !           interdata%curves(ic)%surf(1+mod(is,2))%ptr, &
             !           interdata%curves(ic)%polyline%uv(1:2,is,ip) )
             !      err(is) = norm2(xyz(1:3,is) - &
             !           interdata%curves(ic)%polyline%xyz(1:3,ip))
             !   end do
             !   PRINT *, err
             do is = 1,2
                uvc(1:2,is) = interdata%curves(ic)%polyline%uv(1:2,is,ip)
             end do
             IF ( .TRUE. ) THEN
                call projection_surface( &
                     interdata%curves(ic)%surf(2)%ptr, &
                     interdata%curves(ic)%polyline%xyz(1:3,ip), &
                     uvc(1:2,2), &
                     lowerb, &
                     upperb, &
                     stat, &
                     xyzc )
             ELSE
                call simultaneous_point_inversions( &
                     interdata%curves(ic)%surf, &
                     lowerb, &
                     upperb, &
                     stat, &
                     uvc, &
                     xyzc )
             END IF
          
             if ( stat == 0 ) then
                do is = 1,2
                   call eval( &
                        xyz(1:3,is), &
                        interdata%curves(ic)%surf(is)%ptr, &
                        interdata%curves(ic)%polyline%uv(1:2,is,ip) )
                end do
                PRINT *,'CONVERGE, ERR=',NORM2(XYZ(1:3,1) - XYZ(1:3,2))
                interdata%curves(ic)%polyline%xyz(1:3,ip) = 0.5_FP*SUM(XYZ,2)
             else
                PRINT *,'NON CONVERGE'
             end if
          end if
       end do
    end do

  end subroutine check_intersection_polylines

  

end program demo_EoS_brep
