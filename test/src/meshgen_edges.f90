program meshgen_edges

  use mod_math
  use mod_options
  use mod_util
  use mod_diffgeom
  use mod_polynomial
  use mod_types_intersection
  use mod_intersection
  use mod_types_brep
  use mod_brep
  use mod_init
  
  implicit none

  type(type_options)                            :: options_inter, options_mesh
  integer                                       :: fid
  character(3)                                  :: strnum3
  
  type(type_surface), allocatable, target       :: surf(:)
  integer                                       :: nsurf, isurf
  type(type_intersection_data), target          :: interdata
  integer                                       :: icurv
  integer, allocatable                          :: curvetype(:)
  integer                                       :: i

  type(type_brep)                               :: brep

  type(type_intersection_polyline), allocatable :: mgedges(:)
  

  !! read options file
  call read_options( &
       'meshgen_edges/intersections.opt', &
       options_inter )

  !! import surfaces
  call get_free_unit(fid)
  open( &
       unit = fid, &
       file = 'meshgen_edges/init/surftag.dat', &
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
          'meshgen_edges/init/coef/c_' // strnum3 // '.cheb', &
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
       'meshgen_edges/init/curves.dat', &
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
       options_inter%chord_err, &
       options_inter%hmin, &
       options_inter%hmax )

   ! --->>> DEBUG 
  call write_intersection_data( &
       interdata, &
       'meshgen_edges/brep/intersection_points.dat', &
       'meshgen_edges/brep/intersection_curves.dat' )
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
       'meshgen_edges/brep/verts.dat', &
       'meshgen_edges/brep/edges.dat', &
       'meshgen_edges/brep/faces.dat' )
  call write_brep_edges_geometry( &
       brep, &
       'meshgen_edges/brep/edges_xyz.dat', &
       'meshgen_edges/brep/edges_uv.dat' )
  ! <<<---



  !! read options file
  call read_options( &
       'meshgen_edges/mesh.opt', &
       options_mesh )

  allocate(mgedges(brep%ne))

  call make_all_meshgen_edges( &
       brep, &
       options_mesh, &
       mgedges )



  call write_mgedges( &
       mgedges(1:brep%ne), &
       brep%ne, &
       'meshgen_edges/mgedges/xyz.dat', &
       'meshgen_edges/mgedges/uv.dat' )







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











  subroutine write_mgedges( &
       mgedges, &
       ne, &
       filexyz, &
       fileuv )
    use mod_util
    implicit none
    integer,                          intent(in) :: ne
    type(type_intersection_polyline), intent(in) :: mgedges(ne)
    character(*),                     intent(in) :: filexyz, fileuv
    integer                                      :: iedge, ivert, fid, fid2

    call get_free_unit(fid)
    open(unit=fid, file=filexyz, action='write')
    call get_free_unit(fid2)
    open(unit=fid2, file=fileuv, action='write')
    write (fid,*) ne
    write (fid2,*) ne
    do iedge = 1,ne
       write (fid,*) mgedges(iedge)%np
       write (fid2,*) mgedges(iedge)%np
       do ivert = 1,mgedges(iedge)%np
          write (fid,*) mgedges(iedge)%xyz(1:3,ivert)
          write (fid2,*) mgedges(iedge)%uv(1:2,1:2,ivert)
       end do
    end do
    close(fid)
    close(fid2)

  end subroutine write_mgedges





  
end program meshgen_edges
