program import_geometry

  use mod_util
  use mod_diffgeom
  use mod_import
  use mod_geometry
  use mod_types_brep
  use mod_types_intersection
  use mod_halfedge
  use mod_brep
  use mod_intersection
  use mod_hypergraph
  use mod_mesh
  use mod_geometry
  
  implicit none

  real(kind=fp), parameter                  :: TOLang = 10._fp ! degrees
  real(kind=fp), parameter                  :: TOLcosang = cos(CSTpi * TOLang / 180._fp)
  
  !type(type_surface), allocatable           :: surf(:)
  !integer                                   :: nsurf
  type(type_surface_mesh)                   :: mesh
  type(type_brep)                           :: brep
  type(type_intersection_data), target      :: interdata
  type(type_hypergraph)                     :: hypergraph

  character(100)                            :: filename
  integer                                   :: fid
  integer, allocatable                      :: ref2bf(:)
  integer, allocatable                      :: refs(:)
  integer                                   :: nrefs
  integer, allocatable                      :: ibrepverts(:)
  integer, allocatable                      :: bv2hf(:,:), nbv2hf(:)
  integer                                   :: ntmp
  real(kind=fp)                             :: uvpoint(2,0)
  type(ptr_surface)                         :: surfpoint(0)
  integer                                   :: idpoint
  logical, allocatable                      :: visited(:,:)
  logical                                   :: smooth
  real(kind=fp)                             :: nor(3,2)
  integer                                   :: np
  type(type_intersection_polyline), pointer :: polyline => null()
  logical, allocatable                      :: featurepath(:)
  integer                                   :: ivert, jvert, kvert, iface, jface, iref, jref
  integer, dimension(2)                     :: ihedg, jhedg, iendpoints

  if ( command_argument_count() < 1 ) then
     STOP 'you must provide an input mesh filename'
  else
     call get_command_argument(1, filename)
  end if
  
  
  ! Read xyz, tri and face ref
  call read_msh( &
       trim(filename), &
       mesh%xyz, &
       mesh%nv, &
       mesh%tri, &
       mesh%ihf, &
       mesh%nt )

  allocate(ref2bf(maxval(mesh%ihf(1:mesh%nt))))
  ref2bf(:) = -1
  brep%nf = 0
  do iface = 1,mesh%nt
     iref = mesh%ihf(iface)
     if ( ref2bf(iref) < 0 ) then
        brep%nf = brep%nf + 1
        ref2bf(iref) = brep%nf
     end if
  end do
  mesh%ihf(1:mesh%nt) = ref2bf(mesh%ihf(1:mesh%nt))
  deallocate(ref2bf)
  
  PRINT *,MESH%NV,' VERTICES,', MESH%NT, ' TRIANGLES'

  PRINT *,BREP%NF,' BREP FACES'
  ! read surfaces
  !nsurf = brep%nf
  !do isurf = 1,nsurf
  !     write (strnum3,'(i3.3)') isurf
  !     call read_polynomial( &
  !          surf(isurf)%x, &
  !          trim(dir) // 'coef/c_' // strnum3 // '.cheb', &
  !          nvar=2, &
  !          base=1 )
  !
  !     call economize2( &
  !          surf(isurf)%x, &
  !          EPSmath )
  !
  !     call compute_deriv1(surf(isurf))
  !     call compute_deriv2(surf(isurf))
  !     call compute_pseudonormal(surf(isurf))
  !     call economize2( &
  !          surf(isurf)%pn, &
  !          EPSmath )
  !  end do

  

  allocate(mesh%uv(2,2,mesh%nv), mesh%ids(mesh%nv), mesh%typ(mesh%nv))
  IF ( .FALSE. ) THEN
     call get_free_unit(fid)
     open( &
          unit = fid, &
          file = 'import_geometry/uv0.dat', &
          action = 'read' )
     do ivert = 1,mesh%nv
        read (fid,*) mesh%uv(1:2,1:2,ivert)
     end do
     close(fid)
  END IF
  mesh%typ(1:mesh%nv) = 2

  call write_inria_mesh( &
       mesh, &
       'import_geometry/inputmesh.mesh' )

  ! Make mesh halfedges
  call make_halfedges( &
       mesh )

  ! make BREP vertices
  allocate(refs(brep%nf))
  allocate(bv2hf(brep%nf,10))
  do ivert = 1,mesh%nv ! <------------------------------------+
     ihedg = mesh%v2h(:,ivert) ! outgoing                     !
     iface = get_face(ihedg)                                  !
     jface = iface                                            !
     nrefs = 0                                                !
     do ! <-----------------------------------------------+   !
        if ( jface < 1 ) then ! <---+                     !   !
           jref = -1                !                     !   !
        else ! ---------------------+                     !   !
           jref = mesh%ihf(jface)   !                     !   !
        end if ! -------------------+                     !   !
        do iref = 1,nrefs ! <-------------------------+   !   !
           if ( jref == refs(iref) ) exit             !   !   !
        end do ! <------------------------------------+   !   !
        !                                                 !   !
        if ( nrefs < 1 .or. iref > nrefs ) then ! <---+   !   !
           nrefs = nrefs + 1                          !   !   !
           refs(nrefs) = jref                         !   !   !
        end if ! <------------------------------------+   !   !
        !if ( nrefs > 2 ) exit                             !   !
        !                                                 !   !
        if ( jface < 1 ) exit                             !   !
        !                                                 !   !
        ihedg = get_prev(ihedg)       ! ingoing           !   !
        ihedg = get_twin(mesh, ihedg) ! outgoing          !   !
        jface = get_face(ihedg)                           !   !
        if ( jface == iface ) exit                        !   !
     end do ! <-------------------------------------------+   !
     !                                                        !
     if ( nrefs > 2 ) then ! <---------+                      !
        call insert_after( &           !                      !
             ibrepverts, &             !                      !
             brep%nv, &                !                      !
             ivert, &                  !                      !
             brep%nv )                 !                      !
        mesh%typ(ivert) = 0            !                      !
        mesh%ids(ivert) = brep%nv      !                      !
        !                              !                      !
        call add_intersection_point( & !                      !
             uvpoint, &                !                      !
             mesh%xyz(1:3,ivert), &    !                      !
             surfpoint, &              !                      !
             0, &                      !                      !
             interdata, &              !                      !
             idpoint )                 !                      !
        !                              !                      !
        ntmp = brep%nv - 1             !                      !
        call insert_column_after( &    !                      !
             bv2hf, &                  !                      !
             brep%nf, &                !                      !
             ntmp, &                   !                      !
             refs, &                   !                      !
             ntmp )                    !                      !
        ntmp = brep%nv - 1             !                      !
        call insert_after( &           !                      !
             nbv2hf, &                 !                      !
             ntmp, &                   !                      !
             nrefs, &                  !                      !
             ntmp )                    !                      !
     end if ! <------------------------+                      !
  end do ! <--------------------------------------------------+
  deallocate(refs)

  PRINT *,BREP%NV,' BREP VERTICES'

  ! make paths/brep (half)edges/intersection curves
  allocate(brep%verts(brep%nv))
  allocate(visited(3,mesh%nt))
  visited(1:3,1:mesh%nt) = .false.
  do jvert = 1,brep%nv ! <-----------------------------------------------------------+
     brep%verts(jvert)%point => interdata%points(jvert)                              !
     ivert = ibrepverts(jvert)                                                       !
     ihedg = mesh%v2h(:,ivert)                                                       !
     iface = get_face(ihedg)                                                         !
     jface = iface                                                                   !
     do ! <----------------------------------------------------------------------+   !
        iref = mesh%ihf(jface)                                                   !   !
        jhedg = get_twin(mesh, ihedg)                                            !   !
        jface = get_face(jhedg)                                                  !   !
        if ( jface < 1 ) then ! <----+                                           !   !
           jref = iref - 1           !                                           !   !
        else ! ----------------------+                                           !   !
           jref = mesh%ihf(jface)    !                                           !   !
        end if ! <-------------------+                                           !   !
        !                                                                        !   !
        if ( .not.visited(ihedg(1),ihedg(2)) .and. iref /= jref ) then ! <---+   !   !
           if ( jface < 1 ) then ! <---------------------------------+       !   !   !
              smooth = .false.                                       !       !   !   !
           else ! ---------------------------------------------------+       !   !   !
              nor(:,1) = triangle_normal( &                          !       !   !   !
                   mesh%xyz(:,mesh%tri(:,iface)), &                  !       !   !   !
                   .true. )                                          !       !   !   !
              nor(:,2) = triangle_normal( &                          !       !   !   !
                   mesh%xyz(:,mesh%tri(:,jface)), &                  !       !   !   !
                   .true. )                                          !       !   !   !
              smooth = dot_product(nor(:,1), nor(:,2)) > TOLcosang   !       !   !   !
           end if ! <------------------------------------------------+       !   !   !
           ! make new path                                                   !   !   !
           call make_feature_path( &                                         !   !   !
                mesh, &                                                      !   !   !
                visited, &                                                   !   !   !
                ihedg )                                                      !   !   !
           mesh%paths(mesh%npaths)%hyperedge = mesh%npaths                   !   !   !
           !                                                                 !   !   !
           np = mesh%paths(mesh%npaths)%nv                                   !   !   !
           mesh%typ(mesh%paths(mesh%npaths)%verts(2:np-1)) = 1               !   !   !
           mesh%ids(mesh%paths(mesh%npaths)%verts(2:np-1)) = mesh%npaths     !   !   !
           kvert = mesh%paths(mesh%npaths)%verts(np)                         !   !   !
           if ( iref > jref ) then ! <-----+                                 !   !   !
              iendpoints = [ivert,kvert]   !                                 !   !   !
           else ! -------------------------+                                 !   !   !
              iendpoints = [kvert,ivert]   !                                 !   !   !
           end if ! <----------------------+                                 !   !   !
           call add_intersection_curve( &                                    !   !   !
                interdata, &                                                 !   !   !
                [0._fp, 0._fp, 0._fp], &                                     !   !   !
                iendpoints, &                                                !   !   !
                spread(spread([-1._fp, 1._fp], 2, 2), 3, 2) )                !   !   !
           allocate(interdata%curves(interdata%nc)%polyline)                 !   !   !
           polyline => interdata%curves(interdata%nc)%polyline               !   !   !
           polyline%np = np                                                  !   !   !
           allocate(polyline%xyz(3,np))                                      !   !   !
           if ( iref > jref ) then ! <-----------------------------------+   !   !   !
              polyline%xyz(1:3,1:np) = &                                 !   !   !   !
                   mesh%xyz(1:3,mesh%paths(mesh%npaths)%verts(1:np))     !   !   !   !
           else ! -------------------------------------------------------+   !   !   !
              polyline%xyz(1:3,1:np) = &                                 !   !   !   !
                   mesh%xyz(1:3,mesh%paths(mesh%npaths)%verts(np:1:-1))  !   !   !   !
           end if ! <----------------------------------------------------+   !   !   !
           interdata%curves(interdata%nc)%isplit(2,1:2) = [1,np]             !   !   !
           interdata%curves(interdata%nc)%smooth = smooth                    !   !   !
           !                                                                 !   !   !
        end if ! <-----------------------------------------------------------+   !   !
        !                                                                        !   !
        ! cyle halfedges                                                         !   !
        ihedg = get_prev(ihedg)       ! ingoing                                  !   !
        ihedg = get_twin(mesh, ihedg) ! outgoing                                 !   !
        jface = get_face(ihedg)                                                  !   !
        if ( jface < 1 .or. jface == iface ) exit                                !   !
     end do ! <------------------------------------------------------------------+   !
  end do ! <-------------------------------------------------------------------------+

  PRINT *,MESH%NPATHS,' FEATURE PATHS'


  ! make wires
  allocate(brep%faces(brep%nf))
  hypergraph%nhf = brep%nf
  allocate(hypergraph%hyperfaces(hypergraph%nhf))
  !nwrires = min(brep%nv, mesh%npaths)
  !allocate(wire2path(min
  do iface = 1,brep%nf
     call insert_after( &
          hypergraph%hyperfaces(iface)%faces, &
          hypergraph%hyperfaces(iface)%nf, &
          iface, &
          hypergraph%hyperfaces(iface)%nf )
     brep%faces(iface)%hyperface = iface
  end do


  ! remove false paths/hyperedges
  !allocate(featurepath(mesh%npaths))
  !do ipath = 1,mesh%npaths
  !   featurepath(ipath) = interdata%curves(ipath)%smooth
  !end do
  

 
  

  call write_mesh_files( &
       mesh, &
       'import_geometry/tri.dat', &
       'import_geometry/xyz.dat', &
       'import_geometry/uv.dat', &
       'import_geometry/idstyp.dat', &
       'import_geometry/paths.dat' )

contains


  subroutine make_feature_path( &
       mesh, &
       visited, &
       first_hedg )
    use mod_mesh
    implicit none
    type(type_surface_mesh), intent(inout) :: mesh
    logical,                 intent(inout) :: visited(3,mesh%nt)
    integer,                 intent(in)    :: first_hedg(2)
    type(type_path)                        :: path
    integer, dimension(2)                  :: ihedg, jhedg
    integer                                :: iface, jface, iref, jref, orig, dest
    
    ihedg = first_hedg
    do ! <-------------------------------------------------------+
       iface = get_face(ihedg)                                   !
       iref = mesh%ihf(iface)                                    !
       jhedg = get_twin(mesh, ihedg)                             !
       jface = get_face(jhedg)                                   !
       if ( jface < 1 ) then ! <---+                             !
          jref = iref - 1          !                             !
       else ! ---------------------+                             !
          jref = mesh%ihf(jface)   !                             !
       end if ! <------------------+                             !
       !                                                         !
       if ( .not.visited(ihedg(1),ihedg(2)) .and. &              !
            iref /= jref ) then ! <--------------------------+   !
          visited(ihedg(1),ihedg(2)) = .true.                !   !
          if ( jface > 0 ) visited(jhedg(1),jface) = .true.  !   !
          !                                                  !   !
          orig = get_orig(mesh, ihedg)                       !   !
          call insert_after( &                               !   !
               path%verts, &                                 !   !
               path%nv, &                                    !   !
               orig, &                                       !   !
               path%nv )                                     !   !
          dest = get_dest(mesh, ihedg)                       !   !
          if ( mesh%typ(dest) == 0 ) then ! <---+            !   !
             call insert_after( &               !            !   !
                  path%verts, &                 !            !   !
                  path%nv, &                    !            !   !
                  dest, &                       !            !   !
                  path%nv )                     !            !   !
             call append_path( &                !            !   !
                  mesh, &                       !            !   !
                  path )                        !            !   !
             return                             !            !   !
          end if ! <----------------------------+            !   !
          !                                                  !   !
          ihedg = get_next(ihedg) ! ingoing                  !   !
          cycle                                              !   !
       end if ! <--------------------------------------------+   !
       !                                                         !
       ! cyle halfedges around current vertex                    !
       ihedg = get_twin(mesh, ihedg) ! ingoing                   !
       ihedg = get_next(ihedg)       ! outgoing                  !
    end do ! <---------------------------------------------------+
    
  end subroutine make_feature_path
  
end program import_geometry
