subroutine generate_brep_mesh( &
     brep, &
     hmin, &
     hmax, &
     tol )
  use mod_brep2
  use mod_mesh
  use mod_util
  implicit none
  type(type_brep),         intent(in)    :: brep
  real(kind=fp),           intent(in)    :: hmin
  real(kind=fp),           intent(in)    :: hmax
  real(kind=fp),           intent(in)    :: tol
  type(type_surface_mesh), intent(inout) :: mesh
  integer                                :: finf, fpts, fedg
  integer                                :: bv2mv(brep%nv), be2mv(2,brep%ne)
  integer, allocatable                   :: loc2glob(:), idsf(:), typf(:)
  real(kind=fp), allocatable             :: uvf(:,:), xyzf(:,:)
  integer, allocatable                   :: trif(:,:)
  integer                                :: npf, np, ntrif
  integer                                :: ifirstedge, ifirstpoint, sens, ifirst, ilast
  integer                                :: iface, iloop, ihedg(2), i

  ! write 'info' file (common to all faces)
  call get_free_unit(finf)
  open(unit=finf, file='tmp/info.dat', action='write')
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
     call write_polynomial(brep%faces(iface)%surface%x, 'tmp/c.cheb')
     !
     ! open boundary points & edges files
     call get_free_unit(fpts)
     open(unit=fpts, file='tmp/bpts.dat', action='write')
     call get_free_unit(fedg)
     open(unit=fedg, file='tmp/bedg.dat', action='write')
     !
     npf = 0
     loops : do iloop = 0,brep%faces(iface)%ninner ! <-----------------------------------+
        ! first halfedge of the loop                                                     !
        if ( iloop == 0 ) then ! <-------------------+                                   !
           ihedg = brep%faces(iface)%outer           !                                   !
        else ! --------------------------------------+                                   !
           ihedg = brep%faces(iface)%inner(:,iloop)  !                                   !
        end if ! <-----------------------------------+                                   !
        ifirstedge = ihedg(1)                                                            !
        !                                                                                !
        ifirstpoint = npf + 1                                                            !
        ! traverse the loop                                                              !
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
           if ( allocated(uvf) .and. size(uvf,2) < np ) deallocate(uv)
           if ( .not.allocated(uvf) ) allocate(uvf(2,np))
           uvf(1:2,1:np) = brep%edges(ihedg(1))%curve%polyline%uv&
                (:,sens,ifirst:ilast:(-1)**sens)
           !
           if ( allocated(xyzf) .and. size(xyzf,2) < np ) deallocate(xyzf)
           if ( .not.allocated(xyzf) ) allocate(xyzf(2,np))
           xyzf(1:2,1:np) = brep%edges(ihedg(1))%curve%polyline%xyz&
                (:,ifirst:ilast:(-1)**sens)
           !                                                                          !  !
           ! get brep vertex index of first point
           ivert = get_orig(brep, ihedg)
           if ( bv2mv(ivert) == 0 ) then ! <-----------------------------------+
              ! the brep vertex is not yet present in the mesh                 !
              bv2mv(ivert) = npf + 1                                           !
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
               call append_vertices( &                                         !
                   mesh, &                                                     !
                   xyzf(1:3,1), &                                              !
                   uvf(1:2,1), &                                               !
                   idsf(npf+1), &                                              !
                   typf(npf+1), &                                              !
                   1 )                                                         !
           end if ! <----------------------------------------------------------+
           loc2glob(npf+1) = bv2mv(ivert)
           !
           if ( be2mv(1,ihedg(1)) == 0 ) then ! <------------------------------+
              ! the brep edge is not yet present in the mesh                   !
              ! set first and last points                                      !
              if ( ihedg(1) == 1 ) then ! <--------+                           !
                 be2mv(1,ihedg(1)) = npf + 2       !                           !
                 be2mv(2,ihedg(1)) = npf + np - 1  !                           !
              else ! ------------------------------+                           !
                 be2mv(1,ihedg(1)) = npf + np - 1  !                           !
                 be2mv(2,ihedg(1)) = npf + 2       !                           !
              end if ! <---------------------------+                           !
              !                                                                !
              if ( feat_edge(ihedg(1)) ) then ! <-----------------------+      !
                 ! the mesh vertices are locked on a hyperedge          !      !
                 idsf(npf+1:npf+np-1) = ihedg(1)                        !      !
                 typf(npf+1:npf+np-1) = 1                               !      !
              else ! ---------------------------------------------------+      !
                 ! the mesh vertices are free to move along a hyperface !      !
                 idsf(npf+1:npf+np-1) = iface                           !      !
                 typf(npf+1:npf+np-1) = 2                               !      !
              end if ! <------------------------------------------------+      !
              !                                                                !
              ! append new vertices to the mesh                                !
              call append_vertices( &                                          !
                   mesh, &                                                     !
                   xyzf(1:3,2:np-1), &                                         !
                   uvf(1:2,2:np-1), &                                          !
                   idsf(npf+1), &                                              !
                   typf(npf+1), &                                              !
                   np-2 )                                                      !
           end if ! <----------------------------------------------------------+
           loc2glob(npf+2:npf+np-1) = [(i, i=be2mv(1,ihedg(1)),be2mv(2,ihedg(1)),(-1)**sens)]
           !
           ! write points...                            
           do i = 1,np-1 ! <---------------------------+
              write (fpts,*) uvf(:,i)                  !
           end do ! <----------------------------------+
           ! ...and edges
           do i = 1,np-2 ! <---------------------------+
              write (fedg,*) npf + i, npf + i + 1      !
           end do ! <----------------------------------+
           !
           ! move on to the next halfedge on the loop
           if ( ihedg(1) == ifirstedge ) then ! <------+
              ! the loop is complete                   !
              write (fedg,*) npf, ifirstpoint          !
              exit                                     !
           else ! -------------------------------------+
              write (fedg,*) npf, npf+1                !
           end if ! <----------------------------------+
           !           
        end do halfedges
     end do loops
     !
     ! close files
     close(fpts)
     close(fedg)
     !
     ! run mesher
     PRINT *,'MESH FACE #',IFACE
     call system('/stck/bandrieu/Bureau/MeshGen/./meshgen.out &
          & tmp/c.cheb &
          & tmp/bpts.dat &
          & tmp/bedg.dat &
          & tmp/info.dat &
          & tmp/tri.dat &
          & tmp/uv.dat &
          & tmp/xyz.dat' )
     !
     ! read face submesh                                                                     !
     call read_triangles( &                                                                  !
          'tmp/tri.dat', &                                                                   !
          trif, &                                                                            !
          ntrif )                                                                            !
     call read_points( &                                                                     !
          'tmp/uv.dat', &                                                                    !
          uvf, &                                                                             !
          2, &                                                                               !
          np )                                                                               !
     call read_points( &                                                                     !
          'tmp/xyz.dat', &                                                                   !
          xyzf, &                                                                            !
          3, &                                                                               !
          np )     
     !
     if ( size(loc2glob) < np ) call reallocate_list(loc2glob, np)
     loc2glob(npf+1:np) = mesh%nv + [(i, i=1,np-npf)]
     do i = 1,3 ! <----------------------------------+
        trif(i,1:ntrif) = loc2glob(trif(i,1:ntrif))  !
     end do ! <--------------------------------------+                                       !
     !                                                                                       !
     if ( size(idsf) < np ) call reallocate_list(idsf, np)                                   !
     idsf(npf+1:np) = iface                                                                  !
     !                                                                                       !
     if ( size(typf) < np ) call reallocate_list(typf, np)                                   !
     typf(npf+1:np) = 2                                                                      !
     !                                                                                       !
     ! add new triangles                                                                     !
     call append_triangles( &                                                                !
          mesh, &                                                                            !
          trif(1:3,1:ntrif), &                                                               !
          ntrif )                                                                            !
     !                                                                                       !
     ! add new mesh vertices                                                                 !
     call append_vertices( &                                                                 !
          mesh, &                                                                            !
          xyzf(1:3,npf+1:npf+np), &                                                          !
          uvf(1:2,npf+1:npf+np), &                                                           !
          idsf(npf+1:npf+np), &                                                              !
          typf(npf+1:npf+np), &                                                              !
          np-npf )
     !
  end do faces
end subroutine generate_brep_mesh
