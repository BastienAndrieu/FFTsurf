module mod_optimmesh

  use mod_math
  
  implicit none

contains


  subroutine optim_jiao( &
       brep, &
       mesh, &
       passmax, &
       frac_conf, &
       hmin, &
       hmax )
    use mod_diffgeom
    use mod_intersection
    use mod_brep2
    use mod_mesh
    use mod_linalg
    implicit none
    type(type_brep),         intent(in)    :: brep
    type(type_surface_mesh), intent(inout) :: mesh
    integer,                 intent(in)    :: passmax
    real(kind=fp),           intent(in)    :: frac_conf
    real(kind=fp),           intent(in)    :: hmin
    real(kind=fp),           intent(in)    :: hmax
    real(kind=fp)                          :: hve(mesh%nv)
    real(kind=fp)                          :: wei(mesh%nt)
    real(kind=fp)                          :: ener(mesh%nt)
    real(kind=fp)                          :: grad(3,mesh%nv)
    real(kind=fp)                          :: hess(3,3,mesh%nv)
    real(kind=fp)                          :: duv(2,2,mesh%nv), ds
    real(kind=fp)                          :: tng(3,2), duv_ds(2,2,2)
    integer                                :: stat
    logical                                :: singular
    real(kind=fp), dimension(4)            :: lowerb, upperb
    real(kind=fp)                          :: uvtmp(2,2), xyztmp(3)
    integer                                :: ipass, ivert, ivar

    lowerb(:) = -2._fp
    upperb(:) = 2._fp
    
    do ipass = 1,passmax
       ! compute target edge lengths at vertices
       hve(:) = 1._fp
       ! ...
       ! compute triangle weights
       wei(:) = 1._fp
       ! ...
       ! compute energy, gradient and hessian
       call energy_jiao( &
            mesh%nt, &
            mesh%tri(1:3,1:mesh%nt), &
            wei, &
            mesh%nv, &
            mesh%xyz(1:3,1:mesh%nv), &
            hve, &
            frac_conf, &
            ener, &
            grad, &
            hess )
       !
       ! compute vertex displacements
       do ivert = 1,mesh%nv
          select case ( mesh%typ(ivert) ) ! <-------------------------------+
          case (0) ! -------------------------------------------------------+
             cycle                                                          !
          case (1) ! -------------------------------------------------------+
             call diffgeom_intersection( &                                  !
                  brep%edges(mesh%ids(ivert))%curve%surf, &                 !
                  mesh%uv(:,:,ivert), &                                     !
                  duv_ds, &                                                 !
                  tng, &                                                    !
                  stat )                                                    !
             ds = - dot_product(tng(:,1), grad(:,ivert)) / &                !
                  dot_product(tng(:,1), matmul(hess(:,:,ivert), tng(:,1)))  !
             duv(:,:,ivert) = ds * duv_ds(:,1,:)                            !
             !
             ! reproject onto exact intersection & handle passing to prev/next edge (w out of bounds..)
             uvtmp = mesh%uv(:,:,ivert) + duv(:,:,ivert)
             IF (.FALSE.) THEN
                call newton_intersection_polyline( &
                     brep%edges(mesh%ids(ivert))%curve%surf, &
                     lowerb, &
                     upperb, &
                     mesh%xyz(:,ivert), &
                     ds**2, &
                     stat, &
                     uvtmp, &
                     xyztmp )
                if ( stat > 0 ) then
                   PRINT *,'failed to reproject onto transversal intersection curve'
                   PAUSE
                end if
                mesh%uv(:,:,ivert) = uvtmp
                mesh%xyz(:,ivert) = xyztmp
             END IF
             !
          case (2) ! -------------------------------------------------------+
             do ivar = 1,2 ! <-------------------------------+              !
                call evald1( &                               !              !
                     tng(:,ivar), &                          !              !
                     brep%faces(mesh%ids(ivert))%surface, &  !              !
                     mesh%uv(:,1,ivert), &                   !              !
                     ivar )                                  !              !
             end do ! <--------------------------------------+              !
             call solve_NxN( &                                              !
                  duv(:,1,ivert), &                                         !
                  matmul(matmul(transpose(tng), hess(:,:,ivert)), tng), &   !
                  -matmul(transpose(tng), grad(:,ivert)), &                 !
                  singular )                                                !
             !
             ! handle passing to an adjacent face (crossing of a smooth edge...)
             mesh%uv(:,1,ivert) = mesh%uv(:,1,ivert) + duv(:,1,ivert)
             call eval( &
                  mesh%xyz(:,ivert), &
                  brep%faces(mesh%ids(ivert))%surface, &
                  mesh%uv(:,1,ivert) )
             !
          end select ! <----------------------------------------------------+


       end do


    end do
    
  end subroutine optim_jiao











  

  subroutine energy_jiao( &
       ntr, &
       tri, &
       wei, &
       nve, &
       xyz, &
       hve, &
       frac_conf, &
       ener, &
       grad, &
       hess )
    implicit none
    integer,       intent(in)  :: ntr            ! number of triangles
    integer,       intent(in)  :: tri(3,ntr)     ! triangles (vertex indices)
    real(kind=fp), intent(in)  :: wei(ntr)       ! triangle weights
    integer,       intent(in)  :: nve            ! number of vertices
    real(kind=fp), intent(in)  :: xyz(3,nve)     ! vertices (xyz coordinates)
    real(kind=fp), intent(in)  :: hve(nve)       ! target edge length at vertices
    real(kind=fp), intent(in)  :: frac_conf      ! fraction of conformal energy
    real(kind=fp), intent(out) :: ener(ntr)      ! triangle combined energy
    real(kind=fp), intent(out) :: grad(3,nve)    ! gradient at vertices
    real(kind=fp), intent(out) :: hess(3,3,nve)  ! Hessian matrix at vertices
    real(kind=fp)              :: frac_isom      ! fraction of isometric energy
    real(kind=fp)              :: Id4(3,3)       ! 4*identity(3,3)
    real(kind=fp)              :: nor_tri(3,ntr) ! triangle normals
    real(kind=fp)              :: area_act(ntr)  ! (twice) triangle area (actual)
    real(kind=fp)              :: area_ide(ntr)  ! (twice) triangle area (ideal)
    real(kind=fp)              :: htri           ! target edge length in a triangle
    real(kind=fp)              :: edg_vec(3,3)   ! edge vectors in a triangle
    real(kind=fp)              :: edg_len(3)     ! squared edge lengths in a triangle
    real(kind=fp)              :: edg_90d(3,3)   ! 90 degree CCW rotated edge vectors
    real(kind=fp)              :: ener_conf      ! conformal energy (angle-preserving)
    real(kind=fp)              :: ener_isom      ! isometric energy (area-preserving)
    real(kind=fp)              :: grad_conf(3)   ! gradient of conformal energy
    real(kind=fp)              :: grad_isom(3)   ! gradient of isometric energy
    real(kind=fp)              :: hess_conf(3,3) ! Hessian of conformal energy
    real(kind=fp)              :: hess_isom(3,3) ! Hessian of isometric energy
    integer                    :: itr, ive, i

    frac_isom = 1._fp - frac_conf
    Id4 = 4._fp * identity_matrix(3)

    ! compute triangle normals, actual and target areas
    do itr = 1,ntr
       ! (pseudo-)normal vector
       nor_tri(:,itr) = cross( &
            xyz(:,tri(2,itr)) - xyz(:,tri(1,itr)), &
            xyz(:,tri(3,itr)) - xyz(:,tri(1,itr)) )
       ! actual triangle area
       area_act(itr) = norm2(nor_tri(:,itr))
       nor_tri(:,itr) = nor_tri(:,itr) / area_act(itr) ! unit normal
       ! target edge length for current triangle
       htri = sum(hve(tri(:,itr))) / 3._fp
       ! ideal triangle area
       area_ide(itr) = htri**2
    end do
    ! normalize ideal triangle areas
    area_ide = area_ide * sum(area_act) / sum(area_ide)

    
    grad(:,:) = 0._fp
    hess(:,:,:) = 0._fp
    do itr = 1,ntr ! <-------------------------------------------------------------------------------------+
       ! edge vectors                                                                                      !
       edg_vec = xyz(:,tri([3,1,2],itr)) - xyz(:,tri([2,3,1],itr))                                         !
       ! squared edge lengths                                                                              !
       edg_len = sum(edg_vec**2, 1)                                                                        !
       ! 90 degree counter-clockwise rotated edge vectors                                                  !
       do i = 1,3 ! <-----------------------------------------+                                            !
          edg_90d(:,i) = cross(nor_tri(:,itr), edg_vec(:,i))  !                                            !
       end do ! <---------------------------------------------+                                            !
       ! triangle conformal energy                                                                         !
       ener_conf = sum(edg_len)/area_act(itr)                                                              !
       ! triangle isometric energy                                                                         !
       ener_isom = area_act(itr) / area_ide(itr)                                                           !
       ener_isom = ener_isom + 1._fp/ener_isom                                                             !
       ! triangle combined energy                                                                          !
       ener(itr) = frac_conf*ener_conf + frac_isom*ener_isom                                               !
       !                                                                                                   !
       ! vertex-based gradient and Hessian matrix of total energy                                          !
       do i = 1,3 ! <----------------------------------------------------------------------------------+   !
          ive = tri(i,itr)                                                                             !   !
          ! gradient of conformal energy                                                               !   !
          grad_conf = 2._fp*(edg_vec(:,1+mod(i,3)) - edg_vec(:,1+mod(i+1,3))) - ener_conf*edg_90d(:,i) !   !
          grad_conf = grad_conf / area_act(itr)                                                        !   !
          ! gradient of isometric energy                                                               !   !
          grad_isom = (area_act(itr)**2 - area_ide(itr)**2)*edg_90d(:,i)                               !   !
          grad_isom = grad_isom / ( area_ide(itr) * area_act(itr)**2 )                                 !   !
          ! gradient of total energy                                                                   !   !
          grad(:,ive) = grad(:,ive) + wei(itr) * ( frac_conf*grad_conf + frac_isom*grad_isom )         !   !
          !                                                                                            !   !
          ! Hessian matrix of conformal energy                                                         !   !
          hess_conf = outer_product(grad_conf, edg_90d(:,i))                                           !   !
          hess_conf = hess_conf + transpose(hess_conf)                                                 !   !
          hess_conf = (Id4 - hess_conf) / area_act(itr)                                                !   !
          ! Hessian matrix of isometric energy                                                         !   !
          hess_isom = 2._fp * area_ide(itr) * outer_product(edg_90d(:,i), edg_90d(:,i))                !   !
          hess_isom = hess_isom / area_act(itr)**3                                                     !   !
          ! Hessian matrix of combined energy                                                          !   !
          hess(:,:,ive) = hess(:,:,ive) + wei(itr) * ( frac_conf*hess_conf + frac_isom*hess_isom )     !   !
       end do ! <--------------------------------------------------------------------------------------+   !
       !                                                                                                   !
    end do ! <---------------------------------------------------------------------------------------------+

  end subroutine energy_jiao


end module mod_optimmesh
