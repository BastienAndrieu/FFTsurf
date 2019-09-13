program singular_surface_surface_intersection

  use mod_math
  use mod_polynomial
  use mod_diffgeom
  use mod_linalg
  use mod_util
  use mod_obb
  use mod_intersection

  integer, parameter :: itmax = 20
  real(kind=fp)      :: perturbation, tolerance, scale
  character(99)      :: arg
  type(ptr_surface)  :: surf(2)
  type(type_obb)     :: obb(2)
  real(kind=fp)      :: uv0(2,2), uv(2,2), duvr(2,2)
  integer            :: isurf
  real(kind=fp)      :: beta, histo(3,itmax)
  integer            :: it, method
  real(kind=fp)      :: duv_ds(2,2,2), dxyz_ds(3,2), dxyz_duv(3,2), nor(3,2)
  integer            :: stat
  integer            :: fid
  character          :: strnum



  if ( command_argument_count() < 1 ) then
     perturbation = 0._fp
  else
     call get_command_argument(1, arg)
     read (arg, *) perturbation
  end if

  if ( command_argument_count() < 2 ) then
     scale = 1._fp
  else
     call get_command_argument(2, arg)
     read (arg, *) scale
  end if

  if ( command_argument_count() < 3 ) then
     tolerance = real(1.d-12, kind=fp)
  else
     call get_command_argument(3, arg)
     read (arg, *) tolerance
  end if
  
  
  beta = huge(1._fp)
  do isurf = 1,2
     write (strnum, '(i1)') isurf
     allocate(surf(isurf)%ptr)
     call read_polynomial( &
          surf(isurf)%ptr%x, &
          'singular_surface_surface_intersection/surf'//strnum//'.cheb', &
          nvar=2, &
          base=1 )
     call economize2( &
          surf(isurf)%ptr%x, &
          epsilon(1._fp) )
     surf(isurf)%ptr%x%coef = scale*surf(isurf)%ptr%x%coef
     
     call compute_deriv1(surf(isurf)%ptr)
     call compute_deriv2(surf(isurf)%ptr)
     call compute_pseudonormal(surf(isurf)%ptr)

     call chebOBB2( &
          surf(isurf)%ptr%x%coef, &
          surf(isurf)%ptr%x%degr, &
          obb(isurf) )
     beta = min(beta, 2._fp*norm2(obb(isurf)%rng))
  end do


  call get_free_unit(fid)
  open(unit=fid, file='singular_surface_surface_intersection/uv0.dat', action='read')
  do isurf = 1,2
     read (fid,*) uv0(1:2,isurf)
  end do
  close(fid)


  !uv0(:,:) = 0._fp


  call random_number(duvr)
  uv0 = min(1._fp, max(-1._fp, uv0 + perturbation*(2._fp*duvr - 1._fp)))

  print *, 'perturbation uv =', perturbation
  print *, 'scale xyz =', scale
  print *, 'tolerance =', tolerance
  print *, 'uv0 =', uv0
  print *, ''

  do method = 1,3
     print *, ''
     print *, '*********************************************'
     print *, ''
     uv = uv0
     select case (method)
     case (1)
        print *, 'Newton (scaled, beta =', beta, ')'
        call newton_scaled( &
             uv, &
             surf, &
             beta, &
             tolerance, &
             itmax, &
             histo, &
             it )
     case (2)
        print *, 'Newton (basic)'
        call newton_basic( &
             uv, &
             surf, &
             tolerance, &
             itmax, &
             histo, &
             it )
     case (3)
        print *, 'Simult. pt. inv.'
        call simult_pt_inv( &
             uv, &
             surf, &
             tolerance, &
             itmax, &
             histo, &
             it )
     end select

     print *, '--------------------------+-------------------------+------------------------'
     print *, '         |S1 - S2|        |          |DUV1|         |         |DUV2|         '
     print *, '--------------------------+-------------------------+------------------------'
     call print_mat(transpose(histo(:,1:min(it,itmax))))
     print *, ''
     print *, 'uv =', uv
     print *, '|uv - uv0| =', &
          norm2(uv(1:2,1) - uv0(1:2,1)), &
          norm2(uv(1:2,2) - uv0(1:2,2)), &
          sqrt(sum((uv - uv0)**2))

     !if ( .true. ) then
     do isurf = 1,2
        do ivar = 1,2
           call evald1( &
                dxyz_duv(1:3,ivar), &
                surf(isurf)%ptr, &
                uv(1:2,isurf), &
                ivar )
        end do
        nor(1:3,isurf) = cross(dxyz_duv(1:3,1), dxyz_duv(1:3,2))
        nor(1:3,isurf) = nor(1:3,isurf)/norm2(nor(1:3,isurf))
     end do
     !print *, 'n1.n2 = ', dot_product(nor(1:3,1), nor(1:3,2))
     print *, '|1 - |n1.n2|| =', abs(1._fp - abs(dot_product(nor(1:3,1), nor(1:3,2))))
     !else
     call diffgeom_intersection( &
          surf, &
          uv, &
          duv_ds, &
          dxyz_ds, &
          stat, &
          normalize_scale=.true. )
     print *, 'stat =', stat
     if ( stat > 0 .and. stat < 3 ) then
        print *,'dxyz_ds ='
        call print_mat(transpose(dxyz_ds(1:3,1:stat)))
     end if
     !end if
  end do

contains

  subroutine newton_scaled( &
       uv, &
       surf, &
       beta, &
       tolerance, &
       itmax, &
       histo, &
       it )
    implicit none
    real(kind=fp),     intent(inout) :: uv(2,2)
    type(ptr_surface), intent(in)    :: surf(2)
    real(kind=fp),     intent(in)    :: beta
    real(kind=fp),     intent(in)    :: tolerance
    integer,           intent(in)    :: itmax
    real(kind=fp),     intent(out)   :: histo(3,itmax)
    integer,           intent(out)   :: it
    real(kind=fp)                    :: xyz(3,2), dxyz_duv(3,2,2), d2xyz_duv2(3,3,2)
    real(kind=fp)                    :: norm_dxyz_duv(2,2)
    real(kind=fp)                    :: tangents(3,2,2), n1(3)
    real(kind=fp)                    :: dn1_duv(3,2), dt2_duv(3,2,2)
    real(kind=fp)                    :: f(5), jac(5,4)
    real(kind=fp)                    :: duv(4)
    integer                          :: ivar, jvar, kvar, isurf

    jac(1:5,1:4) = 0._fp
    do it = 1,itmax

       do isurf = 1,2
          call eval( &
               xyz(1:3,isurf), &
               surf(isurf)%ptr, &
               uv(1:2,isurf) )

          do ivar = 1,2
             call evald1( &
                  dxyz_duv(1:3,ivar,isurf), &
                  surf(isurf)%ptr, &
                  uv(1:2,isurf), &
                  ivar )
             norm_dxyz_duv(ivar,isurf) = norm2(dxyz_duv(1:3,ivar,isurf))
             tangents(1:3,ivar,isurf) = dxyz_duv(1:3,ivar,isurf)/norm_dxyz_duv(ivar,isurf)
          end do

          do ivar = 1,3
             call evald2( &
                  d2xyz_duv2(1:3,ivar,isurf), &
                  surf(isurf)%ptr, &
                  uv(1:2,isurf), &
                  ivar )
          end do
       end do

       n1(1:3) = cross(tangents(1:3,1,1), tangents(1:3,1,2))

       ! Residual vector
       f(1:3) = xyz(1:3,1) - xyz(1:3,2)
       histo(1,it) = norm2(f(1:3))
       f(4) = beta*dot_product(n1, tangents(1:3,1,2))
       f(5) = beta*dot_product(n1, tangents(1:3,2,2))

       ! Jacobian matrix
       ! FORMULE FAUSSE POUR DN1_DUV MAIS FONCTIONNE !?!?
       do ivar = 1,2
         dn1_duv(1:3,ivar) = cross( &
                    (d2xyz_duv2(1:3,ivar,1) - &!tangents(1:3,1,1)*&
                    dot_product(d2xyz_duv2(1:3,ivar,1), tangents(1:3,1,1)))/norm_dxyz_duv(1,1), &
                    tangents(1:3,2,1) &
                    ) +  cross( &
                    tangents(1:3,1,1), &
                    (d2xyz_duv2(1:3,ivar+1,1) - &!tangents(1:3,2,1)*&
                    dot_product(d2xyz_duv2(1:3,ivar+1,1), tangents(1:3,2,1)))/norm_dxyz_duv(2,1) &
                    )
       end do

       do jvar = 1,2
         do ivar = 1,2
            kvar = ivar + jvar - 1
            dt2_duv(1:3,ivar,jvar) = ( &
            d2xyz_duv2(1:3,kvar,2) - &
            tangents(1:3,ivar,2)*dot_product(d2xyz_duv2(1:3,kvar,2), tangents(1:3,ivar,2)) &
            )/norm_dxyz_duv(ivar,2)
         end do
       end do

       jac(1:3,1:2) = dxyz_duv(1:3,1:2,1)
       jac(1:3,3:4) = -dxyz_duv(1:3,1:2,2)
       do jvar = 1,2
         do ivar = 1,2
            jac(3+ivar,jvar) = beta*dot_product(dn1_duv(1:3,jvar), tangents(1:3,ivar,2))
            jac(3+ivar,2+jvar) = beta*dot_product(n1, dt2_duv(1:3,ivar,jvar))
         end do
       end do

       ! solve
       call linsolve_svd( &
            duv, &
            jac, &
            -f, &
            5, &
            4, &
            1, &
            tol=tolerance )

       histo(2,it) = norm2(duv(1:2))
       histo(3,it) = norm2(duv(3:4))

       uv(1:2,1) = uv(1:2,1) + duv(1:2)
       uv(1:2,2) = uv(1:2,2) + duv(3:4)

       if ( max(histo(2,it), histo(3,it)) < 1.d-12 ) return
    end do

  end subroutine newton_scaled






  subroutine newton_basic( &
       uv, &
       surf, &
       tolerance, &
       itmax, &
       histo, &
       it )
    implicit none
    real(kind=fp),     intent(inout) :: uv(2,2)
    type(ptr_surface), intent(in)    :: surf(2)
    real(kind=fp),     intent(in)    :: tolerance
    integer,           intent(in)    :: itmax
    real(kind=fp),     intent(out)   :: histo(3,itmax)
    integer,           intent(out)   :: it
    real(kind=fp)                    :: xyz(3,2), dxyz_duv(3,2,2), d2xyz_duv2(3,3,2)
    real(kind=fp)                    :: n1(3), dn1_duv(3,2)
    real(kind=fp)                    :: f(5), jac(5,4)
    real(kind=fp)                    :: duv(4)
    integer                          :: ivar, isurf

    do it = 1,itmax

       do isurf = 1,2
          call eval( &
               xyz(1:3,isurf), &
               surf(isurf)%ptr, &
               uv(1:2,isurf) )

          do ivar = 1,2
             call evald1( &
                  dxyz_duv(1:3,ivar,isurf), &
                  surf(isurf)%ptr, &
                  uv(1:2,isurf), &
                  ivar )
          end do

          do ivar = 1,3
             call evald2( &
                  d2xyz_duv2(1:3,ivar,isurf), &
                  surf(isurf)%ptr, &
                  uv(1:2,isurf), &
                  ivar )
          end do
       end do
       n1(1:3) = cross(dxyz_duv(1:3,1,1), dxyz_duv(1:3,1,2))


       ! Residual vector
       f(1:3) = xyz(1:3,1) - xyz(1:3,2)
       histo(1,it) = norm2(f(1:3))
       f(4) = dot_product(n1, dxyz_duv(1:3,1,2))
       f(5) = dot_product(n1, dxyz_duv(1:3,2,2))

       ! Jacobian matrix
       dn1_duv(1:3,1) = &
            cross(d2xyz_duv2(1:3,1,1), dxyz_duv(1:3,2,1)) + &
            cross(dxyz_duv(1:3,1,1), d2xyz_duv2(1:3,2,1))
       dn1_duv(1:3,2) = &
            cross(d2xyz_duv2(1:3,2,1), dxyz_duv(1:3,2,1)) + &
            cross(dxyz_duv(1:3,1,1), d2xyz_duv2(1:3,3,1))
            
       jac(1:5,1:4) = 0._fp
       jac(1:3,1:2) = dxyz_duv(1:3,1:2,1)
       jac(1:3,3:4) = -dxyz_duv(1:3,1:2,2)
       jac(4,1) = dot_product(dn1_duv(1:3,1), dxyz_duv(1:3,1,2))
       jac(4,2) = dot_product(dn1_duv(1:3,2), dxyz_duv(1:3,1,2))
       jac(4,3) = dot_product(n1, d2xyz_duv2(1:3,1,2))
       jac(4,4) = dot_product(n1, d2xyz_duv2(1:3,2,2))
       jac(5,1) = dot_product(dn1_duv(1:3,1), dxyz_duv(1:3,2,2))
       jac(5,2) = dot_product(dn1_duv(1:3,2), dxyz_duv(1:3,2,2))
       jac(5,3) = dot_product(n1, d2xyz_duv2(1:3,2,2))
       jac(5,4) = dot_product(n1, d2xyz_duv2(1:3,3,2))
       
       ! solve
       call linsolve_svd( &
            duv, &
            jac, &
            -f, &
            5, &
            4, &
            1, &
            tol=tolerance )

       histo(2,it) = norm2(duv(1:2))
       histo(3,it) = norm2(duv(3:4))

       uv(1:2,1) = uv(1:2,1) + duv(1:2)
       uv(1:2,2) = uv(1:2,2) + duv(3:4)

       if ( max(histo(2,it), histo(3,it)) < 1.d-12 ) return
       
    end do
  end subroutine newton_basic


  


  subroutine simult_pt_inv( &
       uv, &
       surf, &
       tolerance, &
       itmax, &
       histo, &
       it )
    implicit none
    real(kind=fp),     intent(inout) :: uv(2,2)
    type(ptr_surface), intent(in)    :: surf(2)
    real(kind=fp),     intent(in)    :: tolerance
    integer,           intent(in)    :: itmax
    real(kind=fp),     intent(out)   :: histo(3,itmax)
    integer,           intent(out)   :: it
    real(kind=fp)                    :: xyz(3,2), dxyz_duv(3,2)
    real(kind=fp)                    :: r(3), rhs(2), mat(2,2)
    real(kind=fp)                    :: duv(2)
    integer                          :: ivar, isurf

    do it = 1,itmax
       do isurf = 1,2
          call eval( &
               xyz(1:3,isurf), &
               surf(isurf)%ptr, &
               uv(1:2,isurf) )
       end do

       r = 0.5_fp*(xyz(1:3,1) - xyz(1:3,2))

       do isurf = 1,2
          do ivar = 1,2
             call evald1( &
                  dxyz_duv(1:3,ivar), &
                  surf(isurf)%ptr, &
                  uv(1:2,isurf), &
                  ivar )
          end do
          mat = matmul(transpose(dxyz_duv), dxyz_duv)
          rhs = matmul(transpose(dxyz_duv), real((-1)**isurf, kind=fp)*r)
          call linsolve_svd( &
               duv, &
               mat, &
               rhs, &
               2, &
               2, &
               1, &
               tol=tolerance )

          histo(1+isurf,it) = norm2(duv)
          uv(1:2,isurf) = uv(1:2,isurf) + duv(1:2)
       end do

       if ( max(histo(2,it), histo(3,it)) < 1.d-12 ) return
    end do

  end subroutine simult_pt_inv
















  !subroutine diffgeom_singular_intersection( &
  !     surf, &
  !     uv, &
  !     duv_ds, &
  !     dxyz_ds, &
  !     stat )
  !  ! Computes tangent direction(s) and curvature of the intersection curve between
  !  ! two surfaces at a given intersection point.
  !  ! Returns: stat =-1 if one surface is singular at that point (normal = zero vector)
  !  !               = 0 if the point is on a single transveral intersection curve
  !  !               = 1 if the point is on a single tangential intersection curve
  !  !               = 2 if the point is a branch point (junction of two intersection curves)
  !  !               = 3 if the point is an isolated tangential contact point (degenerate curve)
  !  !               = 4 if the point is an high-order tangential contact point (possible degeneracy)
  !  use mod_math
  !  use mod_diffgeom
  !end subroutine diffgeom_singular_intersection










  


  
end program singular_surface_surface_intersection
