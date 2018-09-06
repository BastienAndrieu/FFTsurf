subroutine diffgeom_intersection( &
     surf, &
     uv, &
     duv_ds, &
     dxyz_ds, &
     stat, &
     curvature )
  ! Computes tangent direction(s) and curvature of the intersection curve between
  ! two surfaces at a given intersection point.
  ! Returns: stat =-1 if one surface is singular at that point (normal = zero vector)
  !               = 0 if the point is on a single transveral intersection curve
  !               = 1 if the point is on a single tangential intersection curve
  !               = 2 if the point is a branch point (junction of two intersection curves)
  !               = 3 if the point is an isolated tangential contact point (degenerate curve)
  !               = 4 if the point is an high-order tangential contact point (possible degeneracy)
  use mod_math
  use mod_diffgeom
  ! ( see "Shape interrogation for computer aided design and manufacturing", &
  ! Patrikalakis et al. (2009), pp.166-175, 
  ! and "Tracing surface intersections with validated ODE system solver", Mukundan et al (2004)
  implicit none
  real(kind=fp), parameter             :: EPS = 1.d-12
  type(ptr_surface),       intent(in)  :: surf(2)
  real(kind=fp),           intent(in)  :: uv(2,2)
  real(kind=fp),           intent(out) :: duv_ds(2,2,2) ! u/v, #branch, #surf
  real(kind=fp),           intent(out) :: dxyz_ds(3,2)
  integer,                 intent(out) :: stat
  real(kind=fp), optional, intent(out) :: curvature(2)
  real(kind=fp)                        :: dxyz_duv(3,2,2), n(3,2), d2xyz_duv2(3)
  real(kind=fp)                        :: dotn, sgndotn, aux(3), w
  real(kind=fp), dimension(3,2)        :: EFG, LMN
  real(kind=fp), dimension(2)          :: detEFG, invsqrtdetEFG
  real(kind=fp)                        :: A(2,2), B(3), discr, invdenom, kn(2)
  integer                              :: isurf, jsurf, ivar, i, j

  do isurf = 1,2 ! <------------------------------------------------------------------+
     ! tangent vectors                                                                !
     do ivar = 1,2 ! <-----------------------------+                                  !
        call evald1( dxyz_duv(:,ivar,isurf), &     !                                  !
             surf(isurf)%ptr, &                    !                                  !
             uv(:,isurf), &                        !                                  !
             ivar )                                !                                  !
     end do ! <------------------------------------+                                  !
     ! First Fundamental Form coefficients                                            !
     EFG(1,isurf) = dot_product(dxyz_duv(:,1,isurf), dxyz_duv(:,1,isurf))             !
     EFG(2,isurf) = dot_product(dxyz_duv(:,1,isurf), dxyz_duv(:,2,isurf))             !
     EFG(3,isurf) = dot_product(dxyz_duv(:,2,isurf), dxyz_duv(:,2,isurf))             !
     detEFG(isurf) = EFG(1,isurf)*EFG(3,isurf) - EFG(2,isurf)**2                      !
     invsqrtdetEFG(isurf) = sqrt(1._fp / detEFG(isurf))                               !
     !                                                                                !
     ! unit normal vector                                                             !
     n(:,isurf) = cross(dxyz_duv(:,1,isurf), dxyz_duv(:,2,isurf))                     !
     !PRINT *,'N =',N(:,ISURF)
     !n(:,isurf) = n(:,isurf) * invsqrtdetEFG(isurf)
     !n(:,isurf) = n(:,isurf) / norm2(n(:,isurf))
  end do ! <--------------------------------------------------------------------------+
  !PRINT *,'DXYZ_DUV(1) ='
  !CALL PRINT_MAT(DXYZ_DUV(:,:,1))
  !PRINT *,'EFG ='
  !CALL PRINT_MAT(EFG)
  !PRINT *,'DETEFG =',DETEFG
  !PRINT *,'INVSQRTDETEFG =',INVSQRTDETEFG
  !PRINT *,'N ='
  !CALL PRINT_MAT(N)
  dotn = dot_product(n(:,1), n(:,2))
  !PRINT *,'DOTN =',DOTN

  if ( abs(dotn) > (1._fp - EPS)*product(sqrt(detEFG)) ) then ! <---------------------+
     !! tangential intersection                                                       !
     sgndotn = sign(1._fp, dotn)                                                      !
     ! Second Fundamental Form coefficients                                           !
     do isurf = 1,2 ! <-----------------------------------------------+               !
        do ivar = 1,3 ! <-----------------------------------------+   !               !
           call evald2( d2xyz_duv2, &                             !   !               !
                surf(isurf)%ptr, &                                !   !               !
                uv(:,isurf), &                                    !   !               !
                ivar )                                            !   !               !
           LMN(ivar,isurf) = dot_product(n(:,isurf), d2xyz_duv2)  !   !               !
        end do ! <------------------------------------------------+   !               !
     end do ! <-------------------------------------------------------+               !
     where( abs(LMN) < EPSmath ) LMN = 0._fp
     do isurf = 1,2
        LMN(:,isurf) = LMN(:,isurf) * invsqrtdetEFG(isurf)
     end do
     !PRINT *,'DETLMN =',LMN(1,:)*LMN(3,:) - LMN(2,:)**2
     !PRINT *,'LMN ='
     !CALL PRINT_MAT(LMN)
     jsurf = minloc(detEFG, 1)
     isurf = 1+mod(jsurf,2)
     ! linear realtion duvj_ds = A * duvi_ds                                          !
     do i = 1,2 ! <---------------------------------------------------+               !
        aux = real((-1)**i, kind=fp) * &                              !               !
             cross(n(:,jsurf), dxyz_duv(:,1+mod(i,2),isurf))          !               !
        !PRINT *,'AUX =',AUX
        do j = 1,2 ! <-------------------------------------+          !               !
           A(i,j) = dot_product(aux, dxyz_duv(:,j,jsurf))  !          !               !
        end do ! <-----------------------------------------+          !               !
     end do ! <-------------------------------------------------------+               !
     !A = sign(invsqrtdetEFG(jsurf), dotn)*A
     where( abs(A) < EPSmath ) A = 0._fp
     A = A / dotn
     !PRINT *,'A ='
     !CALL PRINT_MAT(A)
     ! quadratic equation for duv1_ds                                                 !
     B(1) = LMN(1,isurf)*A(1,1)**2 + 2._fp*LMN(2,isurf)*A(1,1)*A(2,1) + LMN(3,isurf)*A(2,1)**2
     B(2) = LMN(1,isurf)*A(1,1)*A(1,2) + LMN(2,isurf)*(A(1,2)*A(2,1) + A(1,1)*A(2,2)) + LMN(3,isurf)*A(2,1)*A(2,2)
     B(3) = LMN(1,isurf)*A(1,2)**2 + 2._fp*LMN(2,isurf)*A(1,2)*A(2,2) + LMN(3,isurf)*A(2,2)**2
     !B = B / detEFG(jsurf)
     B = LMN(:,jsurf) - sign(1._fp, dotn)*B
     where( abs(B) < EPSmath ) B = 0._fp
     !PRINT *,'B =',B
     discr = B(2)**2 - B(1)*B(3)                                                      !
     !A = sign(invsqrtdetEFG(jsurf), dotn)*A
     !PRINT *,'DISCR =',DISCR
     !                                                                                !
     if ( all(abs(B) < EPS) ) then ! <-------------------------------------+          !
        ! high-order contact point (dxyz_ds cannot be determined)          !          !
        stat = 4                                                           !          !
        return                                                             !          !
        !                                                                  !          !
     elseif ( discr < -EPS ) then ! ---------------------------------------+          !
        ! isolated tangential contact point (dxyz_ds undefined)            !          !
        stat = 3                                                           !          !
        return                                                             !          !
        !                                                                  !          !
     else ! ---------------------------------------------------------------+          !
        if ( discr < EPS ) then ! <------------------------------+         !          !
           ! tangential intersection curve (one single dxyz_ds)  !         !          !
           stat = 1                                              !         !          !
           discr = 0._fp                                         !         !          !
        else ! --------------------------------------------------+         !          !
           ! branch point (two disctinct dxyz_ds)                !         !          !
           stat = 2                                              !         !          !
           discr = sqrt(discr)                                   !         !          !
        end if ! <-----------------------------------------------+         !          !
        !                                                                  !          !
        ! tangent direction(s) to the intersection branch(es)              !          !
        ! (the number of branches is equal to stat)                        !          !
        if ( abs(B(1)) > abs(B(3)) ) then ! <---------------------------+  !          !
           do i = 1,stat ! <-----------------------------------------+  !  !          !
              w = (discr*real((-1)**i, kind=fp) - B(2)) / B(1)       !  !  !          !
              duv_ds(:,i,jsurf) = [w, 1._fp] / &                     !  !  !          !
                   sqrt(EFG(1,jsurf)*w**2 + &                        !  !  !          !
                   2._fp*EFG(2,jsurf)*w + &                          !  !  !          !
                   EFG(3,jsurf))                                     !  !  !          !
           end do ! <------------------------------------------------+  !  !          !
        else ! ---------------------------------------------------------+  !          !
           do i = 1,stat ! <-----------------------------------------+  !  !          !
              w = (discr*real((-1)**i, kind=fp) - B(2)) / B(3)       !  !  !          !
              duv_ds(:,i,jsurf) = [1._fp, w] / &                     !  !  !          !
                   sqrt(EFG(1,jsurf) + &                             !  !  !          !
                   2._fp*EFG(2,jsurf)*w + &                          !  !  !          !
                   EFG(3,jsurf)*w**2)                                !  !  !          !
           end do ! <------------------------------------------------+  !  !          !
        end if ! <------------------------------------------------------+  !          !
        do i = 1,stat ! <--------------------------------------------+     !          !
           duv_ds(:,i,isurf) = matmul(A, duv_ds(:,i,jsurf))          !     !          !
           dxyz_ds(:,i) = matmul(dxyz_duv(:,:,jsurf), &              !     !          !
                duv_ds(:,i,jsurf))                                   !     !          !
        end do ! <---------------------------------------------------+     !          !
        !                                                                  !          !
     end if ! <------------------------------------------------------------+          !
     !                                                                                !
     !                                                                                !
  else ! -----------------------------------------------------------------------------+
     !! transversal intersection                                                      !
     stat = 0                                                                         !
     invdenom = 1._fp * sqrt(1._fp - dotn**2)                                         !
     dxyz_ds(:,1) = cross(n(:,1), n(:,2))                                             !
     dxyz_ds(:,1) = dxyz_ds(:,1) / norm2(dxyz_ds(:,1))                                !
     !                                                                                !
     do isurf = 1,2 ! <-----------------------------------------------+               !
        aux(1) = dot_product(dxyz_ds(:,1), dxyz_duv(:,1,isurf))       !               !
        aux(2) = dot_product(dxyz_ds(:,1), dxyz_duv(:,2,isurf))       !               !
        duv_ds(1,1,isurf) = aux(1)*EFG(3,isurf) - aux(2)*EFG(2,isurf) !               !
        duv_ds(2,1,isurf) = aux(2)*EFG(1,isurf) - aux(1)*EFG(2,isurf) !               !
        duv_ds(:,1,isurf) = duv_ds(:,1,isurf) / detEFG(isurf)         !               !
     end do ! <-------------------------------------------------------+               !
     !                                                                                !
     if ( present(curvature) ) then ! <--------------------------------------------+  !
        do isurf = 1,2 ! <--------------------------------------------------+      !  !
           ! Second Fundamental Form coefficients                           !      !  !
           do ivar = 1,3 ! <-----------------------------------------+      !      !  !
              call evald2( d2xyz_duv2, &                             !      !      !  !
                   surf(isurf)%ptr, &                                !      !      !  !
                   uv(:,isurf), &                                    !      !      !  !
                   ivar )                                            !      !      !  !
              LMN(ivar,isurf) = dot_product(n(:,isurf), d2xyz_duv2)  !      !      !  !
           end do ! <------------------------------------------------+      !      !  !
           LMN(:,isurf) = LMN(:,isurf) * invsqrtdetEFG(isurf)               !      !  !
           !                                                                !      !  !
           ! normal curvature                                               !      !  !
           kn(isurf) =  &                                                   !      !  !
                LMN(1,isurf)*duv_ds(1,1,isurf)**2 + &                       !      !  !
                2._fp*LMN(2,isurf)*duv_ds(1,1,isurf)*duv_ds(2,1,isurf) + &  !      !  !
                LMN(3,isurf)*duv_ds(2,1,isurf)**2                           !      !  !
        end do ! <----------------------------------------------------------+      !  !
        !                                                                          !  !
        curvature(1) = sqrt(kn(1)**2 + kn(2)**2 - 2._fp*kn(1)*kn(2)*dotn)*invdenom !  !
     end if ! <--------------------------------------------------------------------+  !
  end if ! <--------------------------------------------------------------------------+
  
end subroutine diffgeom_intersection
