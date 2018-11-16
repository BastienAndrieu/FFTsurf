module mod_propagation

  use mod_options
  
  implicit none

contains

  subroutine velocity_vectors_unit_normal_speed( &
       xyz, &
       m, &
       n, &
       vvel )
    use mod_math
    use mod_polynomial
    use mod_chebyshev
    implicit none
    integer,       intent(in)       :: m, n
    real(kind=fp), intent(in)       :: xyz(m,n,3)
    real(kind=fp), intent(out)      :: vvel(m,n,3)
    type(type_polynomial)           :: p, pu, pv
    real(kind=fp), dimension(m,n,3) :: xu, xv
    integer                         :: i, j, k

    ! direct transform -> position vector polynomial
    call reset_polynomial(p, 2, 1, [m,n]-1, 3)
    call fcht2(xyz, p%coef, m, n, 3, EPSmath)

    ! compute derivative (tangents) polynomials
    call diff2(p, pu, pv)

    ! evaluate derivatives (tangents) at CGL nodes
    call ifcht2(pu%coef, xu, m, n, 3)
    call ifcht2(pv%coef, xv, m, n, 3)

    ! pseudonormal vector
    do i = 1,3
       j = 1 + mod(i,3)
       k = 1 + mod(j,3)
       vvel(1:m,1:n,i) = xu(1:m,1:n,j)*xv(1:m,1:n,k) - xu(1:m,1:n,k)*xv(1:m,1:n,j)
    end do
    ! unitize
    vvel = vvel / spread(sqrt(sum(vvel**2,3)), dim=3, ncopies=3)
    
  end subroutine velocity_vectors_unit_normal_speed




  
  
  subroutine velocity_vectors_rotation_z( &
       xyz, &
       m, &
       n, &
       vvel )
    use mod_math
    implicit none
    integer,       intent(in)  :: m, n
    real(kind=fp), intent(in)  :: xyz(m,n,3)
    real(kind=fp), intent(out) :: vvel(m,n,3)

    vvel(1:m,1:n,1) = -xyz(1:m,1:n,2)
    vvel(1:m,1:n,2) =  xyz(1:m,1:n,1)
    vvel(1:m,1:n,3) =  0._fp
    vvel = vvel / CSTpi

  end subroutine velocity_vectors_rotation_z

  
  
  subroutine velocity_vectors_vortex( &
       xyz, &
       m, &
       n, &
       time, &
       vvel )
    use mod_math
    implicit none
    real(kind=fp), parameter   :: period = 1._fp
    integer,       intent(in)  :: m, n
    real(kind=fp), intent(in)  :: xyz(m,n,3)
    real(kind=fp), intent(in)  :: time
    real(kind=fp), intent(out) :: vvel(m,n,3)
    integer                    :: i, j, k
    
    do i = 1,3
       j = 1 + mod(i,3)
       k = 1 + mod(j,3)
       vvel(1:m,1:n,i) = sin(CSTpi * xyz(1:m,1:n,i))**2 * &
            ( sin(2._fp * CSTpi * xyz(1:m,1:n,k)) - sin(2._fp * CSTpi * xyz(1:m,1:n,j)) )
    end do
    vvel = vvel * cos(CSTpi * time/period)

    !IF ( TIME > 0.15 ) VVEL = 0._FP

  end subroutine velocity_vectors_vortex





  subroutine propagation_step_RK4( &
       surf, &
       proplaw, &
       timestep, &
       time )
    use mod_diffgeom
    use mod_chebyshev
    implicit none
    real(kind=fp), parameter          :: fracstep(4) = [1._fp, 0.5_fp, 0.5_fp, 1._fp]
    type(type_surface), intent(inout) :: surf
    integer,            intent(in)    :: proplaw
    real,               intent(in)    :: timestep
    real,               intent(in)    :: time
    real(kind=fp), allocatable        :: vvel(:,:,:,:), xyz(:,:,:)
    integer                           :: m, n
    real(kind=fp)                     :: dti, ti, dt6
    integer                           :: istep

    m = size(surf%x%coef,1)
    n = size(surf%x%coef,2)

    ! position vector (inverse transform)
    allocate(xyz(m,n,3))
    call ifcht2( &
         surf%x%coef(1:m,1:n,1:3), &
         xyz, &
         m, &
         n, &
         3 )

    ! velocity vector
    allocate(vvel(m,n,3,0:4))
    vvel(1:m,1:n,1:3,0) = 0._fp
    
    ! intermediate steps
    IF ( .true. ) THEN
       do istep = 1,4
          dti = real(timestep, kind=fp) * fracstep(istep)
          ti = real(time, kind=fp) + dti
          select case ( proplaw )
          case (1) ! unit normal speed
             call velocity_vectors_unit_normal_speed( &
                  xyz + dti*vvel(1:m,1:n,1:3,istep-1), &
                  m, &
                  n, &
                  vvel(1:m,1:n,1:3,istep) )
             !
          case (2) ! vortex
             call velocity_vectors_vortex( &
                  xyz + dti*vvel(1:m,1:n,1:3,istep-1), &
                  m, &
                  n, &
                  ti, &
                  vvel(1:m,1:n,1:3,istep) )
             !
          case (3) ! rotation around z-axis
             call velocity_vectors_rotation_z( &
                  xyz + dti*vvel(1:m,1:n,1:3,istep-1), &
                  m, &
                  n, &
                  vvel(1:m,1:n,1:3,istep) )
             !
          end select
       end do

       ! time-integration
       dt6 = real(timestep, kind=fp) / 6._fp
       do istep = 1,4
          xyz(1:m,1:n,1:3) = xyz(1:m,1:n,1:3) + (dt6/fracstep(istep)) * vvel(1:m,1:n,1:3,istep)
       end do
    ELSE
       call velocity_vectors_vortex( &
            xyz, &
            m, &
            n, &
            real(time, kind=fp), &
            vvel(1:m,1:n,1:3,1) )
       xyz(1:m,1:n,1:3) = xyz(1:m,1:n,1:3) + real(timestep, kind=fp) * vvel(1:m,1:n,1:3,1)
    END IF
    deallocate(vvel)
    
    ! new Chebyshev coefficients (direct transform)
    call fcht2( &
         xyz, &
         surf%x%coef(1:m,1:n,1:3), &
         m, &
         n, &
         3, &
         EPSmath )
    deallocate(xyz)

    call economize2( &
            surf%x, &
            EPSmath )

    ! new derivatives
    call compute_deriv1(surf)
    call compute_deriv2(surf)
    call compute_pseudonormal(surf)
    call economize2( &
            surf%pn, &
            EPSmath )

  end subroutine propagation_step_RK4







  
  
  
  
  subroutine eos_surface( &
       xyz, &
       speed, &
       m, &
       n, &
       timestep, &
       sens_normal, &
       w )
    use mod_math
    use mod_polynomial
    use mod_chebyshev
    implicit none
    integer,       intent(in)       :: m, n
    real(kind=fp), intent(in)       :: xyz(m,n,3)
    real(kind=fp), intent(in)       :: speed(m,n)
    real(kind=fp), intent(in)       :: timestep
    real(kind=fp), intent(out)      :: w(m,n,3) ! unit pointwise direction to envelope of sphere
    integer,       intent(in)       :: sens_normal ! +-1
    type(type_polynomial)           :: px, pxu, pxv, ps, psu, psv
    real(kind=fp), dimension(m,n,3) :: xu, xv, nor
    real(kind=fp), dimension(m,n)   :: su, sv, wsqr
    real(kind=fp)                   :: dtsqr, wij(3), mat(3,3), rhs(3)
    logical                         :: singular
    integer                         :: i, j, k

    call reset_polynomial(px, 2, 1, [m,n]-1, 3)
    call fcht2(xyz, px%coef, m, n, 3, EPSmath)
    call diff2(px, pxu, pxv)
    call ifcht2(pxu%coef, xu, m, n, 3)
    call ifcht2(pxv%coef, xv, m, n, 3)

    do i = 1,3
       j = 1 + mod(i,3)
       k = 1 + mod(j,3)
       nor(1:m,1:n,i) = xu(1:m,1:n,j)*xv(1:m,1:n,k) - xu(1:m,1:n,k)*xv(1:m,1:n,j)
    end do
    nor = nor / spread(sqrt(sum(nor**2,3)), dim=3, ncopies=3)

    call reset_polynomial(ps, 2, 1, [m,n]-1, 1)
    call fcht2(speed, ps%coef, m, n, 1, EPSmath)
    call ifcht2(psu%coef, su, m, n, 1)
    call ifcht2(psv%coef, sv, m, n, 1)

    dtsqr = timestep**2
    do j = 1,n
       do i = 1,m
          mat(1:3,1) = xu(i,j,1:3)
          mat(1:3,2) = xv(i,j,1:3)
          mat(1:3,3) = nor(i,j,1:3)
          rhs = -[su(i,j), sv(i,j), 0._fp]

          call solve_NxN( &
               w(i,j,1:3), &
               mat, &
               rhs, &
               singular )

          w(i,j,1:3) = wij(1:3)
          wsqr(i,j) = sum(wij**2)

          if ( dtsqr*wsqr(i,j) > 1._fp ) then
             STOP 'envelope_of_spheres : dt * |w| > 1  ==> no real solution.'
          end if
       end do
    end do

    w = timestep * w + &                                    ! tangential component
         nor * real(sens_normal, kind=fp) * &               ! normal...
         spread(sqrt(1._fp - dtsqr*wsqr), dim=3, ncopies=3) ! ...component

  end subroutine eos_surface


end module mod_propagation
