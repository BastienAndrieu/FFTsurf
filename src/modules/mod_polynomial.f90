module mod_polynomial

  use mod_math
  use mod_chebyshev2
  use mod_bernstein2

  implicit none

  type type_polynomial
     integer                        :: nvar = 0
     integer                        :: base = 0
     integer                        :: degr(2) = -1
     integer                        :: dim = 0
     real(kind=fp), allocatable     :: coef(:,:,:)
  end type type_polynomial

  type ptr_polynomial
     type(type_polynomial), pointer :: ptr => null()
  end type ptr_polynomial

  type(type_matrix), allocatable    :: transformation_matrices(:)


contains



  subroutine free_polynomial( poly )
    implicit none
    type(type_polynomial), intent(inout) :: poly
    if ( allocated(poly%coef) ) deallocate( poly%coef )
  end subroutine free_polynomial



  subroutine reset_polynomial( &
       poly, &
       nvar, &
       base, &
       degr, &
       dim )
    implicit none
    type(type_polynomial), intent(inout) :: poly
    integer,               intent(in)    :: nvar
    integer,               intent(in)    :: base
    integer,               intent(in)    :: degr(nvar)
    integer,               intent(in)    :: dim
    logical                              :: dealloc
    integer                              :: i

    if ( nvar < 1 .or. nvar > 2 ) STOP 'reset_polynomial : nvar < 1 ou > 2'
    if ( any(degr < 0) )          STOP 'reset_polynomial : degr < 0'
    if ( dim < 0 )                STOP 'reset_polynomial : dim < 0'

    if ( allocated(poly%coef) ) then
       dealloc = .false.
       do i = 1,nvar
          if ( size(poly%coef,i) <= degr(i) ) dealloc = .true.
       end do
       if ( size(poly%coef,nvar+1) < dim ) dealloc = .true.

       if ( dealloc ) deallocate( poly%coef )
    end if

    if ( .not.allocated(poly%coef) ) then
       if ( nvar == 1 ) then
          allocate( poly%coef(degr(1)+1,dim,1) )
       else
          allocate( poly%coef(degr(1)+1,degr(2)+1,dim) )
       end if
    end if

    poly%coef = 0._fp
    poly%degr(1:nvar) = degr
    poly%base = base
    poly%nvar = nvar
    poly%dim = dim

  end subroutine reset_polynomial







  subroutine read_polynomial( &
       poly, &
       filename, &
       nvar, &
       base )
    use mod_util
    implicit none
    type(type_polynomial), intent(inout) :: poly
    character(*),          intent(in)    :: filename
    integer,               intent(in)    :: nvar
    integer, optional,     intent(in)    :: base
    integer                              :: fileunit
    integer                              :: m, n, p, b, i, j, k
    data p /0/

    if ( present(base) ) then
       b = base
    else
       b = 0
    end if

    call get_free_unit( fileunit )
    if ( fileunit == 0 ) STOP "read_polynomial : could not find free unit"

    open( &
         unit = fileunit, &
         file = filename, &
         action = "read" )

    if ( nvar == 1 ) then

       read (unit=fileunit, fmt=*) m, n
       p = 1
       call reset_polynomial( &
            poly, &
            1, &
            b, &
            [m-1], &
            n )

    elseif ( nvar == 2 ) then

       read (unit=fileunit, fmt=*) m, n, p
       call reset_polynomial( &
            poly, &
            2, &
            b, &
            [m-1, n-1], &
            p )

    end if

    do k = 1,p
       do j = 1,n
          do i = 1,m
             read (unit=fileunit, fmt=*) poly%coef(i,j,k)
          end do
       end do
    end do

    close( fileunit )

  end subroutine read_polynomial






  subroutine write_polynomial( &
       poly, &
       filename )    
    use mod_util
    implicit none
    type(type_polynomial), intent(in) :: poly
    character(*),          intent(in) :: filename
    integer                           :: fileunit
    integer                           :: i, j, k

    call get_free_unit( fileunit )
    if ( fileunit == 0 ) STOP "write_polynomial : could not find free unit"

    open( &
         unit = fileunit, &
         file = filename, &
         action = "write" )
    write (unit=fileunit, fmt=*) poly%degr(1:poly%nvar)+1, poly%dim
    if ( poly%nvar == 1 ) then
       do j = 1,poly%dim
          do i = 1,poly%degr(1)+1
             write (fileunit,fmt="(ES22.15)") poly%coef(i,j,1)
          end do
       end do
    elseif ( poly%nvar == 2 ) then
       do k = 1,poly%dim
          do j = 1,poly%degr(2)+1
             do i = 1,poly%degr(1)+1
                write (fileunit,fmt="(ES22.15)") poly%coef(i,j,k)
             end do
          end do
       end do
    end if

    close( fileunit )

  end subroutine write_polynomial








  !subroutine diff( poly, d1, d2 )
  !  implicit none
  !  type(type_polynomial),           intent(in)  :: poly
  !  type(type_polynomial), optional, intent(out) :: d1
  !  type(type_polynomial), optional, intent(out) :: d2
  !  if ( poly%nvar == 1 ) then
  !     if ( present(d1) ) call diff1( poly, d1 )
  !  else
  !     if ( present(d1) .and. present(d2) ) then
  !        call diff2( poly, d1, d2 )
  !     elseif ( present(d1) ) then
  !        call diff2( poly, du=d1 )
  !     else
  !        call diff2( poly, dv=d2 )
  !     end if
  !  end if
  !end subroutine diff




  subroutine diff1( poly, deriv )
    implicit none
    type(type_polynomial), intent(in)  :: poly
    type(type_polynomial), intent(out) :: deriv

    call reset_polynomial( &
         deriv, &
         1, &!poly%nvar, &
         poly%base, &
         [ max( poly%degr(1)-1, 0 ) ], &
         poly%dim )

    select case ( poly%base )
    case (1) ! Chebyshev basis
       call chebdiff( &
            poly%coef(1:poly%degr(1)+1,:,1), &
            deriv%coef(1:max(poly%degr(1),1),:,1), &
            poly%degr(1), &
            poly%dim )

    case (2) ! Bernstein basis
       call berndiff( &
            poly%coef(1:poly%degr(1)+1,:,1), &
            deriv%coef(1:max(poly%degr(1),1),:,1), &
            poly%degr(1), &
            poly%dim )

    case default
       STOP 'diff1 : polynomial in unknown basis'  
    end select

  end subroutine diff1




  subroutine diff2( poly, du, dv )
    implicit none
    type(type_polynomial),           intent(in)  :: poly
    type(type_polynomial), optional, intent(out) :: du
    type(type_polynomial), optional, intent(out) :: dv
    real(kind=fp)                                :: tmp(max(poly%degr(2),1),max(poly%degr(1)+1,1))
    integer                                      :: k

    if ( present(du) ) then
       !! 1st partial derivative
       call reset_polynomial( &
            du, &
            2, &!poly%nvar, &
            poly%base, &
            max( poly%degr - [1,0], 0 ), &
            poly%dim )

       select case ( poly%base )
       case (1) ! Chebyshev basis
          do k = 1,poly%dim
             call chebdiff( &
                  poly%coef(1:poly%degr(1)+1,1:poly%degr(2)+1,k), &
                  du%coef(1:max(poly%degr(1),1),1:poly%degr(2)+1,k), &
                  poly%degr(1), &
                  poly%degr(2)+1 )
          end do

       case (2) ! Bernstein basis
          do k = 1,poly%dim
             call berndiff( &
                  poly%coef(1:poly%degr(1)+1,1:poly%degr(2)+1,k), &
                  du%coef(1:max(poly%degr(1),1),1:poly%degr(2)+1,k), &
                  poly%degr(1), &
                  poly%degr(2)+1 )
          end do

       case default
          STOP 'diff2 : polynomial in unknown basis'   
       end select

    end if


    if ( present(dv) ) then
       !! 2nd partial derivative
       call reset_polynomial( &
            dv, &
            2, &!poly%nvar, &
            poly%base, &
            max( poly%degr - [0,1], 0 ), &
            poly%dim )

       select case ( poly%base )
       case (1) ! Chebyshev basis
          do k = 1,poly%dim
             call chebdiff( &
                  transpose( poly%coef(1:poly%degr(1)+1,1:poly%degr(2)+1,k) ), &
                  tmp, &
                  poly%degr(2), &
                  poly%degr(1)+1 )
             dv%coef(1:size(tmp,2),1:size(tmp,1),k) = transpose(tmp)
          end do

       case (2) ! Bernstein basis
          do k = 1,poly%dim
             call berndiff( &
                  transpose( poly%coef(1:poly%degr(1)+1,1:poly%degr(2)+1,k) ), &
                  tmp, &
                  poly%degr(2), &
                  poly%degr(1)+1 )
             dv%coef(1:size(tmp,2),1:size(tmp,1),k) = transpose(tmp)
          end do

       case default
          STOP 'diff2 : polynomial in unknown basis'       
       end select

    end if

  end subroutine diff2







  subroutine bivar2univar( &
       poly2, &
       poly1, &
       ivar, &
       ival )
    ! Extracts an extremal iso-parametric line of a bivariate polynomial as a univariate polynomial
    implicit none
    type(type_polynomial), intent(in)    :: poly2
    type(type_polynomial), intent(inout) :: poly1
    integer,               intent(in)    :: ivar, ival

    if ( poly2%nvar /= 2 ) STOP 'bivar2univar : poly%nvar /= 2'

    call reset_polynomial( &
         poly=poly1, &
         nvar=1, &
         base=poly2%base, &
         degr=poly2%degr(1+mod(ivar,2)), &
         dim=poly2%dim )

    select case ( poly2%base )
    case (1) ! Chebyshev basis
       call chebbivar2univar( &
            poly2%coef(1:poly2%degr(1)+1,1:poly2%degr(2)+1,1:poly2%dim), &
            poly2%degr(1)+1, &
            poly2%degr(2)+1, &
            poly2%dim, &
            poly1%coef(1:poly1%degr(1)+1,1:poly1%dim,1), &
            ivar, &
            ival )
    case (2) ! Bernstein basis
       call bernbivar2univar( &
            poly2%coef(1:poly2%degr(1)+1,1:poly2%degr(2)+1,1:poly2%dim), &
            poly2%degr(1)+1, &
            poly2%degr(2)+1, &
            poly2%dim, &
            poly1%coef(1:poly1%degr(1)+1,1:poly1%dim,1), &
            ivar, &
            ival )
    case default
       STOP 'bivar2univar : polynomial in unknown basis'       
    end select

  end subroutine bivar2univar










  subroutine polyval1( &
       f, &
       poly, &
       x, &
       n )
    implicit none
    type(type_polynomial), intent(in)  :: poly
    integer,               intent(in)  :: n
    real(kind=fp),         intent(in)  :: x(n)
    real(kind=fp),         intent(out) :: f(n,poly%dim)
    integer                            :: i

    select case ( poly%base )
    case (1) ! Chebyshev basis
       call clenshaw( &
            f, &
            poly%coef(1:poly%degr(1)+1,1:poly%dim,1), &
            x, &
            n, &
            poly%degr(1), &
            poly%dim )

    case (2) ! Bernstein basis
       do i = 1,n
          call de_casteljau( &
               poly%coef(1:poly%degr(1)+1,1:poly%dim,1), &
               x(i), &
               poly%degr(1), &
               poly%dim, &
               f=f(i,1:poly%dim) )
       end do

    case default
       STOP 'polyval1 : polynomial in unknown basis'       
    end select

  end subroutine polyval1







  subroutine polyval2( &
       f, &
       poly, &
       uv, &
       n )
    implicit none
    type(type_polynomial), intent(in)  :: poly
    integer,               intent(in)  :: n
    real(kind=fp),         intent(in)  :: uv(2,n)
    real(kind=fp),         intent(out) :: f(poly%dim,n)
    real(kind=fp)                      :: a(poly%degr(2)+1,poly%dim)
    integer                            :: i, k

    select case ( poly%base )
    case (1) ! Chebyshev basis
       do i = 1,n
          do k = 1,poly%dim
             call clenshaw( &
                  a(:,k), &
                  poly%coef(1:poly%degr(1)+1,1:poly%degr(2)+1,k), &
                  uv(1,i), &
                  1, &
                  poly%degr(1), &
                  poly%degr(2)+1 )
          end do
          call clenshaw( &
               f(:,i), &
               a, &
               uv(2,i), &
               1, &
               poly%degr(2), &
               poly%dim )
       end do

    case (2) ! Bernstein basis
       STOP 'polyval2 : Bernstein polynomials still not supported'       

    case default
       STOP 'polyval2 : polynomial in unknown basis'       
    end select

  end subroutine polyval2








  subroutine subdiv_bezier1( &
       b, &
       t, &
       bl, &
       br )
    implicit none
    type(type_polynomial), intent(in)  :: b
    real(kind=fp),         intent(in)  :: t
    type(type_polynomial), intent(out) :: bl
    type(type_polynomial), intent(out) :: br

    if ( b%base /= 2 ) STOP 'subdiv_bezier1 : polynomial not in Bernstein basis'

    call reset_polynomial( poly=bl, nvar=1, base=2, degr=b%degr(1), dim=b%dim )
    call reset_polynomial( poly=br, nvar=1, base=2, degr=b%degr(1), dim=b%dim )

    call de_casteljau( &
         b%coef(1:b%degr(1)+1,1:b%dim,1), &
         t, &
         b%degr(1), &
         b%dim, &
         bl=bl%coef(1:b%degr(1)+1,1:b%dim,1), &
         br=br%coef(1:b%degr(1)+1,1:b%dim,1) )

  end subroutine subdiv_bezier1





  subroutine subdiv_bezier2( &
       b, &
       uv, &
       bsw, &
       bse, &
       bnw, &
       bne )
    implicit none
    type(type_polynomial), intent(in)                 :: b
    real(kind=fp),         intent(in)                 :: uv(2)
    type(type_polynomial), intent(out), optional      :: bsw
    type(type_polynomial), intent(out), optional      :: bse
    type(type_polynomial), intent(out), optional      :: bnw
    type(type_polynomial), intent(out), optional      :: bne
    real(kind=fp), dimension(b%degr(1)+1,b%degr(2)+1) :: bw, be
    real(kind=fp), dimension(b%degr(2)+1,b%degr(1)+1) :: bsT, bnT
    integer                                           :: m, n, k

    if ( b%base /= 2 ) STOP 'subdiv_bezier1 : polynomial not in Bernstein basis'

    m = b%degr(1)+1
    n = b%degr(2)+1

    if ( present(bsw) ) call reset_polynomial( bsw, 2, 2, b%degr, b%dim )
    if ( present(bse) ) call reset_polynomial( bse, 2, 2, b%degr, b%dim )
    if ( present(bnw) ) call reset_polynomial( bnw, 2, 2, b%degr, b%dim )
    if ( present(bne) ) call reset_polynomial( bne, 2, 2, b%degr, b%dim )

    do k = 1,b%dim
       call de_casteljau( &
            b%coef(1:m,1:n,k), &
            uv(1), &
            b%degr(1), &
            n, &
            bl=bw, &
            br=be )

       call de_casteljau( &
            transpose(bw), &
            uv(2), &
            b%degr(2), &
            m, &
            bl=bsT, &
            br=bnT )
       if ( present(bsw) ) bsw%coef(1:m,1:n,k) = transpose(bsT)
       if ( present(bnw) ) bnw%coef(1:m,1:n,k) = transpose(bnT)

       call de_casteljau( &
            transpose(be), &
            uv(2), &
            b%degr(2), &
            m, &
            bl=bsT, &
            br=bnT )
       if ( present(bse) ) bse%coef(1:m,1:n,k) = transpose(bsT)
       if ( present(bne) ) bne%coef(1:m,1:n,k) = transpose(bnT)     

    end do

  end subroutine subdiv_bezier2




  subroutine subdiv_bezier2_only_u( &
       b, &
       u, &
       bw, &
       be )
    implicit none
    type(type_polynomial), intent(in)  :: b
    real(kind=fp),         intent(in)  :: u
    type(type_polynomial), intent(out) :: bw
    type(type_polynomial), intent(out) :: be
    integer                            :: m, n, k

    m = b%degr(1)+1
    n = b%degr(2)+1

    call reset_polynomial( bw, 2, 2, b%degr, b%dim )
    call reset_polynomial( be, 2, 2, b%degr, b%dim )

    do k = 1,b%dim
       call de_casteljau( &
            b%coef(1:m,1:n,k), &
            u, &
            b%degr(1), &
            n, &
            bl=bw%coef(1:m,1:n,k), &
            br=be%coef(1:m,1:n,k) )
    end do

  end subroutine subdiv_bezier2_only_u





  subroutine subdiv_bezier2_only_v( &
       b, &
       v, &
       bs, &
       bn )
    implicit none
    type(type_polynomial), intent(in)                 :: b
    real(kind=fp),         intent(in)                 :: v
    type(type_polynomial), intent(out)                :: bs
    type(type_polynomial), intent(out)                :: bn
    real(kind=fp), dimension(b%degr(2)+1,b%degr(1)+1) :: bsT, bnT
    integer                                           :: m, n, k

    m = b%degr(1)+1
    n = b%degr(2)+1

    call reset_polynomial( bs, 2, 2, b%degr, b%dim )
    call reset_polynomial( bn, 2, 2, b%degr, b%dim )

    do k = 1,b%dim
       call de_casteljau( &
            transpose( b%coef(1:m,1:n,k) ), &
            v, &
            b%degr(2), &
            m, &
            bl=bsT, &
            br=bnT )
       bs%coef(1:m,1:n,k) = transpose(bsT)
       bn%coef(1:m,1:n,k) = transpose(bnT)
    end do

  end subroutine subdiv_bezier2_only_v






  subroutine cheb2bern_matrix( A, degr )
    ! Transformation matrix from Chebyshev to Bernstein polynomial basis 
    ! "Transformation of Chebyshev–Bernstein polynomial basis", Rababah (2003) -- p.8 (615)
    implicit none
    integer,       intent(in)  :: degr
    real(kind=fp), intent(out) :: A(1:degr+1,1:degr+1)
    real(kind=fp)              :: fnumerator, fdegr, s
    integer                    :: i, j, k, m, n, imin, p

    if ( degr == 0 ) then
       A(1,1) = 1._fp
       return
    end if

    m = degr / 2
    n = degr+1

    A(1:m+1,1) = 1._fp
    A(1,2:n:2) = -1._fp
    A(1,3:n:2) = 1._fp

    A(2:m+1,2:n) = 0._fp
    do k = 1,degr
       fnumerator = factln(2*k) + factln(degr-k)
       do j = 1,m
          imin = max(0,j+k-degr)
          s = real( (-1)**(k-imin), kind=fp )
          do i = imin,min(j,k)
             A(1+j,1+k) = A(1+j,1+k) + s * exp( fnumerator - &
                  (factln(2*(k-i)) + factln(2*i) + factln(degr-k-j+i) + factln(j-i)) &
                  )
             s = -s
          end do
       end do
    end do

    fdegr = factln( degr )
    do j = 1,m
       A(1+j,2:n) = A(1+j,2:n) / exp( fdegr - factln(degr-j) - factln(j) )
    end do

    p = m
    if ( mod(degr,2) /= 0 ) p = p + 1

    A(n:p+1:-1,2:n:2) = -A(1:m+1,2:n:2)
    A(n:p+1:-1,1:n:2) = A(1:m+1,1:n:2)

  end subroutine cheb2bern_matrix




  
  subroutine get_cheb2bern_mat_from_collection( &
       collection, &
       N, &
       c2b )
    implicit none
    integer,            intent(in)    :: N
    type(type_matrix),  intent(inout) :: collection(:)
    real(kind=fp)                     :: c2b(N,N)

    if ( N > size(collection) ) then
       call cheb2bern_matrix( c2b, N-1 )
    else
       if ( .not.allocated(collection(N)%mat) ) then
          allocate( collection(N)%mat(N,N) )
          call cheb2bern_matrix( collection(N)%mat, N-1 )
       end if
       c2b = collection(N)%mat
    end if

  end subroutine get_cheb2bern_mat_from_collection








  subroutine economize1( &
     poly, &
     tol )
    implicit none
    type(type_polynomial), intent(inout) :: poly
    real(kind=fp),         intent(in)    :: tol
    real(kind=fp)                        :: tolsqr
    integer                              :: i, m

    if ( poly%base /= 1 ) STOP 'economize1 : polynomial not in Chebyshev basis'
    if ( poly%nvar /= 1 ) STOP 'economize1 : not a univariate polynomial'

    tolsqr = tol**2

    m = poly%degr(1)+1
    do i = m,1,-1
       if ( sum( poly%coef(i:m,1:poly%dim,1)**2 ) > tolsqr ) exit
    end do
    poly%degr(1) = i-1

  end subroutine economize1

  
  
  
  subroutine economize2( &
       poly, &
       tol )
    implicit none
    type(type_polynomial), intent(inout) :: poly
    real(kind=fp),          intent(in)   :: tol
    real(kind=fp)                        :: tolsqr
    real(kind=fp), dimension(poly%dim)   :: s, r
    integer                              :: i, j

    if ( poly%base /= 1 ) STOP 'economize2 : polynomial not in Chebyshev basis'
    if ( poly%nvar /= 2 ) STOP 'economize2 : not a bivariate polynomial'

    tolsqr = tol**2

    s = sum( abs( poly%coef(1:poly%degr(1)+1,1:poly%degr(2)+1,1:poly%dim) ) )

    do i = poly%degr(1)+1,1,-1
       r = sum( abs( poly%coef(1:i-1,1:poly%degr(2)+1,1:poly%dim) ) )
       if ( sum( (s - r)**2 ) > tolsqr ) exit
    end do
    poly%degr(1) = i-1

    do j = poly%degr(2)+1,1,-1
       r = sum( abs( poly%coef(1:i,1:j-1,1:poly%dim) ) )
       if ( sum( (s - r)**2 ) > tolsqr ) exit
    end do
    poly%degr(2) = j-1

  end subroutine economize2





  
  subroutine chgvar1( &
       polya, &
       polyb, &
       x0, &
       x1 )
    implicit none
    type(type_polynomial), intent(in)    :: polya
    type(type_polynomial), intent(inout) :: polyb    
    real(kind=fp),         intent(in)    :: x0, x1

    if ( polya%base /= 1 ) STOP 'chgvar1 : polynomial not in Chebyshev basis'
    if ( polya%nvar /= 1 ) STOP 'chgvar1 : not a univariate polynomial'
    
    call reset_polynomial( &
         polyb, &
         1, &
         1, &
         polya%degr(1), &
         polya%dim )

    call chgvar( &
         polya%coef(1:polya%degr(1)+1,1:polya%dim,1), &
         polyb%coef(1:polya%degr(1)+1,1:polya%dim,1), &
         x0, &
         x1, &
         polya%degr(1), &
         polya%dim )

  end subroutine chgvar1





  subroutine chgvar2( &
       polya, &
       polyb, &
       xy0, &
       xy1 )
    implicit none
    type(type_polynomial),       intent(in)    :: polya
    type(type_polynomial),       intent(inout) :: polyb    
    real(kind=fp), dimension(2), intent(in)    :: xy0, xy1
    real(kind=fp)                              :: tmp(polya%degr(2)+1,polya%degr(1)+1)
    integer                                    :: k

    if ( polya%base /= 1 ) STOP 'chgvar2 : polynomial not in Chebyshev basis'
    if ( polya%nvar /= 2 ) STOP 'chgvar2 : not a bivariate polynomial'

    call reset_polynomial( &
         poly=polyb, &
         nvar=2, &
         base=1, &
         degr=polya%degr, &
         dim=polya%dim )

    do k = 1,polya%dim
       call chgvar( &
            polya%coef(1:polya%degr(1)+1,1:polya%degr(2)+1,k), &
            polyb%coef(1:polya%degr(1)+1,1:polya%degr(2)+1,k), &
            xy0(1), &
            xy1(1), &
            polya%degr(1), &
            polya%degr(2)+1 )

       call chgvar( &
            transpose( polyb%coef(1:polya%degr(1)+1,1:polya%degr(2)+1,k) ), &
            tmp, &
            xy0(2), &
            xy1(2), &
            polya%degr(2), &
            polya%degr(1)+1 )
       polyb%coef(1:polya%degr(1)+1,1:polya%degr(2)+1,k) = transpose(tmp)
       
    end do

  end subroutine chgvar2












  subroutine cheb2bern( &
       c, &
       b )
    ! Transforms a polynomial from Chebyshev basis into Bernstein basis
    implicit none
    type(type_polynomial), intent(in)    :: c
    type(type_polynomial), intent(inout) :: b
    real(kind=fp)                        :: au(c%degr(1)+1,c%degr(1)+1)
    real(kind=fp)                        :: av(c%degr(2)+1,c%degr(2)+1)
    integer                              :: i

    if ( c%base /= 1 ) STOP 'cheb2bern : input polynomial not in Chebyshev basis'

    call reset_polynomial( poly=b, nvar=c%nvar, base=2, degr=c%degr(1:c%nvar), dim=c%dim )

    if ( .not.allocated(transformation_matrices) ) allocate( transformation_matrices(32) )

    call get_cheb2bern_mat_from_collection( &
         transformation_matrices, &
         c%degr(1)+1, &
         au )

   
    if ( c%nvar == 1 ) then
       ! univariate polynomial
       b%coef(1:c%degr(1)+1,1:b%dim,1) = matmul( au, c%coef(1:c%degr(1)+1,1:c%dim,1) )

    elseif ( c%nvar == 2 ) then
       ! bivariate polynomial
       call get_cheb2bern_mat_from_collection( &
            transformation_matrices, &
            c%degr(2)+1, &
            av )
       av = transpose( av )

       do i = 1,c%dim
          b%coef(1:c%degr(1)+1,1:c%degr(2)+1,i) = &
               matmul( matmul( au, c%coef(1:c%degr(1)+1,1:c%degr(2)+1,i) ), av )
       end do

    else
       ! error
       STOP 'cheb2bern : nvar /= 1,2'

    end if

  end subroutine cheb2bern



end module mod_polynomial
