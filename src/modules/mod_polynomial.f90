module mod_polynomial

  use mod_constants
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
    implicit none
    type(type_polynomial), intent(in)    :: poly2
    type(type_polynomial), intent(inout) :: poly1
    integer,               intent(in)    :: ivar, ival

    if ( poly2%nvar /= 2 ) STOP 'bivar2univar : poly%nvar /= 2'

    call reset_polynomial( &
         poly1, &
         1, &
         poly2%base, &
         poly2%degr(1+mod(ivar,2)), &
         poly2%dim )

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













end module mod_polynomial
