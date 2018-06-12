module mod_bernstein
  use mod_math

  implicit none


  type type_bernstein_series1
     integer                        :: degr = 0
     integer                        :: dim = 0
     real(kind=MATHpr), allocatable :: coef(:,:)
  end type type_bernstein_series1

  type type_bernstein_series2
     integer                        :: degr(2) = 0
     integer                        :: dim = 0
     real(kind=MATHpr), allocatable :: coef(:,:,:)
  end type type_bernstein_series2


  type ptr_bernstein_series1
     type(type_bernstein_series1), pointer :: ptr => null()
  end type ptr_bernstein_series1

  type ptr_bernstein_series2
     type(type_bernstein_series2), pointer :: ptr => null()
  end type ptr_bernstein_series2


contains



subroutine write_bernstein_series2( &
     b, &
     filename )
  use mod_util
  implicit none
  character(*),                 intent(in) :: filename
  type(type_bernstein_series2), intent(in) :: b
  integer                                  :: fileunit, icoef, jcoef, idim

  call get_free_unit( fileunit )
  if ( fileunit == 0 ) STOP "write_bernstein_series2 : could not find free unit"

  open( unit = fileunit, file = filename, action = "write" )
  write (fileunit,*) b%degr+1, b%dim

  do idim = 1,b%dim
     do jcoef = 1,b%degr(2)+1
        do icoef = 1,b%degr(1)+1
           write (fileunit,"(ES22.15)") b%coef(icoef,jcoef,idim)
        end do
     end do
  end do
  close(fileunit)
end subroutine write_bernstein_series2



subroutine de_casteljau( &
     b, &
     t, &
     degr, &
     dim, &
     f, &
     bl, &
     br )
  use mod_math
  implicit none
  integer,           intent(in)            :: degr, dim
  real(kind=MATHpr), intent(in)            :: b(degr+1,dim)
  real(kind=MATHpr), intent(in)            :: t
  real(kind=MATHpr), intent(out), optional :: f(dim)
  real(kind=MATHpr), intent(out), optional :: bl(degr+1,dim)
  real(kind=MATHpr), intent(out), optional :: br(degr+1,dim)
  real(kind=MATHpr), dimension(degr+1,dim) :: bj, bjm1
  integer                                  :: i, j, n
  
  bjm1 = b
  bj(:,:) = 0._MATHpr
  n = degr + 1
  if ( present(bl) ) bl(1,:) = bjm1(1,:)
  if ( present(br) ) br(n,:) = bjm1(n,:)
  
  do j = 2,n
     do i = 1,n-j+1
        bj(i,:) = bjm1(i,:) * (1._MATHpr - t) + bjm1(i+1,:) * t
     end do
     if ( present(bl) ) bl(j,:) = bj(1,:)
     if ( present(br) ) br(n-j+1,:) = bj(n-j+1,:)
     if ( j > degr ) then
        if ( present(f) ) f = bj(1,:)
        return
     end if
     bjm1 = bj
  end do

end subroutine de_casteljau



subroutine subdiv2_along_u( &
     b, &
     u, &
     bw, &
     be )
  implicit none
  type(type_bernstein_series2), intent(in)  :: b
  real(kind=MATHpr),            intent(in)  :: u
  type(type_bernstein_series2), intent(out) :: bw
  type(type_bernstein_series2), intent(out) :: be
  integer                                   :: m, n, k

  m = b%degr(1)+1
  n = b%degr(2)+1
  
  call reset_bernstein_series2( bw, b%degr, b%dim )
  call reset_bernstein_series2( be, b%degr, b%dim )

  do k = 1,b%dim
     call de_casteljau( &
          b%coef(1:m,1:n,k), &
          u, &
          b%degr(1), &
          n, &
          bl=bw%coef(1:m,1:n,k), &
          br=be%coef(1:m,1:n,k) )
  end do

end subroutine subdiv2_along_u



subroutine berndiff( &
     b, &
     d, &
     degr, &
     dim )
  implicit none
  integer,           intent(in)  :: degr, dim
  real(kind=MATHpr), intent(in)  :: b(degr+1,dim)
  real(kind=MATHpr), intent(out) :: d(max(degr,1),dim)
  
  if ( degr < 1 ) then
     d(:,:) = 0._MATHpr
  else
     d = real( degr, kind=MATHpr ) * ( b(2:degr+1,:) - b(1:degr,:) )
  end if

end subroutine berndiff



subroutine subdiv2( &
     b, &
     uv, &
     bsw, &
     bse, &
     bnw, &
     bne )
  use mod_math
  implicit none
  type(type_bernstein_series2), intent(in)              :: b
  real(kind=MATHpr),            intent(in)              :: uv(2)
  type(type_bernstein_series2), intent(out), optional   :: bsw
  type(type_bernstein_series2), intent(out), optional   :: bse
  type(type_bernstein_series2), intent(out), optional   :: bnw
  type(type_bernstein_series2), intent(out), optional   :: bne
  real(kind=MATHpr), dimension(b%degr(1)+1,b%degr(2)+1) :: bw, be
  real(kind=MATHpr), dimension(b%degr(2)+1,b%degr(1)+1) :: bsT, bnT
  integer                                               :: m, n, k

  m = b%degr(1)+1
  n = b%degr(2)+1

  if ( present(bsw) ) call reset_bernstein_series2( bsw, b%degr, b%dim )
  if ( present(bse) ) call reset_bernstein_series2( bse, b%degr, b%dim )
  if ( present(bnw) ) call reset_bernstein_series2( bnw, b%degr, b%dim )
  if ( present(bne) ) call reset_bernstein_series2( bne, b%degr, b%dim )

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

end subroutine subdiv2



subroutine reset_bernstein_series1( &
     b, &
     degr, &
     dim )
  implicit none
  integer,                      intent(in)  :: degr
  integer,                      intent(in)  :: dim
  type(type_bernstein_series1), intent(out) :: b

  if ( allocated(b%coef) ) then
     if ( size(b%coef,1) <= degr .or. &
          size(b%coef,2) < dim ) deallocate(b%coef)
  end if
  if ( .not.allocated(b%coef) ) allocate(b%coef(1:degr+1, dim))

  b%degr = degr
  b%dim = dim
  b%coef(:,:) = 0._MATHpr

end subroutine reset_bernstein_series1



subroutine subdiv1( &
     b, &
     t, &
     bl, &
     br )
  use mod_math
  implicit none
  type(type_bernstein_series1), intent(in)  :: b
  real(kind=MATHpr),            intent(in)  :: t
  type(type_bernstein_series1), intent(out) :: bl
  type(type_bernstein_series1), intent(out) :: br

  call reset_bernstein_series1( bl, b%degr, b%dim )
  call reset_bernstein_series1( br, b%degr, b%dim )

  call de_casteljau( &
       b%coef(1:b%degr+1,:), &
       t, &
       b%degr, &
       b%dim, &
       bl=bl%coef(1:b%degr+1,:), &
       br=br%coef(1:b%degr+1,:) )

end subroutine subdiv1



subroutine berndiff1( &
     b, &
     d )
  implicit none
  type(type_bernstein_series1), intent(in)  :: b
  type(type_bernstein_series1), intent(out) :: d

  call reset_bernstein_series1( d, max( b%degr-1, 0 ), b%dim )
  call berndiff( b%coef(1:b%degr+1,:), d%coef(1:max(b%degr,1),:), b%degr, b%dim )

end subroutine berndiff1



subroutine subdiv2_along_v( &
     b, &
     v, &
     bs, &
     bn )
  implicit none
  type(type_bernstein_series2), intent(in)              :: b
  real(kind=MATHpr),            intent(in)              :: v
  type(type_bernstein_series2), intent(out)             :: bs
  type(type_bernstein_series2), intent(out)             :: bn
  real(kind=MATHpr), dimension(b%degr(2)+1,b%degr(1)+1) :: bsT, bnT
  integer                                               :: m, n, k

  m = b%degr(1)+1
  n = b%degr(2)+1
  
  call reset_bernstein_series2( bs, b%degr, b%dim )
  call reset_bernstein_series2( bn, b%degr, b%dim )

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

end subroutine subdiv2_along_v



subroutine reset_bernstein_series2( &
     b, &
     degr, &
     dim )
  implicit none
  integer,                      intent(in)  :: degr(2)
  integer,                      intent(in)  :: dim
  type(type_bernstein_series2), intent(out) :: b

  if ( allocated(b%coef) ) then
     if ( size(b%coef,1) <= degr(1) .or. &
          size(b%coef,2) <= degr(2) .or. &
          size(b%coef,3) < dim ) deallocate(b%coef)
  end if
  if ( .not.allocated(b%coef) ) allocate( b%coef(1:degr(1)+1, 1:degr(2)+1, dim) )

  b%degr = degr
  b%dim = dim
  b%coef(:,:,:) = 0._MATHpr

end subroutine reset_bernstein_series2



subroutine berndiff2( &
     b, &
     du, &
     dv )
  implicit none
  type(type_bernstein_series2), intent(in)            :: b
  type(type_bernstein_series2), intent(out)           :: du
  type(type_bernstein_series2), intent(out), optional :: dv
  real(kind=MATHpr)                                   :: tmp(max(b%degr(2),1),max(b%degr(1)+1,1))
  integer                                             :: k

  call reset_bernstein_series2( du, max( [b%degr(1)-1, b%degr(2)], 0 ), b%dim )
  do k = 1,b%dim
     call berndiff( &
          b%coef(1:b%degr(1)+1,1:b%degr(2)+1,k), &
          du%coef(1:max(b%degr(1),1),1:b%degr(2)+1,k), &
          b%degr(1), &
          b%degr(2)+1 )
  end do

  if ( present(dv) ) then
     call reset_bernstein_series2( dv, max( [b%degr(1), b%degr(2)-1], 0 ), b%dim )
     do k = 1,b%dim
        call berndiff( &
             transpose( b%coef(1:b%degr(1)+1,1:b%degr(2)+1,k) ), &
             tmp, &
             b%degr(2), &
             b%degr(1)+1 )
        dv%coef(1:size(tmp,2),1:size(tmp,1),k) = transpose(tmp)
     end do
  end if

end subroutine berndiff2



subroutine write_bernstein_series1( &
     b, &
     filename )
  use mod_util
  implicit none
  character(*),                 intent(in) :: filename
  type(type_bernstein_series1), intent(in) :: b
  integer                                  :: fileunit, icoef, idim

  call get_free_unit( fileunit )
  if ( fileunit == 0 ) STOP "write_bernstein_series1 : could not find free unit"

  open( unit = fileunit, file = filename, action = "write" )
  write (fileunit,*) b%degr+1, b%dim

  do idim = 1,b%dim
     do icoef = 1,b%degr+1
        write (fileunit,"(ES22.15)") b%coef(icoef,idim)
     end do
  end do
  close(fileunit)
end subroutine write_bernstein_series1
end module mod_bernstein
