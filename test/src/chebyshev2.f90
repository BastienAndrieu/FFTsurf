program chebyshev2

  use mod_util
  use mod_constants
  use mod_chebyshev2
  use mod_polynomial

  implicit none
  
  integer, parameter         :: nv = 1000
  type(type_polynomial)      :: c1(2), d1(2), c2(2), du2(2), dv2(2)
  real(kind=fp), allocatable :: f1cgl(:,:), f1(:,:), t(:), f(:,:)
  real(kind=fp), allocatable :: f2cgl(:,:,:), f2ij(:,:,:), f2(:,:), uv(:,:)
  integer                    :: m, n, dim, p
  integer                    :: fileunit
  integer                    :: i, j, k,ival, ivar


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!        UNIVARIATE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call read_polynomial( c1(1), "chebyshev2/c1.cheb", nvar=1, base=1 )

  call get_free_unit( fileunit )
  open( unit=fileunit, file="chebyshev2/f1cgl.dat", action="read" )
  read (fileunit,*) m, dim
  allocate( f1cgl(m,dim) )
  do i = 1,m
     read (fileunit,*) f1cgl(i,1:dim)
  end do
  close(fileunit)


  print *,'fcht1'
  call reset_polynomial( poly=c1(2), nvar=1, base=1, degr=[m-1], dim=dim )
  call fcht1( &
       f1cgl, &
       c1(2)%coef, &
       m, &
       dim, &
       epsilon=0._fp )
  print *, maxval( abs(c1(1)%coef - c1(2)%coef) )


  print *,'ifcht1'
  allocate( f(m,dim) )
  call ifcht1( &
       c1(1)%coef(1:m,1:dim,1), &
       f, &
       m, &
       dim )
  print *, maxval( abs(f1cgl - f) )
  deallocate( f1cgl, f )


  print *,'diff1'
  call read_polynomial( d1(1), "chebyshev2/d1.cheb", nvar=1, base=1 )
  call diff1( c1(1), d1(2) )
  print *, maxval( abs(d1(1)%coef - d1(2)%coef) )


  print *,'polyval1'
  call get_free_unit( fileunit )
  open( unit=fileunit, file="chebyshev2/tf1.dat", action="read" )
  read (fileunit,*) p, dim
  allocate( t(p), f1(p,dim) )
  do i = 1,p
     read (fileunit,*) t(i), f1(i,1:dim)
  end do
  close(fileunit)
  
  allocate( f(p,dim) )
  call polyval1( &
       f, &
       c1(1), &
       t, &
       p )
  print *, sqrt( maxval( sum( (f1 - f)**2, 2 ) ) )
  deallocate( t, f, f1 )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!        BIVARIATE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call read_polynomial( c2(1), "chebyshev2/c2.cheb", nvar=2, base=1 )

  call get_free_unit( fileunit )
  open( unit=fileunit, file="chebyshev2/f2cgl.dat", action="read" )
  read (fileunit,*) m, n, dim
  allocate( f2cgl(m,n,dim) )
  do k = 1,dim
     do j = 1,n
        do i = 1,m
           read (fileunit,*) f2cgl(i,j,k)
        end do
     end do
  end do
  close(fileunit)


  print *,'fcht2'
  call reset_polynomial( poly=c2(2), nvar=2, base=1, degr=[m-1,n-1], dim=dim )
  call fcht2( &
       f2cgl, &
       c2(2)%coef, &
       m, &
       n, &
       dim, &
       epsilon=0._fp )
  print *, maxval( abs(c2(1)%coef - c2(2)%coef) )


  print *,'ifcht1'
  allocate( f2ij(m,n,dim) )
  call ifcht2( &
       c2(1)%coef, &
       f2ij, &
       m, &
       n, &
       dim )
  print *, maxval( abs(f2cgl - f2ij) )
  deallocate( f2cgl, f2ij )


  print *,'diff2'
  call read_polynomial( du2(1), "chebyshev2/du2.cheb", nvar=2, base=1 )
  call read_polynomial( dv2(1), "chebyshev2/dv2.cheb", nvar=2, base=1 )

  call diff2( c2(1), du=du2(2), dv=dv2(2) )
  print *, maxval( abs(du2(1)%coef - du2(2)%coef) )
  print *, maxval( abs(dv2(1)%coef - dv2(2)%coef) )



  print *,'polyval2'
  call get_free_unit( fileunit )
  open( unit=fileunit, file="chebyshev2/uvf2.dat", action="read" )
  read (fileunit,*) p, dim
  allocate( uv(2,p), f2(dim,p) )
  do i = 1,p
     read (fileunit,*) uv(:,i), f2(:,i)
  end do
  close(fileunit)

  allocate( f(dim,p) )
  call polyval2( &
       f, &
       c2(1), &
       uv, &
       p )
  print *, sqrt( maxval( sum( (f2 - f)**2, 1 ) ) )
  deallocate( f, uv )




  print *,'biv2univ'
  allocate( t(nv), uv(2,nv), f1(nv,c2(1)%dim), f(c2(1)%dim,nv) )
  t = -1._fp + 2._fp * [( real(i-1,kind=fp)/real(nv-1,kind=fp) , i=0,nv-1 )]
  do ivar = 1,2
     uv(1+mod(ivar,2),:) = t
     do ival = 1,2
        uv(ivar,:) = real( (-1)**ival, kind=fp )
        call bivar2univar( &
             c2(1), &
             c1(1), &
             ivar, &
             ival )

        call polyval1( &
             f1, &
             c1(1), &
             t, &
             nv )
        call polyval2( &
             f, &
             c2(1), &
             uv, &
             nv )
        print *,sqrt( maxval( sum( (transpose(f1) - f)**2, 1 ) ) )

     end do
  end do




end program chebyshev2
