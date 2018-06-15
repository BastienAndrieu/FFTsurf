program tree
  
  use mod_util
  use mod_math
  use mod_regiontree
  use mod_obb

  implicit none

  integer                    :: fileunit, dim, npts
  real(kind=fp), allocatable :: points(:,:), uvboxroot(:)
  type(type_region)          :: root
  integer*8                  :: tic, toc, countrate
  integer                    :: ipt, idim

  call get_free_unit( fileunit )

  open( &
       unit=fileunit, &
       file='tree/points.dat', &
       action='read' )
  read (fileunit,*) npts, dim
  allocate( points(npts,dim), uvboxroot(2*dim) )
  do ipt = 1,npts
     read (fileunit,*) points(ipt,:)
  end do
  close( fileunit )


  do idim = 1,dim
     uvboxroot(2*idim-1) = minval( points(:,idim) )
     uvboxroot(2*idim)   = maxval( points(:,idim) )
  end do
  PRINT *,'UVBOXROOT =',UVBOXROOT

  call init_region( &
       root, &
       dim, &
       uvboxroot )
  deallocate( uvboxroot )
  
  call append_n( &
       root%ipts, &
       root%npts, &
       [ (ipt,ipt=1,npts) ], &
       npts )

  call system_clock( tic, countrate )
  call test_tree( &
       points, &
       npts, &
       dim, &
       root )
  call system_clock( toc )
  print *,'ELAPSED =',real(toc - tic) / real(countrate)
  

  call export_region_tree( &
       root, &
       'tree/tree.dat' )

  call free_region_tree( root )
  deallocate( points )


contains



  recursive subroutine test_tree( &
       points, &
       npts, &
       dim, &
       region )
    ! Test subroutine : given a set of points in R^dim, constructs a region tree
    ! such that each leaf region contains no more than one point
    use mod_util
    implicit none
    logical, parameter               :: subdivide_at_median = .true.
    integer,           intent(in)    :: npts, dim
    real(kind=fp),     intent(in)    :: points(npts,dim)
    type(type_region), intent(inout) :: region
    logical, allocatable             :: mask(:)
    real(kind=fp), allocatable       :: x(:)
    real(kind=fp)                    :: uvm(dim)
    integer                          :: stat
    integer                          :: jpt, ipt, idim, ichild

    if ( associated(region%parent) ) then
       if ( region%parent%npts > 0 ) then
          allocate( mask(region%parent%npts) )
          mask(1:region%parent%npts) = .true.
          outer : do jpt = 1,region%parent%npts
             ipt = region%parent%ipts(jpt)
             do idim = 1,dim
                if ( points(ipt,idim) < region%uvbox(2*idim-1) .or. &
                     points(ipt,idim) > region%uvbox(2*idim) ) then
                   mask(jpt) = .false.
                   cycle outer
                end if
             end do
          end do outer

          call append_n( &
               region%ipts, &
               region%npts, &
               pack( region%parent%ipts(1:region%parent%npts), mask ), &
               count(mask), &
               unique=.true. )
          
          deallocate( mask )
       end if
    end if

    if ( region%npts < 2 ) return

    if ( subdivide_at_median ) then
       ! pick median point along each dimension
       allocate( x(region%npts) )
       do idim = 1,dim
          x = points( region%ipts(1:region%npts), idim )
          call bubblesort_double( x )
          if ( mod(region%npts,2) == 0 ) then
             uvm(idim) = 0.5_fp * sum( x(region%npts/2+[0,1]) )
          else
             uvm(idim) =  x(region%npts/2+1)
          end if
       end do
       deallocate( x )
    else
       uvm = 0.5_fp * ( region%uvbox(1::2) + region%uvbox(2::2) )
    end if


    
    call subdiv_region( &
         region, &
         uvm, &
         stat )

    if ( stat >= 0 .and. stat < dim ) then
       do ichild = 1,size(region%child)
          call test_tree( &
               points, &
               npts, &
               dim, &
               region%child(ichild) )
       end do
    end if

  end subroutine test_tree




  subroutine print_region( region )
    implicit none
    type(type_region), intent(in) :: region
    
    print *,'-------------------------------------'
    print *,'dim =',region%dim
    if ( allocated( region%uvbox ) ) then
       print *,'uvbox =',region%uvbox
    else
       print *,'uvbox = N/A'
    end if
    if ( associated( region%child ) ) then
       print *,size(region%child),' children'
    else
       print *,0,' children'
    end if
    print *,'-------------------------------------'

  end subroutine print_region

  
end program tree
