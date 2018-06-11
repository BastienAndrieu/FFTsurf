program pre_intersection
  
  use mod_math
  use mod_bernstein
  use mod_chebyshev
  use mod_diffgeom
  use mod_types_intersection
  !use mod_intersection
  
  implicit none
 
  logical, parameter                             :: ECONOMIZE = .true.
  type(type_parametric_surface), target, allocatable :: surf(:)
  type(type_surface_region), target, allocatable :: root(:)
  character                                      :: strnum
  type(type_intersection_data)                   :: interdat
  type(ptr_surface_region)                       :: region(2)
  type(ptr_parametric_surface)                   :: surfroot(2)
  type(ptr_bernstein_series2)                    :: bs(2), bpn(2)
  integer                                        :: isurf
  integer*8                                      :: tic, toc, count_rate


  allocate( surf(2) )
  allocate( root(2) )

  do isurf = 1,2
     write (strnum,'(I1)') isurf
     call read_chebyshev_series2( &
          surf(isurf)%s, &
          '/stck/bandrieu/Bureau/coeffstest/C' // strnum // '_test16.txt' )
          !'/home/bastien/Bureau/coeffstest/C' // strnum // '_test07.txt' )
          !'../Bureau/coeffstest/C' // strnum // '_test01.txt' )

     if ( ECONOMIZE ) then
        PRINT *,'     S DEGR =',SURF(ISURF)%S%DEGR
        call economization2( surf(isurf)%s, MATHeps )
        PRINT *,' S_eco DEGR =',SURF(ISURF)%S%DEGR
     end if
     call write_chebyshev_series2( surf(isurf)%s, 'pre_intersection/c_'//strnum//'.cheb' )

     call compute_first_derivatives( surf(isurf) )
     call compute_second_derivatives( surf(isurf) )
     call compute_pseudonormal( surf(isurf) )

     if ( ECONOMIZE ) then
        PRINT *,'    PN DEGR =',SURF(ISURF)%PN%DEGR
        call economization2( surf(isurf)%pn, MATHeps )
        PRINT *,'PN_eco DEGR =',SURF(ISURF)%PN%DEGR
     end if

     call init_surface_region( &
          root(isurf), &
          spread( real( [-1,1], kind=MATHpr ), 2, 2 ) )
  end do

  do isurf = 1,2
     region(isurf)%ptr => root(isurf)
     surfroot(isurf)%ptr => surf(isurf)

     allocate( bs(isurf)%ptr, bpn(isurf)%ptr )
     call cs2bs( surf(isurf)%s, bs(isurf)%ptr )
     call cs2bs( surf(isurf)%pn, bpn(isurf)%ptr )

     write (strnum,'(I1)') isurf
     call write_bernstein_series2( bs(isurf)%ptr, 'pre_intersection/bs_'//strnum//'.bern' )
     call write_bernstein_series2( bpn(isurf)%ptr, 'pre_intersection/bpn_'//strnum//'.bern' )

  end do

  !allocate( interdat%points(100), interdat%curves(100) )
  call system_clock( tic, count_rate )
  call intersect_surface_surface_bezier( &
       surfroot, &
       bs, &
       bpn, &
       region, &
       interdat ) 
  call system_clock( toc )
  PRINT *,''
  PRINT *,'ELAPSED :',real( toc - tic ) / real( count_rate )


  do isurf = 1,size(root)
     write (strnum,'(I1)') isurf
     call export_surface_region_tree( &
          root(isurf), &
          'pre_intersection/tree_' // strnum // '_bezier.dat' )
     call free_surface_region_tree( root(isurf) )
  end do


contains 

  subroutine cs2bs( c, b )
    implicit none
    type(type_chebyshev_series2), intent(in)  :: c
    type(type_bernstein_series2), intent(out) :: b
    real(kind=MATHpr)                         :: au(c%degr(1)+1,c%degr(1)+1)
    real(kind=MATHpr)                         :: av(c%degr(2)+1,c%degr(2)+1)
    integer                                   :: i
    
    call reset_bernstein_series2( b, c%degr, c%dim )

    call ch2be_matrix( au, c%degr(1) )
    call ch2be_matrix( av, c%degr(2) )
    av = transpose( av )

    do i = 1,c%dim
       b%coef(1:c%degr(1)+1,1:c%degr(2)+1,i) = matmul( matmul( au, c%coef(1:c%degr(1)+1,1:c%degr(2)+1,i) ), av )
    end do    
        
  end subroutine cs2bs

end program pre_intersection
