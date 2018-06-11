program pre_intersection
  
  use mod_math
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
  type(ptr_chebyshev_series2)                    :: surfc(2), surfpn(2)
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

     call write_chebyshev_series2( surf(isurf)%pn, 'pre_intersection/pn_'//strnum//'.cheb' )
     if ( ECONOMIZE ) then
        PRINT *,'    PN DEGR =',SURF(ISURF)%PN%DEGR
        call economization2( surf(isurf)%pn, MATHeps )
        PRINT *,'PN_eco DEGR =',SURF(ISURF)%PN%DEGR
     end if
     call write_chebyshev_series2( surf(isurf)%pn, 'pre_intersection/pn_eco_'//strnum//'.cheb' )

     call init_surface_region( &
          root(isurf), &
          spread( real( [-1,1], kind=MATHpr ), 2, 2 ) )
  end do

  do isurf = 1,2
     region(isurf)%ptr => root(isurf)
     surfroot(isurf)%ptr => surf(isurf)
     surfc(isurf)%ptr => surf(isurf)%s
     surfpn(isurf)%ptr => surf(isurf)%pn
  end do

  !allocate( interdat%points(100), interdat%curves(100) )
  call system_clock( tic, count_rate )
  call intersect_surface_surface( &
       surfroot, &
       surfc, &
       surfpn, &
       ![surf(1)%s, surf(2)%s], &
       ![surf(1)%pn, surf(2)%pn], &
       region, &
       interdat ) 
  !CALL TEST( &
  !     SURFROOT, &!SURF(1:2), &
  !     SURFC([1,2]), & ![SURF(1)%S, SURF(2)%S], &
  !     REGION, &
  !     INTERDAT )
  call system_clock( toc )
  PRINT *,''
  PRINT *,'ELAPSED :',real( toc - tic ) / real( count_rate )


  do isurf = 1,size(root)
     write (strnum,'(I1)') isurf
     call export_surface_region_tree( &
          root(isurf), &
          'pre_intersection/tree_' // strnum // '.dat' )
     call free_surface_region_tree( root(isurf) )

     !deallocate( surf(isurf)%s%coef, &
     !     surf(isurf)%su%coef, &
     !     surf(isurf)%sv%coef, &
     !     surf(isurf)%suu%coef, &
     !     surf(isurf)%suv%coef, &
     !     surf(isurf)%svv%coef, &
     !     surf(isurf)%pn%coef )
  end do
 
  !deallocate( interdat%points, interdat%curves )


CONTAINS
  
  SUBROUTINE TEST( &
       SURF, &
       SURFC, &
       REGION, &
       INTERDAT )
    IMPLICIT NONE
    TYPE(PTR_PARAMETRIC_SURFACE), INTENT(IN)   :: SURF(2)
    !TYPE(TYPE_PARAMETRIC_SURFACE), INTENT(IN)   :: SURF(2)
    TYPE(PTR_CHEBYSHEV_SERIES2), INTENT(IN)    :: SURFC(2)
    !TYPE(TYPE_CHEBYSHEV_SERIES2), INTENT(IN)    :: SURFC(2)
    TYPE(PTR_SURFACE_REGION),     INTENT(INOUT) :: REGION(2)
    TYPE(TYPE_INTERSECTION_DATA), INTENT(INOUT) :: INTERDAT


    !PRINT *,KIND(SURFC%COEF) * SIZE(SURFC%COEF) ! + KIND(SURF%DEGR)*SIZE(SURFC%DEGR) + KIND(SURF%DIM)
    PRINT *,SIZEOF(SURFC)
    !PRINT *,SIZEOF(SURFC%COEF) + SIZEOF(SURFC%DEGR) + SIZEOF(SURFC%DIM)

    PRINT *,SIZEOF(SURF)

    IF ( ASSOCIATED(REGION(1)%PTR%CHILD) ) PRINT *,'HAS CHILDREN'
    IF ( INTERDAT%NP > 0 ) PRINT *,'NP > 0'

  END SUBROUTINE TEST

end program pre_intersection
