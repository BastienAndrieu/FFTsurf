module mod_errors_intersection

  use mod_util
  implicit none


  integer, parameter :: ERR = 1




  character(200) :: working_directory_errors = ''
  
contains

  subroutine error_intersection_point( &
       msg, &
       surf, &
       uv, &
       xyz, &
       logfile )
    use mod_math
    use mod_diffgeom
    use mod_polynomial
    !...
    implicit none
    character(*),       intent(in)  :: msg
    type(type_surface), intent(in)  :: surf(2)
    real(kind=MATHpr),  intent(in)  :: uv(2,2)
    real(kind=MATHpr),  intent(in)  :: xyz(3)
    character(*),       intent(out) :: logfile
    character(16)                   :: name ! dd_mm_yyyy_xxhxx
    character(200)                  :: filename
    character                       :: strnum
    integer                         :: file_unit, isurf

    ! >>> common to all error-reporting subroutines ------
    call generate_name_from_date( name )
    call system( 'mkdir -p ' // trim(working_directory_errors) // name )
    call get_free_unit( file_unit )
    logfile = trim(working_directory_errors) // name // '/error.log'
    open( &
         unit = file_unit, &
         file = logfile, &
         action = 'write' )
    ! ----------------------- <<<

    write ( file_unit, '(A10,1x,A2,1x,A5)' ) name(1:2),'/',name(4:5),'/',name(7:10),'at',name(12:16)
    write ( file_unit, * ) 'Error message:'
    write ( file_unit, * ) trim( msg )
    write ( file_unit, * ) 

    write ( file_unit, * ) 'Intersection point:'
    write ( file_unit, * ) 'uv  ='
    write ( file_unit, '(ES22.15)' ) uv
    write ( file_unit, * ) 'xyz ='
    write ( file_unit, '(ES22.15)' ) xyz
    write ( file_unit, * ) 
    write ( file_unit, * ) 'Surface Chebyshev coefficients:'
    do isurf = 1,2
       write ( strnum, '(I1)' ) isurf
       filename = trim(working_directory_errors) // name // '/surface_' // strnum // '.cheb'
       call write_polynomial(surf(isurf)%x, trim(filename))
       write ( file_unit, * ) trim(filename)
    end do

    close( file_unit )

  end subroutine error_intersection_point

end module mod_errors_intersection
