program subdiv2_bezier

  use mod_util
  use mod_math
  use mod_bernstein
  use mod_chebyshev
  use mod_obb

  implicit none
  
  type(type_chebyshev_series2) :: c
  type(type_bernstein_series2) :: b, bsw, bse, bnw, bne
  real(kind=MATHpr)            :: uv(2)
  type(type_obb)               :: box
  
  call read_chebyshev_series2( &
       c, &
       '/stck/bandrieu/Bureau/coeffstest/C2_test18.txt' )
  call write_chebyshev_series2( c, 'subdiv2_bezier/c.cheb' )
     

  call cs2bs( c, b )

  call write_bernstein_series2( b, 'subdiv2_bezier/b.bern' )
       
  !call random_number( uv )
  uv(:) = 0.5_MATHpr
  call subdiv2( &
       b, &
       uv, &
       bsw, &
       bse, &
       bnw, &
       bne )
  
  call write_bernstein_series2( bsw, 'subdiv2_bezier/bsw.bern' )
  call write_bernstein_series2( bse, 'subdiv2_bezier/bse.bern' )
  call write_bernstein_series2( bnw, 'subdiv2_bezier/bnw.bern' )
  call write_bernstein_series2( bne, 'subdiv2_bezier/bne.bern' )

  
  call bernOBB2( b%coef(1:b%degr(1)+1,1:b%degr(2)+1,:), b%degr, box )
  call write_obb( box, 'subdiv2_bezier/obb.dat' )

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

end program subdiv2_bezier
