program bernstein2

  use mod_math
  use mod_bernstein
  use mod_polynomial

  implicit none
  
  type(type_polynomial) :: c
  type(type_polynomial) :: b, bw, be, bs, bn, bsw, bse, bnw, bne
  real(kind=fp)         :: uv(2)

  call read_polynomial( &
       c, &
       '/stck/bandrieu/Bureau/coeffstest/C2_test18.txt', &
       nvar=2, &
       base=1 )
  call write_polynomial( c, 'bernstein2/c.cheb' )


  call cheb2bern_poly( c, b )


  call write_polynomial( b, 'bernstein2/b.bern' )
  
  !call random_number( uv )
  !uv(:) = 0.5_fp
  uv = [ 0.34_fp, 0.87_fp ]
  call subdiv_bezier2( &
       b, &
       uv, &
       bsw, &
       bse, &
       bnw, &
       bne )

  call subdiv_bezier2_only_u( &
       b, &
       uv(1), &
       bw, &
       be )

  call subdiv_bezier2_only_v( &
       b, &
       uv(2), &
       bs, &
       bn )

  call write_polynomial( bw, 'bernstein2/bw.bern' )
  call write_polynomial( be, 'bernstein2/be.bern' )
  call write_polynomial( bs, 'bernstein2/bs.bern' )
  call write_polynomial( bn, 'bernstein2/bn.bern' )
  call write_polynomial( bsw, 'bernstein2/bsw.bern' )
  call write_polynomial( bse, 'bernstein2/bse.bern' )
  call write_polynomial( bnw, 'bernstein2/bnw.bern' )
  call write_polynomial( bne, 'bernstein2/bne.bern' )
  
contains
  
  subroutine cheb2bern_poly( c, b )
    implicit none
    type(type_polynomial), intent(in)  :: c
    type(type_polynomial), intent(out) :: b
    real(kind=fp)                      :: au(c%degr(1)+1,c%degr(1)+1)
    real(kind=fp)                      :: av(c%degr(2)+1,c%degr(2)+1)
    integer                            :: i
    
    if ( c%base /= 1 ) STOP 'cheb2bern_poly : input polynomial not in Chebyshev basis'

    call reset_polynomial( poly=b, nvar=c%nvar, base=2, degr=c%degr, dim=c%dim )
    
    call cheb2bern_matrix( au, c%degr(1) )
    
    select case (c%nvar)
    case (1)
       b%coef(1:c%degr(1)+1,1:b%dim,1) = matmul( au, c%coef(1:c%degr(1)+1,1:c%dim,1) )

    case (2)
       call cheb2bern_matrix( av, c%degr(2) )
       av = transpose( av )

       do i = 1,c%dim
          b%coef(1:c%degr(1)+1,1:c%degr(2)+1,i) = matmul( &
               matmul( au, c%coef(1:c%degr(1)+1,1:c%degr(2)+1,i) ), &
               av )
       end do
       
    case default
       STOP 'cheb2bern_poly : nvar /= 1,2'

    end select

  end subroutine cheb2bern_poly

end program bernstein2
