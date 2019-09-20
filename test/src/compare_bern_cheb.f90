program compare_bern_cheb

use mod_util
use mod_math
use mod_chebyshev
use mod_bernstein
use mod_polynomial
use mod_obb

implicit none

integer, parameter :: nrep = 1000, ntest = 19
integer                      :: itest, isurf, jtest, irep, ibox
character(2)                 :: stritest
character                    :: strisurf, stribox
type(type_polynomial)        :: c, b
type(type_obb)               :: box(2)
real*8, dimension(2,ntest*2) :: timer, volume
integer*8                    :: count_rate, start, finish
integer                      :: fid

!====== OBBs ======!
! comparer volumes et temps CPU
jtest = 0
do itest = 1,ntest
    write (stritest,'(i2.2)') itest
    do isurf = 1,2
        jtest = jtest + 1
        PRINT *, ''
        PRINT *, jtest, ') TEST #', itest, ', SURF #', isurf

        write (strisurf,'(i1)') isurf
        call read_polynomial( &
            c, &
            'coeffstest/C'//strisurf//'_test'//stritest//'.txt', &
            nvar=2, &
            base=1 &
        )
        call economize2( &
            c, &
            epsilon(1._fp) )

        PRINT *,'    degree = ', c%degr
        call cheb2bern_poly(c, b)

        call system_clock(start, count_rate)
        do irep = 1,nrep
            call chebOBB2( &
                c%coef(1:c%degr(1)+1,1:c%degr(2)+1,1:3), &
                c%degr, &
                box(1) &
            )
        end do
        call system_clock(finish)
        timer(1,jtest) = dble(finish - start)/dble(nrep*count_rate)


        call system_clock(start, count_rate)
        do irep = 1,nrep
            call bernOBB2( &
                b%coef(1:b%degr(1)+1,1:b%degr(2)+1,1:3), &
                b%degr, &
                box(2) &
            )
        end do
        call system_clock(finish)
        timer(2,jtest) = dble(finish - start)/dble(nrep*count_rate)


        do ibox = 1,2
            volume(ibox,jtest) = product(box(ibox)%rng)
            write (stribox, '(i1)') ibox
            call write_obb( &
                box(ibox), &
                'compare_bern_cheb/test'//stritest//'surf'//strisurf//'box'//stribox//'.dat' &
            )
        end do

        PRINT *,'    C volume =', volume(1,jtest)
        PRINT *,'    B volume =', volume(2,jtest)
        PRINT *,'    B/C volume ratio =', volume(2,jtest)/volume(1,jtest)


    end do
end do


call get_free_unit(fid)
open(unit=fid, file='compare_bern_cheb/obb2.dat', action='write')
do jtest = 1,ntest*2
    write (fid,*) timer(1:2,jtest), volume(1:2,jtest)
end do
close(fid)




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


end program compare_bern_cheb