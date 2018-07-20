module mod_bernstein

  use mod_math

  implicit none

contains

  subroutine de_casteljau( &
       b, &
       t, &
       degr, &
       dim, &
       f, &
       bl, &
       br )
    implicit none
    integer,       intent(in)            :: degr, dim
    real(kind=fp), intent(in)            :: b(degr+1,dim)
    real(kind=fp), intent(in)            :: t
    real(kind=fp), intent(out), optional :: f(dim)
    real(kind=fp), intent(out), optional :: bl(degr+1,dim)
    real(kind=fp), intent(out), optional :: br(degr+1,dim)
    real(kind=fp), dimension(degr+1,dim) :: bj, bjm1
    integer                                  :: i, j, n

    bjm1 = b
    bj(:,:) = 0._fp
    n = degr + 1
    if ( present(bl) ) bl(1,:) = bjm1(1,:)
    if ( present(br) ) br(n,:) = bjm1(n,:)

    do j = 2,n
       do i = 1,n-j+1
          bj(i,:) = bjm1(i,:) * (1._fp - t) + bjm1(i+1,:) * t
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








  subroutine berndiff( &
       b, &
       d, &
       degr, &
       dim )
    implicit none
    integer,       intent(in)  :: degr, dim
    real(kind=fp), intent(in)  :: b(degr+1,dim)
    real(kind=fp), intent(out) :: d(max(degr,1),dim)

    if ( degr < 1 ) then
       d(:,:) = 0._fp
    else
       d = real( degr, kind=fp ) * ( b(2:degr+1,:) - b(1:degr,:) )
    end if

  end subroutine berndiff








  subroutine bernbivar2univar( &
       b2, &
       m, &
       n, &
       p, &
       b1, &
       ivar, &
       ival )
    implicit none
    integer,       intent(in)  :: m, n, p
    real(kind=fp), intent(in)  :: b2(m,n,p)
    integer,       intent(in)  :: ivar, ival
    real(kind=fp), intent(out) :: b1(size(b2,1+mod(ivar,2)),p)
    integer                    :: j

    j = 1 + (ival - 1)*size(b2,ivar)
    if ( ivar == 1 ) then
       b1 = b2(j,:,:)
    else
       b1 = b2(:,j,:)
    end if

  end subroutine bernbivar2univar




end module mod_bernstein
