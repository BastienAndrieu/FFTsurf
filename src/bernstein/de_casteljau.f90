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
