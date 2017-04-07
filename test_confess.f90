program tests
   use confess
   implicit none

   integer :: i, j
   real(kind=rp) :: A(nqd,nen), B(nqd,nen,nen)

   do j = 1,nen
      write(*,*) (N(i,j), i = 1,nqd)
   end do
   write(*,*)
!   A = grad(N)
!   do j = 1,nen
!      write(*,*) (A(i,j), i = 1,nqd)
!   end do
!   write(*,*)
!   B = N.otimes.N
!   write(*,*) B
end program tests
