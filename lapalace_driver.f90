module deform
!*********************************************************************
! Definitionn of Form Functions
! Stiffness matrix for Laplace equation
!---------------------------------------------------------------------
! Apr 18, 2017 | M. Furquan
!*********************************************************************
   use inform
   use assembler
   implicit none

   integer :: nn, ne
   real(kind=rp),allocatable :: x(:,:), K(:,:), d(:), f(:), u(:,:)
   integer,      allocatable :: conn(:,:), id(:,:)
   real(kind=rp)             :: x_e(nsd,nen), u_e(ndf,nen)
   integer                   :: id_e(ndf,nen)

   do concurrent (iele = 1:ne)
      do concurrent(ien = 1:nen)
         x_e (:,ien) = x (:,conn(ien,iele))
         u_e (:,ien) = u (:,conn(ien,iele))
         id_e(:,ien) = id(:,conn(ien,iele))
      end do
      call fill_matrix(K_laplace,x_e,u_e,id_e,K)
      call fill_vector()
   end do

   contains

      pure function K_laplac(inod,jnod,xe,ue)
         implicit none
         integer,     intent(in)  :: inod, jnod
         real(kind=rp),intent(in) :: xe(nsd,nen), ue(ndf,nen)
         real(kind=rp)            :: K_laplace(ndf)
         type(vector)             :: x_ele(nen)
         integer                  :: i

         do concurrent (i = 1:nen)
            x_ele(i)%e = xe(:,i)
         end do
         K_laplace = integrate(grad_N(:,inod).dot.grad_N(:,jnod),x_ele)
      end function K_laplac
end module deform
