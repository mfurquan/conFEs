module deform
!*********************************************************************
! Definitionn of Form Functions
! Stiffness matrix for Laplace equation
!---------------------------------------------------------------------
! Apr 18, 2017 | M. Furquan
!*********************************************************************
   use inform
   implicit none

   contains

      pure function K_laplac(inod,jnod,x_ele,u_ele)
         implicit none
         integer,     intent(in)  :: inod, jnod
         type(vector),intent(in)  :: x_ele(nen)
         real(kind=rp),intent(in) :: u_ele(ndf,nen)
         real(kind=rp)            :: K_laplace(ndf)

         K_laplace = integrate(grad_N(:,inod).dot.grad_N(:,jnod),x_ele)
      end function K_laplac
end module deform
