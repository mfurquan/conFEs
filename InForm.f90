module inform
!**********************************************************************
! Ingredients for defining Form Functions
! M. Furquan | Apr 07, 2017
!**********************************************************************
   use control
   use gauss_legendre
   use linalg
   implicit none

   integer,private :: i, j, k, p, q, r

   logical,parameter,dimension(nqd1,nen1,nen1) ::                     &
      mask = spread(reshape([((i/=j, i = 1,nen1), j = 1,nen1)],       &
                            [nen1,nen1]),1,nqd1)
   
   real(kind=rp),parameter,dimension(nqd) ::                          &
      wq   = [(((wg(qi+i)*wg(qi+j)*wg(qi+k),                          &
               i = 0,nqd1-1), j = 0,nqd1-1), k = 0,nqd1-1)]

   real(kind=rp),parameter,dimension(nen1) ::                         &
      xm   = [(-1+2*i/deg, i=0,deg)]                 

   real(kind=rp),parameter,dimension(nqd1,nen1,nen1) ::               &
      xq   = spread(spread(xg(qi:qf),2,nen1),3,nen1),                 &
      xa   = spread(spread(xm,1,nen1),1,nqd1),                        &
      xb   = spread(spread(xm,1,nqd1),3,nen1)        

   real(kind=rp),parameter,dimension(nqd1,nen1) ::                    &
      N1   = product(xq - xa,3,mask)/product(xb - xa,3,mask),         &
      G1   = sum(1./(xq - xa),3,mask)*N1             

   real(kind=rp),parameter,dimension(nqd,nen) ::                      &
      N    = reshape([((((((N1(i,p)*N1(j,q)*N1(k,r),                  &
                     i = 1,nqd1), j = 1,nqd1), k = 1,nqd1),           &
                     p = 1,nen1), q = 1,nen1), r = 1,nen1)],          &
                     [nqd,nen])

   type(vector),parameter,dimension(nqd,nen) ::                       &
      N_xi = reshape([((((((vector[                                   &
                     (G1(i,p)*N1(j,q)*N1(k,r),                        &
                      N1(i,p)*G1(j,q)*N1(k,r),                        &
                      N1(i,p)*N1(j,q)*G1(k,r)]),                      &
                      i = 1,nqd1), j = 1,nqd1), k = 1,nqd1),          &
                      p = 1,nen1), q = 1,nen1), r = 1,nen1)],         &
                      [nqd,nen])

   interface integrate
      module procedure integrate_scalar, integrate_vector,            &
                       integrate_tensor
   end interface integrate

   contains

      elemental function integrate_scalar(expr,x_ele)
         implicit none
         real(kind=rp),intent(in) :: expr(nqd)
         type(vector), intent(in) :: x_ele(nen)
         real(kind=rp)            :: integrate_scalar

         integrate_scalar =  SUM(wq*jac(x_ele)*expr)
      end function integrate_scalar

      elemental function integrate_vector(expr,x_ele)
         implicit none
         type(vector),intent(in)  :: expr(nqd)
         type(vector),intent(in)  :: x_ele(nen)
         type(vector)             :: integrate_vector
         real(kind=rp)            :: jacob(nqd)
         integer                  :: iq

         jacob = jac(x_ele)
         integrate_vector%e =  0._rp
         do concurrent (iq =  1:nqd)
            integrate_vector = integrate_vector
                             + wq(iq)*expr(iq)*jacob(iq)
         end do
      end function integrate_vector

      elemental function integrate_tensor(expr,x_ele)
         implicit none
         type(tensor),intent(in)  :: expr(nqd)
         type(vector),intent(in)  :: x_ele(nen)
         type(tensor)             :: integrate_tensor
         real(kind=rp)            :: jacob(nqd)
         integer                  :: iq

         jacob = jac(x_ele)
         integrate_tensor%e =  0._rp
         do concurrent (iq =  1:nqd)
            integrate_tensor = integrate_tensor
                             + wq(iq)*expr(iq)*jacob(iq)
         end do
      end function integrate_tensor

      pure function x_xi(x_ele)
         implicit none
         type(vector),intent(in)  :: x_ele(nen)
         type(tensor)             :: x_xi(nqd)
         integer                  :: i, j

         x_xi = tensor(reshape([0.,0.,0.,                             &
                                0.,0.,0.,                             &
                                0.,0.,0.],[3,3])
         do concurrent (i = 1:nqd, j = 1:nen)
            x_xi(i) = x_xi(i) + N_xi(i,j).otimes.x_ele(j)
         end do
      end function x_xi

      pure function jac(x_ele)
         implicit none
         type(vector),intent(in)  :: x_ele(nen)
         real(kind=rp)            :: jac(nqd)
         type(tensor)             :: jm(nqd)
         integer                  :: i

         jac = det(x_xi(x_ele))
      end function jac

      pure function grad_N(x_ele)
         implicit none
         type(vector),intent(in)  :: x_ele(nen)
         type(vector),intent(in)  :: grad_N(nqd,nen)
         type(tensor)             :: tmp(nqd)

         grad_N = Tinv(x_xi(x_ele),N_xi)
      end function grad_N

end module inform
