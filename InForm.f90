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
      N_xi = reshape([((((((vector                                    &
                     (G1(i,p)*N1(j,q)*N1(k,r),                        &
                      N1(i,p)*G1(j,q)*N1(k,r),                        &
                      N1(i,p)*N1(j,q)*G1(k,r)),                       &
                      i = 1,nqd1), j = 1,nqd1), k = 1,nqd1),          &
                      p = 1,nen1), q = 1,nen1), r = 1,nen1)],         &
                      [nqd,nen])

   contains

      pure function intgrt_1form(F_i,x_ele)
         implicit none
         type(vector),intent(in)  :: x_ele(nen)
         real(kind=rp)            :: intgrt_1form(nen)
         interface
            pure function F_i(x_elem)
               use control
               use linalg
               type(vector),intent(in) :: x_elem(nen)
               real(kind=rp)           :: F_i(nqd,nen)
            end function F_i
         end interface

         intgrt_1form = SUM(spread(wq,2,nen)*F_i(x_ele)               &
                           *spread(jac(x_ele),2,nen),1)
      end function intgrt_1form

      pure function intgrt_2form(F_ij,x_ele)
         implicit none
         type(vector),intent(in)  :: x_ele(nen)
         real(kind=rp)            :: intgrt_2form(nen)
         interface
            pure function F_ij(x_elem)
               use control
               use linalg
               type(vector),intent(in) :: x_elem(nen)
               real(kind=rp)           :: F_ij(nqd,nen,nen)
            end function F_ij
         end interface

         intgrt_2form = SUM(spread(wq,2,nen)*F_ij(x_ele)              &
                           *spread(jac(x_ele),2,nen),1)
      end function intgrt_2form

      pure function x_xi(x_ele)
         implicit none
         type(vector),intent(in)  :: x_ele(nen)
         type(tensor)             :: x_xi(nqd)
         integer                  :: i, j

         x_xi = tensor(vector(0.,0.,0.),                              &
                       vector(0.,0.,0.),                              &
                       vector(0.,0.,0.))
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
