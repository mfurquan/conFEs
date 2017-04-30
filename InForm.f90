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
      mask  = spread(reshape([((i/=j, i = 1,nen1), j = 1,nen1)],      &
                             [nen1,nen1]),1,nqd1)
   
   real(kind=rp),parameter,dimension(nqd) ::                          &
      wq    = [(((wg(qi+i)*wg(qi+j)*wg(qi+k),                         &
                i = 0,nqd1-1), j = 0,nqd1-1), k = 0,nqd1-1)]

   real(kind=rp),parameter,dimension(nen1) ::                         &
      xm    = [(-1+2*i/deg, i=0,deg)]                 

   real(kind=rp),parameter,dimension(nqd1,nen1,nen1) ::               &
      xq    = spread(spread(xg(qi:qf),2,nen1),3,nen1),                &
      xa    = spread(spread(xm,1,nen1),1,nqd1),                       &
      xb    = spread(spread(xm,1,nqd1),3,nen1)        

   real(kind=rp),parameter,dimension(nqd1,nen1) ::                    &
      N1    = product(xq - xa,3,mask)/product(xb - xa,3,mask),        &
      G1    = sum(1./(xq - xa),3,mask)*N1             

   real(kind=rp),parameter,dimension(nqdf,nenf) ::                    &
      S     = reshape([((((N1(i,p)*N1(j,q),                           &
                      i = 1,nqd1), j = 1,nqd1),                       &
                      p = 1,nen1), q = 1,nen1)],                      &
                      [nqdf,nenf])

   real(kind=rp),parameter,dimension(nqd,nen) ::                      &
      N     = reshape([((((((N1(i,p)*N1(j,q)*N1(k,r),                 &
                      i = 1,nqd1), j = 1,nqd1), k = 1,nqd1),          &
                      p = 1,nen1), q = 1,nen1), r = 1,nen1)],         &
                      [nqd,nen])

   type(vector),parameter,dimension(nqdf,nenf)  ::                    &
      S_xi  = reshape([((((vector([                                   &
                       G1(i,p)*N1(j,q),                               &
                       N1(i,p)*G1(j,q)],                              &
                       N1(i,p)*N1(j,q)*())      &
                       i = 1,nqd1), j = 1,nqd1),                      &
                       p = 1,nen1), q = 1,nen1)],                     &
                       [nqdf,nenf])

   type(vector),parameter,dimension(nqd,nen) ::                       &
      N_xi = reshape([((((((vector([                                  &
                      G1(i,p)*N1(j,q)*N1(k,r),                        &
                      N1(i,p)*G1(j,q)*N1(k,r),                        &
                      N1(i,p)*N1(j,q)*G1(k,r)]),                      &
                      i = 1,nqd1), j = 1,nqd1), k = 1,nqd1),          &
                      p = 1,nen1), q = 1,nen1), r = 1,nen1)],         &
                      [nqd,nen])


   contains

      pure function vol_int(expr,x_ele)
         implicit none
         real(kind=rp),intent(in) :: expr(nqd,ndf)
         type(vector), intent(in) :: x_ele(nen)
         real(kind=rp)            :: vol_int(ndf), jacob(nqd)
         integer                  :: idf

         jacob = jac(x_ele)
         do concurrent (idf = 1,ndf)
            vol_int(idf) =  SUM(wq*jacob*expr(:,idf))
         end do
      end function vol_int

      pure function surfdot_int(expr,x_fac)
         implicit none
         real(kind=rp),intent(in) :: expr(nqdf,ndf)
         type(vector), intent(in) :: x_fac(nenf)
         real(kind=rp)            :: surfdot_int
         type(vector)             :: jacob(nqdf)
         integer                  :: idf

         jacob = jacs(x_fac)
         do concurrent (idf = 1,ndf)
            surfdot_int(idf) = SUM(wq2*(jacob.dot.expr(:,idf)))
         end do
      end function surfdot_int

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

      pure function xs_xi(x_fac)
         implicit none
         type(vector),intent(in)  :: x_fac(nenf)
         type(tensor)             :: xs_xi(nqdf)
         integer                  :: i, j

         xs_xi = tensor(reshape([0.,0.,0.,                            &
                                 0.,0.,0.,                            &
                                 0.,0.,0.],[3,3])
         do concurrent (i = 1:nqdf, j = 1:nenf)
            xs_xi(i) = xs_xi(i) + S_xi(i,j).otimes.x_fac(j)
         end do
      end function xs_xi

      pure function jac(x_ele)
         implicit none
         type(vector),intent(in)  :: x_ele(nen)
         real(kind=rp)            :: jac(nqd)
         integer                  :: i

         jac = det(x_xi(x_ele))
      end function jac

      pure function jacs(x_fac)
         implicit none
         type(vector),intent(in)  :: x_fac(nenf)
         type(vector),intent(in)  :: jacs(nqdf)
         integer                  :: i

         jacs = detz(xs_xi(x_fac))
      end function jacf

      pure function grad_N(x_ele)
         implicit none
         type(vector),intent(in)  :: x_ele(nen)
         type(vector),intent(in)  :: grad_N(nqd,nen)

         grad_N = Tinv(x_xi(x_ele),N_xi)
      end function grad_N

      pure function grad_Nf(x_fac)
         implicit none
         type(vector),intent(in)  :: x_fac(nenf)
         type(vector),intent(in)  :: grad_Nf(nqdf,nenf)

         grad_Nf = Tinv(xf_xi(x_fac),Nf_xi)
      end function grad_S

end module inform
