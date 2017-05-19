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

   real(kind=rp),parameter,dimension(nqdf) ::                         &
      wqf   = [(((wg(qi+i)*wg(qi+j),                                  &
                i = 0,nqd1-1), j = 0,nqd1-1)]


   real(kind=rp),parameter,dimension(nen1) ::                         &
      xm    = [(-1+2*i/deg, i=0,deg)]                 

   real(kind=rp),parameter,dimension(nqd1,nen1,nen1) ::               &
      xq    = spread(spread(xg(qi:qf),2,nen1),3,nen1),                &
      xa    = spread(spread(xm,1,nen1),1,nqd1),                       &
      xb    = spread(spread(xm,1,nqd1),3,nen1)        


!   1 D  S H A P E  F U N C T I O N S  A N D   D E R I V A T I V E S
!----------------------------------------------------------------------
   real(kind=rp),parameter,dimension(nqd1,nen1) ::                    &
      N1    = product(xq - xa,3,mask)/product(xb - xa,3,mask),        &
      G1    = sum(1./(xq - xa),3,mask)*N1             

!  3 D   S H A P E   F U N C T I O N S   A N D   D E R I V A T I V E S
!----------------------------------------------------------------------
   real(kind=rp),parameter,dimension(nqd,nen) ::                      &
      N     = reshape([((((((N1(i,p)*N1(j,q)*N1(k,r),                 &
                      i = 1,nqd1), j = 1,nqd1), k = 1,nqd1),          &
                      p = 1,nen1), q = 1,nen1), r = 1,nen1)],         &
                      [nqd,nen])

   type(vector),parameter,dimension(nqd,nen) ::                       &
      N_xi  = reshape([((((((vector([                                 &
                       G1(i,p)*N1(j,q)*N1(k,r),                       &
                       N1(i,p)*G1(j,q)*N1(k,r),                       &
                       N1(i,p)*N1(j,q)*G1(k,r)]),                     &
                       i = 1,nqd1), j = 1,nqd1), k = 1,nqd1),         &
                       p = 1,nen1), q = 1,nen1), r = 1,nen1)],        &
                       [nqd,nen])

!  3 D  S H A P E   F U N C T I O N S   A N D   D E R I V A T I V E S  
!               E V A L U A T E D  O N   A   F A C E
!----------------------------------------------------------------------
   real(kind=rp),parameter,dimension(nqds,nen) ::                     &
      Nf    = reshape([((((((N1(i,p)*N1(j,q),                         &
                      i = 1,nqd1), j = 1,nqd1),                       &
                      p = 1,nen1), q = 1,nen1), r = 1,nen1)],         &
                      [nqds,nen])

   type(vector),parameter,dimension(nqds,nen) ::                      &
      Nf_xi = reshape([(((((vector([                                  &
                       G1(i,p)*N1(j,q),                               &
                       N1(i,p)*G1(j,q),                               &
                       N1(i,p)*N1(j,q)]),                             &
                       i = 1,nqd1), j = 1,nqd1),                      &
                       p = 1,nen1), q = 1,nen1), r = 1,nen1)],        &
                       [nqds,nen])

!  2 D  S H A P E   F U N C T I O N S  A N D   D E R I V A T I V E S 
!       F O R   F I N D I N G  S U R F A C E   J A C O B I A N
!----------------------------------------------------------------------
   real(kind=rp),parameter,dimension(nqds,nens) ::                    &
      S     = reshape([((((N1(i,p)*N1(j,q),                           &
                      i = 1,nqd1), j = 1,nqd1),                       &
                      p = 1,nen1), q = 1,nen1)],                      &
                      [nqds,nens])

   type(vector),parameter,dimension(nqds,nens)  ::                    &
      S_xi  = reshape([((((vector([                                   &
                       G1(i,p)*N1(j,q),                               &
                       N1(i,p)*G1(j,q),                               &
                       0.0_rp]),                                      &
                       i = 1,nqd1), j = 1,nqd1),                      &
                       p = 1,nen1), q = 1,nen1)],                     &
                       [nqds,nens])

!			          +-----------+
!			          |\           \
!	^ z	          | \     6     \
!	|		          |  \           \
!	|	    y	       |   +-----------+ nen^3->(nen1,nen1,nen1)
!	+------>        | 3 |           |
!	 \		           \  |           |
!	  \	            \ |     2     |
!	   \ x             \|           |
!	                    +-----------+ nen^2->(nen1,nen1,1)
   integer,parameter ::                                               &
      facmap(nens,6) = reshape(                                       &
                       [(i, i = 1,nen+1-nen1,nen1),                   &
                        (i, i = nen1,nen,nen1),                       &
                        ((i+j, i = 0,nens*(nsd-1,nens), j = 1,nen1)), &
                        ((i+j, i = ))],[nens,6])


   contains

!  V O L U M E   I N T E G R A L ,   S H A P E   F U N C T I O N
!          D E R I V A T I V E S  A N D  J A C O B I A N
!----------------------------------------------------------------------
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
         integer                  :: i

         jac = det(x_xi(x_ele))
      end function jac

      pure function grad_N(x_ele)
         implicit none
         type(vector),intent(in)  :: x_ele(nen)
         type(vector),intent(in)  :: grad_N(nqd,nen)

         grad_N = Tinv(x_xi(x_ele),N_xi)
      end function grad_N

!   S U R F A C E   I N T E G R A L ,   S H A P E   F U N C T I O N
!          D E R I V A T I V E S  A N D  J A C O B I A N
!----------------------------------------------------------------------
      pure function surfdot_int(expr,x_fac)
         implicit none
         real(kind=rp),intent(in) :: expr(nqds,nds)
         type(vector), intent(in) :: x_fac(nens)
         real(kind=rp)            :: surfdot_int
         type(vector)             :: jacob(nqds)
         integer                  :: idf

         jacob = jacs(x_fac)
         do concurrent (idf = 1,ndf)
            surfdot_int(idf) = SUM(wq2*(jacob.dot.expr(:,idf)))
         end do
      end function surfdot_int

      pure function xs_xi(x_fac)
         implicit none
         type(vector),intent(in)  :: x_fac(nens)
         type(tensor)             :: xs_xi(nqds)
         integer                  :: i, j

         xs_xi = tensor(reshape([0.,0.,0.,                            &
                                 0.,0.,0.,                            &
                                 0.,0.,0.],[3,3])
         do concurrent (i = 1:nqds, j = 1:nens)
            xs_xi(i) = xs_xi(i) + S_xi(i,j).otimes.x_fac(j)
         end do
      end function xs_xi

      pure function jacs(x_fac)
         implicit none
         type(vector),intent(in)  :: x_fac(nens)
         type(vector),intent(in)  :: jacs(nqds)
         integer                  :: i

         jacs = detz(xs_xi(x_fac))
      end function jacs

      pure function grad_Nf(x_fac)
         implicit none
         type(vector),intent(in)  :: x_fac(nens)
         type(vector),intent(in)  :: grad_Nf(nqds,nens)

         grad_Nf = Tinv(xf_xi(x_fac),Nf_xi)
      end function grad_S

end module inform
