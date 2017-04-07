module confess
   use gauss_legendre
   implicit none

!            A D J U S T A B L E   P A R A M E T E R S
!----------------------------------------------------------------------
   integer,parameter ::            &
      rp   = 8,                    & ! floating-point precision
      deg  = 2,                    & ! degree of interpolation
      nqd1 = 2                       ! no. of gauss-quadrature points

!     F I X E D  ( H A R D - C O D E D )  P A R A M E T E R S
!----------------------------------------------------------------------
   integer,parameter :: nsd = 3

!                 T Y P E  D E F I N I T I O N S
!----------------------------------------------------------------------
   type :: vector
      real(kind=rp) :: x, y, z
   end type vector

   type :: tensor
      type(vector) :: x, y, z
   end type tensor

   interface operator (.otimes.)
      module procedure tensor_prod
   end interface operator (.otimes.)

   interface operator (*)
      module procedure scalar_prod
   end interface operator (*)

!      S E L F - E V A L U A T I N G    P A R A M E T E R S
!----------------------------------------------------------------------
   integer,parameter ::                                               &
      qi   = (nqd1 - 1)*nqd1/2 + 1,                                   &
      qf   = (nqd1 + 1)*nqd1/2,                                       &
      nen1 = deg + 1,                                                 &
      nen  = nen1**nsd,                                               &
      nqd  = nqd1**nsd

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

      pure function x_xi(x_ele)
         implicit none
         type(vector),intent(in) :: x_ele(nen)
         type(tensor)            :: x_xi(nqd)
         integer                 :: i

         x_xi = SUM(N_xi.ocross.spread(x_ele,1,nqd),2)
      end function x_xi

      pure function jac(x_ele)
         implicit none
         type(vector),intent(in) :: x_ele(nen)
         real(kind=rp)           :: jac(nqd)
         type(tensor)            :: jm(nqd)
         integer                 :: i

         jac = det(x_xi(x_ele))
      end function jac

      pure function grad_N(x_ele)
         implicit none
         type(vector),intent(in) :: x_ele(nen)
         type(vector),intent(in) :: grad_N(nqd,nen)
         type(tensor)            :: tmp(nqd)

         grad_N = Tinv(x_xi(x_ele),N_xi)
      end function grad_N

      elemental function det(s)
         implicit none
         type(tensor),intent(in) :: s
         real(kind=rp)           :: det

         det = s%x%x*(s%y%y*s%z%z - s%z%y*s%y%z)                      &
             - s%x%y*(s%y%x*s%z%z - s%z%x*s%y%z)                      &
             + s%x%z*(s%y%x*s%z%y - s%z%x*s%y%y)                      &
      end function det

      elemental function Tinv(A,b)
         implicit none
         type(tensor),intent(in) :: A
         type(vector),intent(in) :: b
         type(vector)            :: Tinv
         real(kind=rp)           :: fact

         fact   = A%y%x/A%x%x
         A%y%y  = A%y%y - A%x%y*fact
         A%y%z  = A%y%z - A%x%z*fact
         b%y    = b%y   - b%x  *fact

         fact   = A%z%x/A%x%x
         A%z%y  = A%z%y - A%x%y*fact
         A%z%z  = A%z%z - A%x%z*fact
         b%z    = b%z   - b%x  *fact

         fact   = A%z%y/A%y%y
         Tinv%z = (b%z - b%y*fact)/(A%z%z - A%y%z*fact)
         Tinv%y = (b%y - A%y%z*Tinv%z)/A%y%y
         Tinv%x = (b%x - A%x%z*Tinv%z - A%x%y*Tinv%y)/A%x%x
      end function Tinv

      elemental function tensor_prod(v1,v2)
         implicit none
         type(vector),intent(in) :: v1, v2
         type(tensor)            :: tensor_prod

         tensor_prod%x = v1%x*v2
         tensor_prod%y = v1%y*v2
         tensor_prod%z = v1%z*v2
      end function tensor_prod

      elemental function scalar_prod(alpha,v)
         implicit none
         real(kind=rp),intent(in) :: alpha
         type(vector), intent(in) :: v
         type(vector)             :: scalar_prod

         scalar_prod%x = alpha*v%x
         scalar_prod%y = alpha*v%y
         scalar_prod%z = alpha*v%z
      end function scalar_prod

end module confess
