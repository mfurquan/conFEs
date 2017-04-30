module linalg
!**********************************************************************
! Vector/Tensor types
!----------------------------------------------------------------------
! M. Furquan | Apr 14, 2017
!**********************************************************************
   use control
   implicit none

   type :: vector
      real(kind=rp) :: e(nsd)
   end type vector

   type :: tensor
      real(kind=rp) :: e(nsd,nsd)
   end type tensor

   interface operator (+)
      module procedure vector_add, tensor_add
   end interface operator (+)

   interface operator (*)
      module procedure scale_vector, scale_tensor
   end interface operator (*)

   interface operator (.otimes.)
      module procedure tensor_prod
   end interface operator (.otimes.)

   interface operator (.dot.)
      module procedure vec_dot_vec
   end interface operator (.dot.)

   contains

      elemental function det(t)
         implicit none
         type(tensor),intent(in)  :: t
         real(kind=rp)            :: det
         integer                  :: i

         det = t%e(1,1)*(t%e(2,2)*t%e(3,3) - t%e(3,2)*t%e(2,3))
             - t%e(1,2)*(t%e(2,1)*t%e(3,3) - t%e(3,1)*t%e(2,3))
             + t%e(1,3)*(t%e(2,1)*t%e(3,2) - t%e(3,1)*t%e(2,2))
      end function det

      elemental function detz(t)
         implicit none
         type(tensor),intent(in)  :: t
         type(vector),intent(in)  :: detz
         integer                  :: i

         detz = vector([t%e(1,2)*t%e(2,3) - t%e(2,2)*t%e(3,3),
                      - t%e(3,1)*t%e(1,3) - t%e(1,1)*t%e(3,3),
                        t%e(3,1)*t%e(1,2) - t%e(1,1)*t%e(3,2)])
      end function detz

      elemental function Tinv(A,b)
         implicit none
         type(tensor),intent(in)  :: A
         type(vector),intent(in)  :: b
         type(vector)             :: Tinv
         real(kind=rp)            :: k, A22, A23

         k         = A%e(2,1)/A%e(1,1)
         A22       = A%e(2,2) -  k*A%e(1,2)
         A23       = A%e(2,3) -  k*A%e(1,3)
         b2        = b%e(2)   -  k*b%e(1)

         k         = A%e(3,1)/A%e(1,1)
         A32       = A%e(3,2) -  k*A%e(1,2)
         A33       = A%e(3,3) -  k*A%e(1,3)
         b3        = b%e(3)   -  k*b%e(1)

         k         = A%e(3,2)/A%e(2,2)
         Tinv%e(3) = (b3 - k*b2)/(A33 -  k*A23)
         Tinv%e(2) = (b2 - Tinv%e(3)*A23)/A22
         Tinv%e(1) = (b%e(1) - Tinv%e(3)*A%e(1,3)                     &
                             - Tinv%e(2)*A%e(1,2))/A%e(1,1)
      end function Tinv

      elemental function vector_add(v1,v2)
         implicit none
         type(vector),intent(in)  :: v1, v2
         type(vector)             :: vector_add

         vector_add%e = v1%e + v2%e
      end function vector_add

      elemental function tensor_add(t1,t2)
         implicit none
         type(tensor),intent(in)  :: t1, t2
         type(tensor)             :: tensor_add

         tensor_add%e = t1%e + t2%e
      end function tensor_add

      elemental function scale_vector(k,v)
         implicit none
         real(kind=rp),intent(in) :: k
         type(vector),intent(in)  :: v
         type(vector)             :: scale_vector

         scale_vector%e = k*v%e
      end function scale_vector

      elemental function scale_tensor(k,t)
         implicit none
         real(kind=rp),intent(in) :: k
         type(tensor), intent(in) :: t
         type(tensor)             :: scale_tensor

         scale_tensor%e = k*t%e
      end function scale_tensor

      elemental function tensor_prod(v1,v2)
         implicit none
         type(vector),intent(in) :: v1, v2
         type(tensor)            :: tensor_prod

         tensor_prod%e = spread(v1%e,2,nsd)*spread(v2%e,1,nsd)
      end function tensor_prod

      elemental function vec_dot_vec(v1,v2)
         implicit none
         type(vector),intent(in) :: v1, v2
         real(kind=rp)           :: vec_dot_vec

         vec_dot_vec = SUM(v1%e*v1%e)
      end function vec_dot_vec
end module linalg
