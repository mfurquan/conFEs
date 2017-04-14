module deform
!*********************************************************************
! Definitionn of Form Functions
!*********************************************************************
   use inform
   implicit none

   contains

      pure function vol(x_elem)
         implicit none
         type(vector),intent(in) :: x_elem(nen)
         real(kind=rp)           :: vol(nqd,nen)

         vol = 1.0
      end function vol
end module deform
