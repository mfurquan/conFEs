module assembler
   use control
   implicit none
   contains

      subroutine fill_vector(fn_ele,x_ele,u_ele,id_ele,V)
         implicit none
         real(kind=rp),intent(in) :: x_ele(nsd,nen), u_ele(ndf,nen)
         integer,intent(in)       :: id_ele(ndf,nen)
         real(kind=rp)            :: V(:), vloc(ndf)
         integer                  :: ien, idf, k
         interface
            pure function fn_ele(p,xe,ue)
               integer,intent(in)       :: p
               real(kind=rp),intent(in) :: xe(nsd,nen), ue(ndf,nen)
               real(kind=rp)            :: fn_ele(ndf)
            end function fn_ele
         end interface

         do concurrent (ien =  1:nen)
            vloc = fn_ele(ien,x_ele,u_ele)
            do concurrent (idf = 1:ndf, id_ele(idf,ien)>0)
               k =  id_ele(idf,ien)
               V(k) = V(k) + vloc(idf)
            end do
         end do
      end subroutine fill_vector

      subroutine fill_matrix(fn_ele,x_ele,u_ele,id_ele,M)
         implicit none
         real(kind=rp),intent(in) :: x_ele(nsd,nen), u_ele(ndf,nen)
         integer,intent(in)       :: id_ele(ndf,nen)
         real(kind=rp)            :: V(:,:), mloc(ndf)
         integer                  :: ien, jen, idf, k, l
         interface
            pure function fn_ele(p,q,xe,ue)
               integer,intent(in)       :: p, q
               real(kind=rp),intent(in) :: xe(nsd,nen), ue(ndf,nen)
               real(kind=rp)            :: fn_ele(ndf)
            end function fn_ele
         end interface

         do concurrent (ien =  1:nen, jen = 1:nen)
            mloc = fn_ele(ien,jen,x_ele,u_ele)
            do concurrent (idf = 1:ndf,                               &
               id_ele(idf,ien)>0, id_ele(idf,jen)>0)
               k =  id_ele(idf,ien); l = id_ele(idf,jen)
               V(k,l) = V(k,l) + mloc(idf)
            end do
         end do
      end subroutine fill_matrix

end module assembler
