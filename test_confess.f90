program tests
   use inform
   implicit none
   type(vector) :: x(nen)=[vector(-1.,-1.,-1.), &
                           vector( 1.,-1.,-1.), &
                           vector(-1., 1.,-1.), &
                           vector( 1., 1.,-1.), &
                           vector(-1.,-1., 1.), &
                           vector( 1.,-1., 1.), &
                           vector(-1., 1., 1.), &
                           vector( 1., 1., 1.)]

   write(*,*) intgrt_1form(vol,x)
end program tests
