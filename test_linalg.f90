program test
   use linalg
   implicit none
   type(vector)  :: v1, v2
   type(tensor)  :: t1, t2
   real(kind=rp) :: k = 0.5

   v1%e = [1.0, 2.0, 3.0]
   v2   = vector([2.0, 5.0, 7.0])

   t1%e(1,:) = [3.0, 2.0, 6.0]
   t1%e(2,:) = [1.0, 2.0, 5.0]
   t1%e(3,:) = [5.0, 3.0, 2.0]


   write(*,*) 'v1+v2 =',v1+v2
   write(*,*) 't1+t2 =',t1+t2
   write(*,*) 'k.v2 =',k*v2
   write(*,*) 'k.t1 =',k*t1
end program test
