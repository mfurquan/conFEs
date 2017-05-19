module control
!**********************************************************************
! Control parameters
! M. Furquan | Apr 07, 2017
!**********************************************************************
   implicit none

!            A D J U S T A B L E   P A R A M E T E R S
!----------------------------------------------------------------------
   integer,parameter ::            &
      rp   = 8,                    & ! floating-point precision
      deg  = 2,                    & ! degree of interpolation
      nqd1 = 2                       ! no. of gauss-quadrature points

!     F I X E D  ( H A R D - C O D E D )  P A R A M E T E R S
!----------------------------------------------------------------------
   integer,parameter :: nsd = 3, ndf = 1

!      S E L F - E V A L U A T I N G    P A R A M E T E R S
!----------------------------------------------------------------------
   integer,parameter ::                                               &
      qi   = (nqd1 - 1)*nqd1/2 + 1,                                   &
      qf   = (nqd1 + 1)*nqd1/2,                                       &
      nen1 = deg + 1,                                                 &
      nens = nen1**2,                                                 &
      nqds = nqd1**2,                                                 &
      nen  = nen1**3,                                                 &
      nqd  = nqd1**3

end module control
