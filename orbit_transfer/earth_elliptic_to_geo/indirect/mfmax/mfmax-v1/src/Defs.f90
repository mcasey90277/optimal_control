
Module Defs

implicit none

! state dimension
integer, parameter :: n = 8
! control dimension
integer, parameter :: m = 3

! Real parameters
integer, parameter :: lpar = 11
real (kind=8), dimension(lpar), save :: par
! Integer parameter
integer, parameter :: lipar = 1  ! At least 1
integer, dimension(lipar), save :: ipar

! Factor to convert thrust from Newton to kg.Mm.h^-2
real (kind=8), parameter :: tmaxconv = 12.96d0

! Maximum number of prediction between two refinement (with single shooting)
real (kind=8), parameter :: limitd = 20

! Tolerance of scaling (do not scale variables with magnitude less than scaltol)
real (kind=8), parameter :: scaltol = 1.d-4

! Range of scaling (yscale in [10^(-scalrange);10^(1-scalerange)[)
integer, parameter :: scalrange = 1
! Multiplier factor for Pl scaling (if 1 then same range of scaling than others)
real (kind=8), parameter :: scalpL = 1.d-1

! Relative and Absolute error for shooting function integration
real (kind=8), parameter :: relerr = 1.e-10
real (kind=8), parameter :: abserr = 1.e-14

integer, dimension(n), save :: free0
real (kind=8), dimension(2*n), save  :: y0, scale
real (kind=8), save :: L0, Lmin
real (kind=8), dimension(2), save :: l_com
integer, save :: nit

real (kind=8), parameter :: epstol = 1.11022302462516D-16
integer, parameter :: lodeiw = 5, lodew = 3+6*2*N, lnlew = 2*N**2+7*N

integer, save :: trace

! Finite differences step
real (kind=8), save :: jac_step
! Required precision on curve tracking
real (kind=8), parameter :: arcre = 1.d-6
real (kind=8), parameter :: arcae = 1.d-6
! Required precision after zero path tracking
real (kind=8), parameter :: ansre = 1.e-10
real (kind=8), parameter :: ansae = 1.e-10
real (kind=8), dimension(4), save :: sspar = 0.d0

End module Defs
