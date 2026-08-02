
subroutine Sfun(lambda,Z,S)
use Defs
implicit none

!  Function s = Sfun(lambda,y). Shooting function value

!  .. Arguments ..
!  lambda    homotopic parameters
!  Z         Shooting unknown
!  S         Shooting value

real (kind=8), dimension(2), intent(in) :: lambda
real (kind=8), dimension(n), intent(in)  :: Z
real (kind=8), dimension(n), intent(out) :: S

!  .. Local declarations ..
integer :: i
real (kind=8), dimension(2*n) :: y(2*n), yd(2*n), lam(2)

!  .. Subroutines ..
external Pfun
external B2fun

!  .. First executable statement ..

lam(1) = lambda(1)
! force lambda(2) to be less than one
lam(2) = min(1.d0,lambda(2)) 

call Pfun(lam,Z,Y)

YD(:) = Y(:)/scale(:)

call B2fun(lam,YD,S)

return
end subroutine Sfun

!************************************************************
!************************************************************

subroutine Phifun(t,lambda,Y,PHI)
use Defs
implicit none

!  Function phi = phifun(l,lambda,y). Defines the modified (by variable change)
!  second member of the BVP

!  .. Arguments ..
!  t       time
!  lambda  homotopic parameters
!  Y       Modified State and Costate
!  PHI     Modified second member of the BVP

real (kind=8), intent(in) :: t
real (kind=8), dimension(2), intent(in) :: lambda
real (kind=8), dimension(2*n), intent(in) :: Y
real (kind=8), dimension(2*n), intent(out) :: PHI

!  .. Local declarations ..
integer :: i, j
real (kind=8), dimension(n) :: X, P
real (kind=8), dimension(m) :: U
real (kind=8), dimension(2*n) :: F

!  .. Subroutines ..
external Ufun
external Dhfun

!  .. First executable statement ..

! State and Costate descaling
X(:) = Y(1:n)/scale(1:n)
P(:) = Y(n+1:2*n)/scale(n+1:2*n)

! Force ellipse parameter to be stictly greater than 0
x(1) = max(1.d-3,x(1))

! Compute optimal control
call Ufun(t,lambda(2),X,P,U)

! Compute modified state and costate dynamics
call Dhfun(t,lambda(2),X,P,U,F)

PHI(:) = F(1:2*n)*scale(:)

return
end subroutine Phifun

!************************************************************
!************************************************************

subroutine Sfun_i(neq,Z,S,iflag)
use Defs
implicit none

!  Function [s iflag] = sfun_i(neq,z,iflag). Interface for SFUN
!  used by the NLE solver. See also SFUN.

!  .. Arguments ..
!  NEQ       Number of equations
!  Z         Shooting unknown
!  S         Shooting value 
!  IFLAG     Flag from the NLE solver

integer, intent(in) :: neq
real (kind=8), dimension(n), intent(in) :: Z
real (kind=8), dimension(n), intent(out) :: S
integer, intent(inout) :: iflag

!  .. Local declarations ..

!  .. Subroutines ..
external Sfun

!  .. First executable statement ..

call Sfun(l_com,Z,S)

return
end subroutine Sfun_i

!************************************************************
!************************************************************

subroutine Ssolve(lambda,ZI,Z,ifail)
use Defs
implicit none

!  Function z = ssolve(z_i). Solves the shooting equation
!  S(z) = 0 with initial guess z_i.

!  .. Arguments ..
!  lambda  homotopic parameters
!  ZI      Initial guess for z
!  Z       Shooting solution

real (kind=8), dimension(2), intent(in) :: lambda
real (kind=8), dimension(n), intent(in) :: ZI
real (kind=8), dimension(n), intent(out) :: Z
integer, intent (out) :: ifail

!  .. Local declarations ..
integer :: neq, maxfev, mode, i, nprint, nfev, lr
real (kind=8) :: xtol, factor, S(N), NLEW(lnlew)
integer :: ml, mu
real (kind=8) :: epsfcn

!  .. Subroutines ..
external Sfun_i
external hybrd
external sfun

!  .. First executable statement ..

Z = ZI ! initialization of Z

l_com = lambda
IFAIL = -1
NEQ = N
MAXFEV = 800*(NEQ+1) ! maximum number of evaluations
XTOL = SQRT(EPSTOL) ! NLE solver tolerance
ML = NEQ-1
MU = NEQ-1
EPSFCN = 0.0D0
MODE = 2 ! manual scaling
DO I = 1,NEQ
   NLEW(I) = 1.0D0
END DO
FACTOR = 1.0D+2
NPRINT = -1 ! verbose mode (each NPRINT iterations)
LR = (NEQ*(NEQ+1))/2

CALL HYBRD(SFUN_I,NEQ,Z,S,XTOL,MAXFEV,ML,MU,EPSFCN,NLEW(1),MODE,FACTOR, &
&          NPRINT,IFAIL,NFEV,NLEW(1+6*NEQ+LR),NEQ,NLEW(1+6*NEQ),LR,     &
&          NLEW(1+NEQ),NLEW(1+2*NEQ),NLEW(1+3*NEQ),NLEW(1+4*NEQ),NLEW(1+5*NEQ))

return
end subroutine Ssolve

!************************************************************
!************************************************************

subroutine Term(lambda,Z)
use Defs
implicit none

!  Function term(lambda,z). Writes the result into the files ./out.dat .
!  ./out.dat contains all the information required to reproduce the
!  computation plus an evaluation of T, Y, U (solution) at NIT steps.
!  Also create "next.dat" of the same type as "in.dat"

!  .. Arguments ..
!  lambda     homotopic parameters
!  Z          Shooting solution
real (kind=8), dimension(2), intent(in) :: lambda
real (kind=8), dimension(n), intent(in) :: Z

!  .. Local declarations ..
integer :: i, j, ifail, ODEIW(lodeiw)
integer, parameter :: fid = 10
real (kind=8) :: S(n), tin, tout, ODEW(lodew)
real (kind=8) :: step, Y(2*n), YD(2*n), Yt(2*n), U(m)
logical :: failed

!  .. Subroutines ..
real (kind=8) :: dnrm2
external dnrm2
external Sfun
external Phifun
external rkf45
external Ufun
external pfun

!  .. First executable statement ..

call Sfun(lambda,Z,S)

OPEN(FID,FILE='next.dat')
WRITE(FID,'(5(e24.16,1x))') (Z(I)/SCALE(FREE0(I)),I=1,N)
WRITE(FID,'(10(i3,1x))') (FREE0(I),I=1,N)
WRITE(FID,'(5(e24.16,1x))') (Y0(I)/SCALE(I),I=1,2*N)
WRITE(FID,'(e24.16)') t0
WRITE(FID,'(e24.16)') tmin
WRITE(FID,'(2(e24.16,1x))') LAMBDA
WRITE(FID,'(5(e24.16,1x))') (PAR(I),I=1,LPAR)
WRITE(FID,'(i3,1x)') (IPAR(I),I=1,LIPAR)
WRITE(FID,*) NIT
WRITE(FID,*) jac_step    
WRITE(FID,*) TRACE      
WRITE(FID,'(2(e24.16,1x))') SSPAR(1:2)
CLOSE(FID)

OPEN(FID,FILE='out.dat')
WRITE(FID,*) '# trajectory '
WRITE(FID,'(" # FREE0 = ",7(i3,1x))') (FREE0(I),I=1,N)
WRITE(FID,'(" # Y0 = ",5(e24.16,1x))') Y0/SCALE(:)
WRITE(FID,'(" # T0 = ",(e24.16))') t0
WRITE(FID,'(" # TF = ",(e24.16))') tmin
WRITE(FID,'(" # LAMBDA = ",2(e24.16,1x))') LAMBDA
WRITE(FID,'(" # PAR = ",5(e24.16,1x))') PAR
WRITE(FID,'(" # IPAR = ",5(i3,1x))') IPAR
WRITE(FID,*) " # NIT = ", nit
WRITE(FID,'(" # Z = ",5(e24.16e3,1x))') (Z(I)/SCALE(FREE0(I)),I=1,N)
WRITE(FID,'(" # S = ",5(e24.16e3,1x))') S
WRITE(FID,'(" # |S| = ",(e24.16))') DNRM2(N,S,1)
WRITE(FID,*) '# Data'

TIN = t0
STEP = (tmin-t0)*par(1)/NIT
Y = Y0
Y(free0(:)) = Z(1:n)
failed = .false.

! Descaling
DO I = 1,2*N
   YD(I) = Y(I)/SCALE(I)
END DO

CALL UFUN(TIN,lambda(2),YD(1),YD(N+1),U)

WRITE(FID,'(e24.16,2x,(e24.16e3))') TIN,(Y(I)/SCALE(I),I=1,2*N),(U(I),I=1,M)

DO J = 1,NIT
   IFAIL = 1 ! normal call
   TOUT = t0 + J*STEP
   CALL RKF45(lambda,PHIFUN,2*N,Y,TIN,TOUT,RELERR,ABSERR,IFAIL,ODEW,ODEIW)
   IF ((IFAIL.NE.2).and.(.not.failed)) THEN
      print *, 'Solution file may be incorrect (rkf45 error), ifail ', ifail
      failed = .true.
   END IF
   DO I = 1,2*N
      YD(I) = Y(I)/SCALE(I)
   END DO
   CALL UFUN(TOUT,lambda(2),YD(1),YD(N+1),U)

   WRITE(FID,'(e24.16,2x,(e24.16e3))') TOUT,(Y(I)/SCALE(I),I=1,2*N),(U(I),I=1,M)

END DO

CLOSE(FID)

print *, "Final mass (kg) = ", Y(n)/scale(n)

return
end subroutine Term

!************************************************************
!************************************************************

subroutine Init(lambda,ZI)
use Defs
implicit none

!  Function (lambda,zi) = init(). Initializes zi and the global variables
!  (see Defs.f90). Reads file ./in.dat with the following format:
!  ZI FREE0 Y0 L0 LMIN lambda PAR IPAR NIT JAC_STEP TRACE HMIN HMAX

!  .. Arguments ..
!  lambda  homotopic parameters
!  ZI      Initial guess for the shooting unknown

real (kind=8), dimension(2), intent(out) :: lambda
real (kind=8), dimension(n), intent(out) :: ZI

!  .. Local declarations ..
integer, parameter :: fid = 20
integer :: i, k

!  .. Subroutines ..
external rescale

!  .. First executable statement ..

OPEN(FID,FILE='in.dat')

READ(FID,*) (ZI(I),I=1,N)
READ(FID,*) (FREE0(I),I=1,N)
READ(FID,*) (Y0(I),I=1,2*N)
READ(FID,*) t0
READ(FID,*) tmin
READ(FID,*) lambda(:)
READ(FID,*) (PAR(I),I=1,LPAR)
READ(FID,*) (IPAR(I),I=1,LIPAR)
READ(FID,*) NIT
READ(FID,*) jac_step     ! finite differences step
READ(FID,*) TRACE        ! TRACE = 0
READ(FID,*) SSPAR(1:2)   ! min and max steplength for prediction

CLOSE(FID)

! Automatic scaling 
DO I = 1,2*N
   if (abs(Y0(i)) > scaltol) then
      k = int(log10(abs(Y0(i)))) + scalrange
      if (abs(Y0(i)) < 1.d0) k = k-1
      scale(i) = 10.d0**(-k)
   else
      scale(i) = 1.d0
   endif
ENDDO

ZI(:) = ZI(:)*scale(free0(:))
call rescale(ZI)
ZI(:) = ZI(:)/scale(free0(:))

! print some read data
WRITE(*,*) '[INIT] ZI = ', (ZI(I),I=1,N)
WRITE(*,*) '[INIT] FREE0 = ', (FREE0(I),I=1,N)
WRITE(*,*) '[INIT] Y0 = ', (Y0(I),I=1,2*N)
WRITE(*,*) '[INIT] T0 = ', t0
WRITE(*,*) '[INIT] TMIN = ', tmin
WRITE(*,*) '[INIT] PAR = ', (PAR(I),I=1,LPAR)
WRITE(*,*) '[INIT] IPAR = ', (IPAR(I),I=1,LIPAR)
WRITE(*,*) '[INIT] NIT = ',NIT
WRITE(*,*) '[INIT] FINITE DIFF STEP = ', jac_step  
WRITE(*,*) '[INIT] TRACE = ', TRACE

DO I = 1,N
   ZI(I) = SCALE(FREE0(I)) * ZI(I)
END DO
DO I = 1,2*N
   Y0(I) = SCALE(I) * Y0(I)
END DO

return
end subroutine Init


!************************************************************
!************************************************************

SUBROUTINE ELTIME(T)
IMPLICIT NONE
! Written on Mon Jul  2 11:52:16 MET DST 2001
! by Jean-Baptiste Caillau - ENSEEIHT-IRIT, UMR CNRS 5505

! Function t = eltime(). Returns the elapsed time.

! .. Global ..

! .. Arguments ..
! T        (output) REAL (KIND=R8)
!          Elapsed time
real (kind=8) :: T

! .. Local declarations ..
integer, dimension(8) :: aux1
real (kind=8) :: aux2

! .. Subroutines ..

! .. First executable statement ..
call date_and_time(VALUES=aux1)

aux2 = 1.d-3*aux1(8)+1.d0*aux1(7)+60.d0*aux1(6)+3600.d0*aux1(5)

T = aux2-T

RETURN

! .. Non executable statements ..
END SUBROUTINE ELTIME
