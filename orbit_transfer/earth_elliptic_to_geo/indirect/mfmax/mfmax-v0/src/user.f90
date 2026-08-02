subroutine B2fun(lambda,Y,B2)
use Defs
implicit none

!  Function b2 = b2fun(lambda,y). Computes the boundary condition (at
!  final instant) of the modified BVP. **** Must be set by the user****.
!  See also DHFUN and UFUN.

!  .. Arguments ..
!  lambda  homotopic parameters
!  Y       Flow of BVP at final instant
!  B2      Boundary condition

real (kind=8), dimension(2), intent(in) :: lambda
real (kind=8), dimension(2*n), intent(in) :: Y
real (kind=8), dimension(n), intent(out) :: B2

!  .. Local declarations ..

!  .. First declarations ..
B2(1) = Y(1) - par(10)
B2(2) = Y(2) - par(11)
B2(3) = Y(3) - PAR(12)
B2(4) = Y(4) - PAR(13)
B2(5) = Y(5) - PAR(14)
B2(6) = Y(6) - (Y0(6)/scale(6)+par(3)*(par(4)-Y0(6)/scale(6)))
B2(7) = Y(7+N) ! pm = 0, final mass is free

return
end subroutine B2fun

!************************************************************
!************************************************************

subroutine Dhfun(t,lambda,X,P,U,F)
use Defs
implicit none

!  Computation of second member of state and costate dynamics
!  *** Must be set by user ***
!  See also PFUN, UFUN, B2FUN

!  .. Arguments ..
!  T        Time
!  lambda   homotopic parameters
!  X        State
!  P        Costate
!  U        Control
!  F        Second member of the dynamics

real (kind=8), intent(in) :: t
real (kind=8), intent(in) :: lambda
real (kind=8), dimension(n), intent(in) :: X
real (kind=8), dimension(n), intent(in) :: P
real (kind=8), dimension(m), intent(in) :: U
real (kind=8), dimension(2*n), intent(out) :: F

!  . Local declarations
real (kind=8) :: Tmax, Beta, mu0
real (kind=8) :: Pmu0, Pmu0Z, co, si
real (kind=8) :: Z, A, B, XX, H
real (kind=8) :: Zpoint, Apoint, Bpoint, XXpoint, Hpoint
real (kind=8) :: Hu3, p11, p22, p33, p44, p55
real (kind=8) :: ZpsZ, p6P3, p6PH3, Pmu03

!  .. Subroutines ..

!  .. First executable statement ..

Tmax = par(5)*tmaxconv
beta = par(6)
mu0 = par(7)

co = dcos(x(6))
si = dsin(x(6))
Z = 1.d0 + x(2)*co + x(3)*si
A = x(2) + (1.d0+Z)*co
B = x(3) + (1.d0+Z)*si
XX = 1.d0 + x(4)**2+x(5)**2
H = x(4)*si-x(5)*co
        
Pmu0 = dsqrt(x(1)/mu0)/x(7)

Hu3 = H*u(3)
Pmu0Z = Pmu0/Z*Tmax
Pmu03 = dsqrt(mu0/x(1)**3)*Z**2
p6P3 = p(6)*Pmu0Z*u(3)
p6PH3 = p6P3*H
      
F(1) = 2.d0*Pmu0Z*x(1)*u(2)
F(2) = Pmu0Z*(Z*si*u(1)+A*u(2)-x(3)*Hu3)
F(3) = Pmu0Z*(-Z*co*u(1)+B*u(2)+x(2)*Hu3)
F(4) = .5d0*Pmu0Z*XX*co*u(3)
F(5) = .5d0*Pmu0Z*XX*si*u(3)
F(6) = Pmu03+Pmu0Z*Hu3
F(7) = -Beta*Tmax*dsqrt(u(1)**2+u(2)**2+u(3)**2)
      
!  .. Dynamic of costate ..

p11 = p(1)*F(1)
p22 = p(2)*F(2)
p33 = p(3)*F(3)
p44 = p(4)*F(4)
p55 = p(5)*F(5)
   
F(8) = -0.5d0/x(1)*(3.d0*p11+p22+p33+p44+p55-3.d0*Pmu03*p(6)+ p6pH3)
      
Zpoint = co
Apoint = 1.d0 + co**2
Bpoint = co*si
      
ZpsZ = Zpoint/Z

F(9) = Pmu0Z*(-p(2)*((Apoint-A*ZpsZ)*u(2) + ZpsZ*x(3)*Hu3) - &
&      p(3)*((Bpoint-B*ZpsZ)*u(2) + (1.d0-x(2)*ZpsZ)*Hu3)) + &
&      ZpsZ*(p11+p44+p55+p(6)*(Pmu0Z*Hu3-2.d0*Pmu03)) 
    
Zpoint = si
Apoint = si*co
Bpoint = 1.d0+si**2
      
ZpsZ = Zpoint/Z

F(10) = -Pmu0Z*(p(2)*((Apoint-A*ZpsZ)*u(2)-(1.d0-x(3)*ZpsZ)*Hu3) + &
&       p(3)*((Bpoint-B*ZpsZ)*u(2) - ZpsZ*x(2)*Hu3)) + &
&       ZpsZ*(p11+p44+p55-p(6)*(2.d0*Pmu03-Pmu0Z*Hu3))   
      
XXpoint = 2.d0*X(4)
Hpoint = si

F(11) = Pmu0Z*Hpoint*u(3)*(p(2)*x(3) - p(3)*x(2) - p(6)) - &
&       XXpoint/XX*(p44+p55) 
      
XXpoint = 2.d0*X(5)
Hpoint = -co
      
F(12) = Pmu0Z*Hpoint*u(3)*(p(2)*x(3) - p(3)*x(2) - p(6))&
&       - XXpoint/XX*(p44+p55) 

Zpoint = -x(2)*si + x(3)*co
Apoint = -si*(1.d0+Z) + Zpoint*co
Bpoint = co*(1.d0+Z) + Zpoint*si
Hpoint = x(4)*co + x(5)*si
      
ZpsZ = Zpoint/Z

F(13) = Pmu0Z*(p(2)*(Z*co*u(1)+(Apoint-A*ZpsZ)*u(2)-x(3)*(Hpoint-H*ZpsZ)*u(3))+&
&       p(3)*(si*u(1)*Z+(Bpoint-B*ZpsZ)*u(2)+ x(2)*(Hpoint-H*ZpsZ)*u(3))) + &
&       0.5d0*Pmu0Z*XX*u(3)*(p(4)*(-si-co*ZpsZ)+p(5)*(co-si*ZpsZ)) + &
&       ZpsZ*(p(6)*2.d0*Pmu03 - p11) + p6P3*(Hpoint-H*ZpsZ)
F(13) = -F(13)
      
F(14) = (p11 + p22 + p33 + p44 + p55 + p6PH3)/x(7)

return
end subroutine Dhfun

!************************************************************
!************************************************************

subroutine Ufun(t,lambda,X,P,U)
use Defs
implicit none

!  Function u = u(t,lambda,x,p). Computes the control as a function of t,x
!  and p (maximum principle). **** Must be set by the user ****.
!  See also FFUN, B2FUN, PFUN.

!  .. Arguments ..
!  T        Time
!  lambda   homotopic parameters
!  X        State
!  P        Costate
!  U        Control
real (kind=8), intent(in) :: t
real (kind=8), intent(in) :: lambda
real (kind=8), dimension(n), intent(in) :: X
real (kind=8), dimension(n), intent(in) :: P
real (kind=8), dimension(m), intent(out) :: U

!  .. Local declarations ..
real (kind=8) :: Tmax, Beta, mu0
real (kind=8) :: Z, A, B, XX, H, co, si
real (kind=8) :: P3, Pmu0, Pmu0Z, comp(m), norm_comp, coeff, coeffc
!  .. Subroutines ..

!  .. First executable statement ..

Tmax = par(5)*tmaxconv
Beta = par(6)
mu0  = par(7)

co = dcos(x(6))
si = dsin(x(6))
Z = 1.d0 + x(2)*co + x(3)*si
A = x(2) + (1.d0+Z)*co
B = x(3) + (1.d0+Z)*si
XX = 1.d0 + x(4)**2+x(5)**2
H = x(4)*si-x(5)*co

Pmu0 = dsqrt(x(1)/mu0)/x(7)
Pmu0Z = Pmu0*Tmax/Z

comp(1) = Pmu0*Tmax*(p(2)*si - p(3)*co)
comp(2) = Pmu0Z*(2.d0*p(1)*x(1)+p(2)*A + p(3)*B)
comp(3) = Pmu0Z*(H*(-p(2)*x(3) + p(3)*x(2) + p(6)) + &
&         0.5d0*XX*(p(4)*co+p(5)*si))

norm_comp = dsqrt(comp(1)**2 + comp(2)**2 + comp(3)**2)

if (ipar(2) == 1) then  ! Convex homotopy
   coeff = lambda - Beta*Tmax*p(7)
   coeffc = 2.d0*(1.d0-lambda)
   
   if ((norm_comp <= coeff).or.(norm_comp == 0.d0)) then
      u(:) = 0.d0
   elseif (norm_comp >= (coeff+coeffc)) then
      u(:) = - comp(:)/norm_comp
   else
      u(:) = -comp(:)*(norm_comp-coeff)/(coeffc*norm_comp)
   end if
   
else ! power homotopy
   coeff = -Beta*Tmax*P(7)
   coeffc = 2.d0-lambda
   
   if (lambda < 1.d0) then
      if ((norm_comp <= coeff).or.(norm_comp == 0.d0)) then
         u(:) = 0.d0
      elseif (norm_comp >= (coeff+coeffc)) then
         u(:) = -comp(:)/norm_comp
      else
         u = -comp(:)*((norm_comp-coeff)/ &
&            (coeffc*norm_comp**(1.d0-lambda)))**(1.d0/(1.d0-lambda))
      end if
   else
      if ((norm_comp <= coeff+1.d0).or.(norm_comp == 0.d0)) then
         u(:) = 0.d0
      else
         u(:) = -comp(:)/norm_comp
      end if
   end if
   
end if

return
end subroutine Ufun

!************************************************************
!************************************************************

subroutine Pfun(lambda,Z,Y)
use Defs
implicit none

! Function y = pfun(lambda,z). Computes the evaluation at the final
! instant of the maximal flow associated with the modified BVP for the
! initial condition defined by z.

!  .. Arguments ..
!  lambda    homotopic parameters
!  Z         Shooting unknown
!  Y         Flow at final instant

real (kind=8), dimension(2), intent(in) :: lambda
real (kind=8), dimension(n), intent(in) :: Z
real (kind=8), dimension(2*n), intent(inout) :: Y

!  .. Local declarations ..
integer :: i, ifail, odeiw(lodeiw)
real (kind=8) :: tin, tout, odew(lodew)

!  .. Subroutines ..
external phifun
external rkf45

!  .. First executable statement ..
tin = t0
! Final longitude
tout = (tmin-t0)*par(1)+t0

y = y0
y(free0(:)) = z(:)

if (lambda(1) < 1.d0) then ! Initial condition homotopy
y(1:n-2) = (1.d0-lambda(1))*par(10:14)*scale(1:n-2)+lambda(1)*y0(1:n-2)
end if

ifail = 1   ! normal call
call rkf45(LAMBDA,PHIFUN,2*N,Y,TIN,TOUT,RELERR,ABSERR,&
     &          IFAIL,ODEW,ODEIW)
if (ifail.ne.2) then
   write(*,*) ' [PFUN] **** IFAIL (RKF45) =', ifail
endif

return
end subroutine Pfun

!*******************************************************
!*******************************************************

subroutine RHO(A,LAMBDA,X,V)

use Defs
implicit none
      
! Function V = RHO(A,lambda,X) compute the homotopy map

! ... Arguments ...
real (kind=8), intent(in):: A(:),X(:)
real (kind=8), intent(inout):: LAMBDA
real (kind=8), intent(out):: V(:)

! ... Local declaration ...
real (kind=8) :: LAM(2)

! ... Subroutine ...
external Sfun

! Force predicted point to have  LAMBDA .GE. 0  .
IF (lambda < 0.d0)  lambda = 0.d0

! Use global variable l_com to compute correct lambda(2)
if (l_com(1) == 1.d0) then ! main homotopy, lambda(1) = 1
   lam(1) = 1.d0
   lam(2) = lambda
else ! initial condition homotopy, lambda(1) < 1
   lam(1) = lambda
   lam(2) = 0.d0
endif

call Sfun(lam,X,V)

return
end subroutine RHO

!*******************************************************
!*******************************************************

SUBROUTINE JACFUN(A,LAMBDA,X,Q)

use Defs
implicit none

! Function Q = jacfun(A,lambda,X), compute the jacobian of the homotopy map
! Return in the matrix  Q the Jacobian matrix [D RHO/D LAMBDA, D RHO/DX]
! evaluated at the point (LAMBDA, X).

! ... Arguments ...
real (kind=8), intent(in)::A(:), X(:)
real (kind=8), intent(inout):: LAMBDA
real (kind=8), intent(out):: Q(:,:)

! ... Local declarations ...
real (kind=8), dimension(N) :: X_moins, X_plus ,S_moins, S_plus
real (kind=8) :: LAM(2), step
integer :: i

! ... First executable statements ...
      
! Force predicted point to have  LAMBDA .GE. 0
if (lambda < 0.d0) lambda = 0.d0

! Approximate Jacobian with finite center differences

if (lambda > 1.d-2) then
   step = jac_step*lambda
else
   step = jac_step
end if

if (lambda < 1.d0) then
   if (l_com(1) == 1.d0) then
      LAM(1) = 1.d0
      LAM(2) = lambda - step
      CALL SFUN(LAM,X,S_moins)
      LAM(2) = lambda + step
      CALL SFUN(LAM,X,S_plus)
      LAM(2) = lambda
   else ! we are in HomCI
      LAM(2) = 0.d0
      LAM(1) = lambda - step
      CALL SFUN(LAM,X,S_moins)
      LAM(1) = lambda + step
      CALL SFUN(LAM,X,S_plus)
      LAM(1) = lambda
   endif
   Q(:,1) = (S_plus(:) - S_moins(:))/(2.d0*step)
else ! dH/dlambda = 0 for lambda >= 1.
   Q(:,1) = 0.d0
end if

X_moins(:) = X(:)
X_plus(:) = X(:)

do i=1,n

   if (abs(X(i)) > 1.d-2) then
      step = jac_step*abs(X(i))
   else
      step = jac_step
   end if
   
   X_moins(i) = X(i)-step
   call sfun(LAM,X_moins,S_moins)
   X_moins(i) = X(i)
   
   X_plus(i) = X(i)+step
   call sfun(LAM,X_plus,S_plus)
   X_plus(i) = X(i)
   
   Q(:,i+1) = (S_plus(:) - S_moins(:))/(2.d0*step)

enddo

return
end subroutine Jacfun

